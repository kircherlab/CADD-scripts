#!/usr/bin/env python
# -*- coding: ASCII -*-

import sys
import os
import pysam
from optparse import OptionParser
import multiprocessing as mp

def buffer_vcf_by_chromosome(input_stream):
    """Read VCF from input stream and buffer by chromosome"""
    vcf_by_chrom = {}
    header_lines = []
    
    for line in input_stream:
        if line.startswith('#'):
            header_lines.append(line)
            continue
            
        fields = line.strip().split('\t')
        chrom = fields[0]
        
        if chrom not in vcf_by_chrom:
            vcf_by_chrom[chrom] = []
        vcf_by_chrom[chrom].append(line)
            
    return header_lines, vcf_by_chrom

def setup_output_dir(output_base, chrom):
    """Create chromosome-specific output directory"""
    chrom_dir = os.path.join(output_base, chrom)
    try:
        os.makedirs(chrom_dir)
    except OSError:
        if not os.path.isdir(chrom_dir):
            raise
    return chrom_dir

def extract_prescored_chromosome(input_file, output_base, chrom):
    """Extract records for a single chromosome from prescored TSV file"""
    try:
        # Setup output directory
        chrom_dir = setup_output_dir(output_base, chrom)
        input_file_name = os.path.basename(input_file)
        input_file_name_base = input_file_name.rsplit('.', 1)[0]
        output_file = os.path.join(chrom_dir, "{0}.{1}.tsv".format(input_file_name_base, chrom))
        compressed_file = "{0}.gz".format(output_file)
        
        # Check if extraction is needed
        if os.path.exists(compressed_file):
            if os.path.getmtime(compressed_file) > os.path.getmtime(input_file):
                return compressed_file
        
        # Extract records for this chromosome using tabix
        tbx = pysam.TabixFile(input_file)
        with open(output_file, 'w') as f:
            for row in tbx.fetch(chrom):
                f.write("{0}\n".format(row))
        
        # Compress and index the output file
        pysam.tabix_compress(output_file, compressed_file, force=True)
        pysam.tabix_index(compressed_file, 
                         preset=None,
                         force=True,
                         seq_col=0,
                         start_col=1,
                         end_col=1,
                         zerobased=False)
        
        # Remove uncompressed file
        os.remove(output_file)
        return compressed_file
        
    except Exception as e:
        raise Exception("Error extracting prescored chromosome {0}: {1}".format(chrom, str(e)))

def process_chromosome(args):
    """Process a single chromosome"""
    chrom, vcf_lines, prescored_file, temp_dir, fpos, fref, falt = args
    try:
        # First extract prescored records for this chromosome
        prescored_chrom_file = extract_prescored_chromosome(
            prescored_file,
            os.path.join(temp_dir, "prescored"),
            chrom
        )

        # Setup output files for this chromosome
        found_file = os.path.join(temp_dir, "matches", "found.{0}.tmp".format(chrom))
        notfound_file = os.path.join(temp_dir, "matches", "notfound.{0}.tmp".format(chrom))
        try:
            os.makedirs(os.path.dirname(found_file))
        except OSError:
            if not os.path.isdir(os.path.dirname(found_file)):
                raise

        if prescored_chrom_file is None:
            # Create empty found file and output all records to notfound file
            with open(notfound_file, 'w') as f_notfound:
                for line in vcf_lines:
                    f_notfound.write(line)
            return chrom, True
        
        # Open prescored tabix file
        pre_tbx = pysam.TabixFile(prescored_chrom_file)
        
        with open(found_file, 'w') as f_found, open(notfound_file, 'w') as f_notfound:
            # Process each variant
            for line in vcf_lines:
                fields = line.strip().split('\t')
                pos = int(fields[1])
                lref, allele = fields[-2], fields[-1].strip()
                found = False
                
                # Look for matches in prescored file
                for pre_line in pre_tbx.fetch(chrom, pos-1, pos):
                    vfields = pre_line.rstrip().split('\t')
                    if (vfields[fref] == lref) and (vfields[falt] == allele) and (vfields[fpos] == fields[1]):
                        f_found.write(pre_line + '\n')
                        found = True
                        break
                
                if not found:
                    f_notfound.write(line)
        
        return chrom, True
    except Exception as e:
        sys.stderr.write('Error processing chromosome {0}: {1}\n'.format(chrom, str(e)))
        return chrom, False


def main():
    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", help="Path to scored variants.")
    parser.add_option("-i", "--input", dest="input", help="Read variants from vcf file (default stdin)", default=None)
    parser.add_option("--found_out", dest="found_out", help="Write found variants to file (default: stdout)", default=None)
    parser.add_option("--header", dest="header", help="Write full header to output (default none)",
                      default=False, action="store_true")
    (options, args) = parser.parse_args()

    # Setup input stream
    input_stream = sys.stdin
    if options.input and options.input != "-":
        try:
            input_stream = open(options.input, 'r')
        except IOError as e:
            sys.stderr.write("Error opening input file: {0}\n".format(str(e)))
            sys.exit(1)

    # Setup output stream
    found_out = open(options.found_out, 'w') if options.found_out else sys.stdout
    
    try:
        # Create temporary directory
        temp_dir = "temp_extract_scored"
        try:
            os.makedirs(temp_dir)
        except OSError:
            if not os.path.isdir(temp_dir):
                raise

        # Initialize column indices
        fpos, fref, falt = 1, 2, 3
        
        # Check prescored file
        if not (os.path.exists(options.path) and os.path.exists(options.path+".tbi")):
            raise IOError("No valid file with pre-scored variants.\n")
        
        # Get header and column indices from prescored file
        pre_tbx = pysam.TabixFile(options.path, 'r')
        header = list(pre_tbx.header)
        
        # Write headers to output files if requested
        if options.header:
            for line in header:
                found_out.write(line+"\n")
        
        # Get column indices from header
        for line in header:
            try:
                fref = line.split('\t').index('Ref')
                falt = line.split('\t').index('Alt')
            except ValueError:
                pass

        # Buffer VCF data and get chromosomes
        header_lines, vcf_by_chrom = buffer_vcf_by_chromosome(input_stream)
        chromosomes = sorted(vcf_by_chrom.keys())
        
        # Write VCF headers to stdout
        for line in header_lines:
            sys.stdout.write(line)
        
        # Get number of threads from Snakemake
        threads = int(os.environ.get("SNAKEMAKE_THREADS", "1"))
        threads = min(threads, len(chromosomes))
        print("Using {0} threads to extract the scored variants across all chromosomes".format(threads), file=sys.stderr)
        
        # Setup parallel processing args
        process_args = [
            (chrom, vcf_by_chrom[chrom], options.path, temp_dir, fpos, fref, falt)
            for chrom in chromosomes
        ]
        
        # Process chromosomes in parallel
        pool = mp.Pool(threads)
        results = pool.map(process_chromosome, process_args)
        pool.close()
        pool.join()
        
        # Combine results
        for chrom, success in results:
            if success:
                found_file = os.path.join(temp_dir, "matches", "found.{0}.tmp".format(chrom))
                notfound_file = os.path.join(temp_dir, "matches", "notfound.{0}.tmp".format(chrom))
                
                if os.path.exists(found_file):
                    with open(found_file) as f:
                        for line in f:
                            found_out.write(line)
                    os.remove(found_file)
                
                if os.path.exists(notfound_file):
                    with open(notfound_file) as f:
                        for line in f:
                            sys.stdout.write(line)
                    os.remove(notfound_file)
        
        # Cleanup
        try:
            import shutil
            shutil.rmtree(temp_dir)
        except:
            pass

    finally:
        # Close input file if it's not stdin
        if options.input and options.input != "-":
            input_stream.close()
            
        # Close output file if it's not stdout
        if options.found_out:
            found_out.close()

if __name__ == "__main__":
    main()
