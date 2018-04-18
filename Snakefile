# This snakemake pipeline needs to be run with a configfile (like config/CADD1.4-GRCh38.yml)

import glob

# generate output file names for all files in input folder
def inputfolder2outputfolder():

    vcf_files = glob.glob('%s/*.vcf' % config['InputFolder'])
    vcfgz_files = glob.glob('%s/*.vcf.gz' % config['InputFolder'])

    output_files = [vcf.split('/')[-1].replace('.vcf','.tsv.gz') for vcf in vcf_files]
    output_files.extend([vcf.split('/')[-1].replace('.vcf.gz','.tsv.gz') for vcf in vcfgz_files])

    return ['%s/%s' % (config['OutputFolder'], out) for out in output_files]

### get all files from input folder
rule default:
    input:
        inputfolder2outputfolder()

# in case an uncompressed vcf file is loaded
rule compress:
    input:
        '%s/{file}.vcf' % config['InputFolder']
    output:
        temp('%s/{file}.vcf.gz' % config['InputFolder'])
    shell:
        '''
        gzip -c {input} > {output};
        '''

rule prescore:
    input:
        inputfile='%s/{file}.vcf.gz' % config['InputFolder'],
        prescorefiles=glob.glob('%s*.tsv.gz' % config['PrescoredFolder'])
    output:
        prescored=temp('data/pipeline/prescored/{GenomeBuild}/{file}.tsv.gz'),
        novel=temp('data/pipeline/input/{GenomeBuild}/{file}.vcf')
    conda:
        'src/environment.yml'
    shell:
        '''
        zcat {input.inputfile} \
        | python {config[CADDpath]}/src/scripts/VCF2vepVCF.py \
        | sort -k1,1 -k2,2n -k3,3 -k4,4 | uniq > {output.novel};

        echo '# Prescored variant file' | gzip -c > {output.prescored};
        for PRESCORED in {input.prescorefiles}
        do
            cat {output.novel} \
            | python {config[CADDpath]}/src/scripts/extract_scored.py --header \
                -p $PRESCORED --found_out=>( gzip -c >> {output.prescored} ) \
            > {output.novel}.tmp;
            mv {output.novel}.tmp {output.novel};
        done;
        '''

rule annotate:
    input:
        'data/pipeline/input/{GenomeBuild}/{file}.vcf'
    output:
        temp('data/pipeline/annotation/{GenomeBuild}/{file}.tsv.gz')
    conda:
        'src/environment.yml'
    shell:
        '''
        export CADD={config[CADDpath]};
        cat {input} \
        | vep --quiet --cache --buffer 1000 --no_stats --offline --vcf \
            --dir data/vep/{config[GenomeBuild]} \
            --species homo_sapiens --db_version=90 \
            --assembly {config[GenomeBuild]} --regulatory --sift b \
            --polyphen b --per_gene --ccds --domains --numbers --canonical \
            --total_length --force_overwrite --format vcf --output_file STDOUT \
        | python {config[CADDpath]}/src/scripts/annotateVEPvcf.py -c {config[ReferenceConfig]} \
        | gzip -c > {output};
        '''

rule impute:
    input:
        'data/pipeline/annotation/{GenomeBuild}/{file}.tsv.gz'
    output:
        temp('data/pipeline/impute/{GenomeBuild}/{file}.csv.gz')
    conda:
        'src/environment.yml'
    shell:
        """
        zcat {input} \
        | python {config[CADDpath]}/src/scripts/trackTransformation.py -b \
            -c {config[ImputeConfig]} -o {output} --noheader;
        """

rule generate_sparse_matrix:
    input:
        'data/pipeline/impute/{GenomeBuild}/{file}.csv.gz'
    output:
        temp('data/pipeline/sparse/{GenomeBuild}/{file}.npz')
    conda:
        'src/environment.yml'
    shell:
        """
        python {config[CADDpath]}/src/scripts/saveSparseMatrix.py -i {input} -o {output};
        """

rule predict:
    input:
        matrix='data/pipeline/sparse/{GenomeBuild}/{file}.npz',
        anno='data/pipeline/annotation/{GenomeBuild}/{file}.tsv.gz',
        model='data/models/{GenomeBuild}/%s' % config['Model'],
        table='data/models/{GenomeBuild}/%s' % config['ConversionTable']
    output:
        temp('data/pipeline/result/{GenomeBuild}/{file}.tsv.gz')
    conda:
        'src/environment.yml'
    shell:
        """
        python {config[CADDpath]}/src/scripts/predictSKmodel.py \
            -i {input.matrix} -m {input.model} -a {input.anno} \
        | python {config[CADDpath]}/src/scripts/max_line_hierarchy.py --all \
        | python {config[CADDpath]}/src/scripts/appendPHREDscore.py -t {input.table} \
        | gzip -c > {output};
        """

rule join:
    input:
        prescored='data/pipeline/prescored/%s/{file}.tsv.gz' % config['GenomeBuild'],
        novel='data/pipeline/result/%s/{file}.tsv.gz' % config['GenomeBuild']
    output:
        '%s/{file}.tsv.gz' % config['OutputFolder']
    shell:
        '''
        (
            echo "##CADD {config[GenomeBuild]}-v1.4 (c) University of Washington, Hudson-Alpha Institute for Biotechnology and Berlin Institute of Health 2013-2018. All rights reserved.";
            zcat {input.novel} | head -n 1;
            zcat {input.prescored} {input.novel} | grep -v "^#" | sort -k1,1 -k2,2n -k3,3 -k4,4;
            ) | bgzip -c > {output};
        '''
