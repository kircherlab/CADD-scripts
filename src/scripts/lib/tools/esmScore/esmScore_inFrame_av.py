"""
Description: input is vep annotated vcf and a file containing all peptide sequences with Ensemble transcript IDs. 
The script adds a score for frameshifts and stop gains to the info column of the vcf file. In brief, Scores for frameshift or stop gain variants were calculated 
for variants annotated with the Ensembl VEP tools' frameshift
or stop gain consequence annotation. Amino acid sequences of reference alleles were obtained as described above, with the only difference that a window of 250 amino acids was used. 
Log transformed probabilities were calculated and summed up for the entire reference amino acid sequence. 
Calculation of a score for the alternative allele was carried out based on the entire reference amino acid sequence. Here, we summed up log probabilities of the reference
sequences amino acids up to the point where the frameshift or stop gain occurred as obtained from the Ensembl VEP tools' annotations. For every amino acid that is lost due
to the frameshift or stop gain, we used the median of log transformed probabilities calculated from all possible amino acids at each individual position in the remaining sequence
and added them to the sum corresponding to the alternative allele.
The average of logs odds ratios between the reference and alternative sequences from the five models  was than used as a final score. 

Author: Thorben Maass, Max Schubach
Contact: tho.maass@uni-luebeck.de
Year:2023

Refractored by yangyxt (using numpy array instead of list appending, greatly improved the performance when dealing with huge VCF file)
"""

import warnings
import numpy as np
from Bio.bgzf import BgzfReader, BgzfWriter
import torch
from esm import pretrained
import click

# Constants
WINDOW_SIZE = 250
BATCH_SIZE = 20

def read_and_extract_vcf_data(input_file):
    """Reads the VCF file and extracts relevant information."""
    vcf_data = []
    with BgzfReader(input_file, "r") as vcf_file:
        for line in vcf_file:
            vcf_data.append(line)

    info_pos = {}
    for line in vcf_data:
        if line.startswith("##INFO="):
            info = line.split("|")
            for i, item in enumerate(info):
                if item in ("Feature", "Protein_position", "Amino_acids", "Consequence"):
                    info_pos[item] = i
            if len(info_pos) == 4:
                break

    # Preallocate NumPy arrays
    num_variants = sum(1 for line in vcf_data if not line.startswith("#"))
    variant_ids = np.empty(num_variants, dtype=object)
    transcript_ids = np.empty(num_variants, dtype=object)
    oAA = np.empty(num_variants, dtype=object)
    nAA = np.empty(num_variants, dtype=object)
    prot_pos_start = np.empty(num_variants, dtype=int)
    prot_pos_end = np.empty(num_variants, dtype=int)
    cons = np.empty(num_variants, dtype=object)

    idx = 0
    for variant in vcf_data:
        if not variant.startswith("#"):
            variant_entry = variant.split(",")
            for i in range(len(variant_entry)):
                variant_info = variant_entry[i].split("|")
                consequences = variant_info[info_pos["Consequence"]].split("&")
                if ("frameshift_variant" in consequences or "stop_gained" in consequences) and len(variant_info[info_pos["Amino_acids"]].split("/")) == 2:
                    variant_ids[idx] = variant_entry[0].split("|")[0]
                    transcript_ids[idx] = variant_info[info_pos["Feature"]]
                    cons[idx] = consequences
                    oAA[idx] = variant_info[info_pos["Amino_acids"]].split("/")[0]
                    nAA[idx] = variant_info[info_pos["Amino_acids"]].split("/")[1]
                    prot_pos_range = variant_info[info_pos["Protein_position"]].split("/")[0]
                    if "-" in prot_pos_range:
                        start, end = map(int, prot_pos_range.split("-"))
                        prot_pos_start[idx] = start
                        prot_pos_end[idx] = end
                    else:
                        pos = int(prot_pos_range)
                        prot_pos_start[idx] = pos
                        prot_pos_end[idx] = pos
                    idx += 1

    # Trim arrays to actual size
    variant_ids = variant_ids[:idx]
    transcript_ids = transcript_ids[:idx]
    cons = cons[:idx]
    oAA = oAA[:idx]
    nAA = nAA[:idx]
    prot_pos_start = prot_pos_start[:idx]
    prot_pos_end = prot_pos_end[:idx]

    return vcf_data, variant_ids, transcript_ids, oAA, nAA, prot_pos_start, prot_pos_end, cons

def process_transcript_data(transcript_file, transcript_ids, prot_pos_start, prot_pos_end):
    """Processes transcript data and creates aa_seq_ref."""
    with open(transcript_file, "r") as f:
        transcript_info_entries = f.read().split(">")[1:]

    transcript_info = []
    transcript_info_id = []
    for entry in transcript_info_entries:
        parts = entry.split(" ")
        transcript_info.append(parts)
        transcript_id_full = parts[4]
        transcript_id = transcript_id_full.split(".")[0]
        transcript_info_id.append(transcript_id)

    # Preallocate arrays
    num_transcripts = len(transcript_ids)
    aa_seq_ref = np.empty(num_transcripts, dtype=object)
    total_stop_codons = np.zeros(num_transcripts, dtype=int)
    stop_codons_before_mutation = np.zeros(num_transcripts, dtype=int)
    stop_codons_in_indel = np.zeros(num_transcripts, dtype=int)

    for j, transcript_id in enumerate(transcript_ids):
        transcript_found = False
        for i, info_id in enumerate(transcript_info_id):
            if info_id == transcript_id:
                transcript_found = True
                temp_seq = transcript_info[i][-1].replace("\n", "")

                stop_codons_before_mutation[j] = temp_seq[:prot_pos_start[j]].count("*")
                stop_codons_in_indel[j] = temp_seq[prot_pos_start[j]:prot_pos_end[j]].count("*")
                total_stop_codons[j] = temp_seq.count("*")

                aa_seq_ref[j] = temp_seq.replace("*", "")
                break

        if not transcript_found:
            aa_seq_ref[j] = "NA"
            stop_codons_before_mutation[j] = 9999
            total_stop_codons[j] = 9999
            stop_codons_in_indel[j] = 9999

    return aa_seq_ref, total_stop_codons, stop_codons_before_mutation, stop_codons_in_indel

def prepare_data_for_esm(aa_seq, transcript_ids, prot_pos_start, stop_codons_before_mutation):
    """Prepares data for the ESM model."""
    data = []
    prot_pos_mod = []
    for i in range(len(aa_seq)):
        if aa_seq[i] == "NA":
            continue

        adjusted_pos = prot_pos_start[i] - stop_codons_before_mutation[i]
        seq_len = len(aa_seq[i])

        if seq_len < WINDOW_SIZE:
            data.append((transcript_ids[i], aa_seq[i]))
            prot_pos_mod.append(adjusted_pos)
        elif adjusted_pos + 1 + WINDOW_SIZE // 2 <= seq_len and adjusted_pos + 1 - WINDOW_SIZE // 2 >= 1:
            start = adjusted_pos - WINDOW_SIZE // 2
            end = adjusted_pos + WINDOW_SIZE // 2
            data.append((transcript_ids[i], aa_seq[i][start:end]))
            prot_pos_mod.append(WINDOW_SIZE // 2)
        elif seq_len >= WINDOW_SIZE and adjusted_pos + 1 - WINDOW_SIZE // 2 < 1:
            data.append((transcript_ids[i], aa_seq[i][:WINDOW_SIZE]))
            prot_pos_mod.append(adjusted_pos)
        else:
            data.append((transcript_ids[i], aa_seq[i][-WINDOW_SIZE:]))
            prot_pos_mod.append(adjusted_pos - (seq_len - WINDOW_SIZE))

    return data, prot_pos_mod

def calculate_esm_scores(data, prot_pos_mod, modelsToUse, conseq, batch_size=BATCH_SIZE):
    """Runs the ESM model and calculates scores."""
    model_scores = []
    for model_name in modelsToUse:
        torch.cuda.empty_cache()
        model, alphabet = pretrained.load_model_and_alphabet(model_name)
        model.eval()
        batch_converter = alphabet.get_batch_converter()

        if torch.cuda.is_available():
            model = model.cuda()

        seq_scores = []
        for i in range(0, len(data), batch_size):
            batch_data = data[i:i + batch_size]
            batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)

            with torch.no_grad():
                if torch.cuda.is_available():
                    batch_tokens = batch_tokens.cuda()
                token_probs = torch.log_softmax(model(batch_tokens)["logits"], dim=-1).cpu()

            for j, (transcript_id, seq) in enumerate(batch_data):
                idx = i + j
                if conseq[idx] == "FS":
                    score = 0
                    for y, aa in enumerate(seq):
                        if y < prot_pos_mod[idx]:
                            score += token_probs[j, y + 1, alphabet.get_idx(aa)].item()
                        else:
                            aa_scores = [token_probs[j, y + 1, k].item() for k in range(4, 24)]
                            aa_scores.append(token_probs[j, y + 1, 26].item())
                            aa_scores.sort()
                            mid = len(aa_scores) // 2
                            median = (aa_scores[mid] + aa_scores[~mid]) / 2
                            score += median
                    seq_scores.append(score)
                elif conseq[idx] == "NA":
                    seq_scores.append(0)

        model_scores.append(seq_scores)

    return np.array(model_scores)

def annotate_vcf_and_write_output(vcf_data, variant_ids, transcript_ids, np_array_score_diff, modelsToUse, output_file, aa_seq_ref):
    """Adds scores to the VCF file and writes the output."""
    header_end = 0
    for i, line in enumerate(vcf_data):
        if line.startswith("#CHROM"):
            vcf_data[i - 1] += '##INFO=<ID=EsmScoreFrameshift,Number=.,Type=String,Description="esmScore for one submodels. Format: esmScore">\n'
            header_end = i
            break

    vcf_data_modified = vcf_data[:header_end + 1]

    for i in range(header_end + 1, len(vcf_data)):
        line = vcf_data[i]
        new_line = line
        j = 0
        while j < len(variant_ids):
            if line.split("|")[0] == variant_ids[j]:
                num_scores = 0
                for l in range(j, len(variant_ids)):
                    if line.split("|")[0] == variant_ids[l]:
                        num_scores += 1
                    else:
                        break

                new_line = new_line[:-1] + ";EsmScoreFrameshift" + "=" + new_line[-1:]

                for h in range(num_scores):
                    if aa_seq_ref[j + h] != "NA":
                        avg_score = np.mean(np_array_score_diff[:, j + h])
                        new_line = new_line[:-1] + "{0}|{1:.3f}".format(transcript_ids[j + h][11:], avg_score) + new_line[-1:]
                    else:
                        new_line = new_line[:-1] + "{0}|NA".format(transcript_ids[j + h][11:]) + new_line[-1:]

                    if h < num_scores - 1:
                        new_line = new_line[:-1] + "," + new_line[-1:]

                j += num_scores
            else:
                j += 1
        vcf_data_modified.append(new_line)

    with BgzfWriter(output_file, "w") as vcf_file_output:
        for line in vcf_data_modified:
            vcf_file_output.write(line)

@click.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    multiple=False,
    type=click.Path(exists=True, readable=True),
    help="bgzip compressed vcf file",
)
@click.option(
    "--transcripts",
    "transcript_file",
    required=True,
    multiple=False,
    type=click.Path(exists=True, readable=True),
    help="Fasta file of peptides for all transcripts, used for esm score",
)
@click.option(
    "--model-directory",
    "model_directory",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="EMS torch model directory, usedfor torch hub",
)
@click.option(
    "--model",
    "modelsToUse",
    multiple=True,
    type=str,
    default=[
        "esm1v_t33_650M_UR90S_1",
        "esm1v_t33_650M_UR90S_2",
        "esm1v_t33_650M_UR90S_3",
        "esm1v_t33_650M_UR90S_4",
        "esm1v_t33_650M_UR90S_5",
    ],
    help="Models for download, default is all 5 models",
)
@click.option(
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file, vcf file bgzip compresssed",
)
@click.option(
    "--batch-size",
    "batch_size",
    type=int,
    default=BATCH_SIZE,
    help="Batch size for esm model, default is 20",
)
def cli(
    input_file, transcript_file, model_directory, modelsToUse, output_file, batch_size
):
    """Main CLI function."""
    torch.hub.set_dir(model_directory)

    vcf_data, variant_ids, transcript_ids, oAA, nAA, prot_pos_start, prot_pos_end, cons = read_and_extract_vcf_data(input_file)
    aa_seq_ref, total_stop_codons, stop_codons_before_mutation, stop_codons_in_indel = process_transcript_data(
        transcript_file, transcript_ids, prot_pos_start, prot_pos_end
    )

    conseq = []
    aa_seq_alt = []
    for j in range(0, len(aa_seq_ref)):
        if aa_seq_ref[j] == "NA":
            aa_seq_alt.append("NA")
            conseq.append("NA")
        elif "*" in nAA[j] or "X" in nAA[j]:  # stop codon gained or complete frameshift
            aa_seq_alt.append(aa_seq_ref[j])  # add alt seq without stop codon
            conseq.append("FS")
        else:
            aa_seq_alt.append("NA")
            conseq.append("NA")
            warnings.warn(
                "there is a problem with the ensembl data base and vep. The ESMframesift score of this variant will be artificially set to 0. Affected transcript is "
                + str(transcript_ids[j])
            )

    data_ref, prot_pos_mod_ref = prepare_data_for_esm(aa_seq_ref, transcript_ids, prot_pos_start, stop_codons_before_mutation)
    data_alt, prot_pos_mod_alt = prepare_data_for_esm(aa_seq_alt, transcript_ids, prot_pos_start, stop_codons_before_mutation)

    ref_scores = calculate_esm_scores(data_ref, prot_pos_mod_ref, modelsToUse, conseq, batch_size)
    alt_scores = calculate_esm_scores(data_alt, prot_pos_mod_alt, modelsToUse, conseq, batch_size)
    np_array_score_diff = alt_scores - ref_scores

    annotate_vcf_and_write_output(vcf_data, variant_ids, transcript_ids, np_array_score_diff, modelsToUse, output_file, aa_seq_ref)

if __name__ == "__main__":
    cli()