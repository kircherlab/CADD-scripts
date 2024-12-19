"""
Description: input is vep annotated vcf and a file containing all peptide sequences with Ensemble transcript IDs. 
The script adds a score for frameshifts and stop gains to the info column of the vcf file. In brief, scores for inframe InDel variants were calculated for variants annotated with the Ensembl VEP tools'
inframe insertion or inframe deletion consequence annotation. Variants with missense consequence annotations were only used if multiple amino acids are substituted. 
Variants with stop gain, stop lost, and stop retained consequence annotations were explicitly excluded. Amino acid sequences of reference alleles were obtained as described 
above, with the only difference that a window of 250 amino acids was used. As for InDel variants multiple amino acids can be affected, the amino acid sequence corresponding to the 
alternative allele can differ in more than one position from the sequence of the reference allele. To account for this, we generated the entire amino acid sequence of the alternative
allele using the Ensembl VEP tools' annotations and the reference sequence. To calculate scores for inframe InDel variants, log transformed probabilities of the entire reference and 
alternative sequences were added up, respectively, and substracted from each other, yielding log odds ratios. 
The log odds ratios resulting from each of the five models were than averaged and used as final score. 
Author: Thorben Maass, Max Schubach
Contact: tho.maass@uni-luebeck.de
Year:2023

Refractored by yangyxt, replacing list appending with numpy array storing. Improving computation speed by 20x when dealing with large amount of variants.
"""

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
                if ("inframe_insertion" in consequences or "inframe_deletion" in consequences or "missense_variant" in consequences) and len(variant_info[info_pos["Amino_acids"]].split("/")) == 2:
                    variant_ids[idx] = variant_entry[0].split("|")[0]
                    transcript_ids[idx] = variant_info[info_pos["Feature"]].split(".")[0]
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

def prepare_data_for_esm(aa_seq, transcript_ids, prot_pos_start, stop_codons_before_mutation, window_size):
    """Prepares data for ESM, handling windowing and edge cases."""
    data = []
    prot_pos_mod = np.copy(prot_pos_start)

    for i, seq in enumerate(aa_seq):
        if seq == "NA":
            continue

        if len(seq) < window_size:
            data.append((transcript_ids[i], seq))
            prot_pos_mod[i] -= stop_codons_before_mutation[i]
        else:
            start = prot_pos_start[i] - stop_codons_before_mutation[i]
            if start + 1 + window_size // 2 <= len(seq) and start + 1 - window_size // 2 >= 1:
                seq_window = seq[start - window_size // 2 : start + window_size // 2]
                data.append((transcript_ids[i], seq_window))
                prot_pos_mod[i] = len(seq_window) // 2
            elif start + 1 - window_size // 2 < 1:
                data.append((transcript_ids[i], seq[:window_size]))
                prot_pos_mod[i] = start
            else:
                data.append((transcript_ids[i], seq[-window_size:]))
                prot_pos_mod[i] = start - (len(seq) - window_size)

    return data, prot_pos_mod

def calculate_esm_scores(data, prot_pos_mod, modelsToUse, conseq, batch_size):
    """Calculates ESM scores for a given dataset."""
    model_scores = []
    for model_name in modelsToUse:
        model, alphabet = pretrained.load_model_and_alphabet(model_name)
        batch_converter = alphabet.get_batch_converter()
        model.eval()
        if torch.cuda.is_available():
            model.cuda()

        seq_scores = []
        for i in range(0, len(data), batch_size):
            batch_end = min(i + batch_size, len(data))
            batch_data = data[i:batch_end]
            batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)

            with torch.no_grad():
                if torch.cuda.is_available():
                    batch_tokens = batch_tokens.cuda()
                token_probs = torch.log_softmax(model(batch_tokens)["logits"], dim=-1).cpu()

            for j, (_, seq) in enumerate(batch_data):
                score = 0
                if conseq[i + j] in ["inFrame", "MultiMissense"]:
                    for y, aa in enumerate(seq):
                        score += token_probs[j, y + 1, alphabet.get_idx(aa)].item()
                elif conseq[i + j] == "NA":
                    score = 999
                seq_scores.append(score)

        model_scores.append(seq_scores)

    return np.array(model_scores)

def annotate_vcf_and_write_output(vcf_data, variant_ids, transcript_ids, score_diff, modelsToUse, output_file, aa_seq_ref):
    """Annotates the VCF with ESM scores and writes the output."""
    header_end = 0
    for i, line in enumerate(vcf_data):
        if line.startswith("#CHROM"):
            vcf_data[i - 1] += "##INFO=<ID=EsmScoreInFrame,Number=.,Type=String,Description=\"esmScore for one submodels. Format: esmScore\">\n"
            header_end = i
            break

    vcf_output = BgzfWriter(output_file, "w")
    for line in vcf_data[:header_end + 1]:
        vcf_output.write(line)

    for i, line in enumerate(vcf_data[header_end + 1:]):
        j = 0
        while j < len(variant_ids):
            if line.split("|")[0] == variant_ids[j]:
                if aa_seq_ref[j] != "NA":
                    score_string = [str(round(score_diff[m][j], 3)) for m in range(len(modelsToUse))]
                    line = line.strip() + ";EsmScoreInFrame=" + ",".join(score_string) + "\n"
                else:
                    line = line.strip() + ";EsmScoreInFrame=NA\n"
                break
            j += 1
        vcf_output.write(line)
    vcf_output.close()

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
    aa_seq_alt = np.empty(len(aa_seq_ref), dtype=object)
    for j in range(0, len(aa_seq_ref), 1):
        if aa_seq_ref[j] == "NA":
            aa_seq_alt[j] = "NA"
            conseq.append("NA")
        elif len(nAA[j]) == len(oAA[j]) and "-" not in oAA[j] and "-" not in nAA[j]:
            nAA_mod = nAA[j].replace("*", "")
            aa_seq_alt[j] = (
                aa_seq_ref[j][: prot_pos_start[j] - stop_codons_before_mutation[j] - 1]
                + nAA_mod
                + aa_seq_ref[j][
                    prot_pos_end[j]
                    - stop_codons_before_mutation[j]
                    - stop_codons_in_indel[j] :
                ]
            )
            conseq.append("MultiMissense")
        elif len(nAA[j]) >= len(oAA[j]) and "-" in oAA[j]:
            aa_seq_alt[j] = (
                aa_seq_ref[j][: prot_pos_start[j] - stop_codons_before_mutation[j] - 1]
                + nAA[j]
                + aa_seq_ref[j][prot_pos_start[j] - stop_codons_before_mutation[j] - 1 :]
            )
            conseq.append("inFrame")
        elif len(nAA[j]) > len(oAA[j]):
            nAA_mod = nAA[j].replace("*", "")
            aa_seq_alt[j] = (
                aa_seq_ref[j][: prot_pos_start[j] - stop_codons_before_mutation[j] - 1]
                + nAA_mod
                + aa_seq_ref[j][
                    prot_pos_end[j]
                    - stop_codons_before_mutation[j]
                    - stop_codons_in_indel[j] :
                ]
            )
            conseq.append("inFrame")
        elif len(nAA[j]) <= len(oAA[j]) and "-" in nAA[j]:
            aa_seq_alt[j] = (
                aa_seq_ref[j][: prot_pos_start[j] - stop_codons_before_mutation[j] - 1]
                + aa_seq_ref[j][
                    prot_pos_end[j]
                    - stop_codons_in_indel[j]
                    - stop_codons_before_mutation[j] :
                ]
            )
            conseq.append("inFrame")
        elif len(nAA[j]) < len(oAA[j]):
            nAA_mod = nAA[j].replace("*", "")
            aa_seq_alt[j] = (
                aa_seq_ref[j][: prot_pos_start[j] - stop_codons_before_mutation[j] - 1]
                + nAA_mod
                + aa_seq_ref[j][
                    prot_pos_end[j]
                    - stop_codons_before_mutation[j]
                    - stop_codons_in_indel[j] :
                ]
            )
            conseq.append("inFrame")

    data_ref, prot_pos_mod_ref = prepare_data_for_esm(aa_seq_ref, transcript_ids, prot_pos_start, stop_codons_before_mutation, WINDOW_SIZE)
    data_alt, prot_pos_mod_alt = prepare_data_for_esm(aa_seq_alt, transcript_ids, prot_pos_start, stop_codons_before_mutation, WINDOW_SIZE)

    ref_scores = calculate_esm_scores(data_ref, prot_pos_mod_ref, modelsToUse, conseq, batch_size)
    alt_scores = calculate_esm_scores(data_alt, prot_pos_mod_alt, modelsToUse, conseq, batch_size)
    np_array_score_diff = ref_scores - alt_scores

    annotate_vcf_and_write_output(vcf_data, variant_ids, transcript_ids, np_array_score_diff, modelsToUse, output_file, aa_seq_ref)

if __name__ == "__main__":
    cli()