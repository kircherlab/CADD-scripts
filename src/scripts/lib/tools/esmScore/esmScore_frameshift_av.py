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
"""

import warnings
import numpy as np
from Bio.bgzf import BgzfReader, BgzfWriter
import torch
from esm import pretrained
import click


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
    default=20,
    help="Batch size for esm model, default is 20",
)
def cli(
    input_file, transcript_file, model_directory, modelsToUse, output_file, batch_size
):
    torch.hub.set_dir(model_directory)

    # get information from vcf file with SNVs and write them into lists (erstmal Bsp, später automatisch aus info zeile extrahieren)
    vcf_file_data = BgzfReader(input_file, "r")  # TM_example.vcf.gz
    vcf_data = []
    for line in vcf_file_data:
        vcf_data.append(line)

    info_pos_Feature = False  # TranscriptID
    info_pos_ProteinPosition = False  # resdidue in protein that is mutated
    info_pos_AA = False  # mutation from aa (amino acid) x to y
    info_pos_consequence = False
    # identify positions of annotations importnat for esm score
    for line in vcf_data:
        if line[0:7] == "##INFO=":
            info = line.split("|")
            for i in range(0, len(info), 1):
                if info[i] == "Feature":
                    info_pos_Feature = i
                if info[i] == "Protein_position":
                    info_pos_ProteinPosition = i
                if info[i] == "Amino_acids":
                    info_pos_AA = i
                if info[i] == "Consequence":
                    info_pos_consequence = i
            break

    # extract annotations important for esm score, "NA" for non-coding variants
    variant_ids = []
    transcript_id = []
    oAA = []
    nAA = []
    protPosStart = []
    protPosEnd = []
    protPos_mod = []
    cons = []

    for variant in vcf_data:
        if variant[0:1] != "#":
            variant_entry = variant.split(",")
            for i in range(0, len(variant_entry), 1):
                variant_info = variant_entry[i].split("|")
                consequences = variant_info[info_pos_consequence].split("&")
                if (
                    "frameshift_variant" in consequences
                    or "stop_gained" in consequences
                ) and len(variant_info[info_pos_AA].split("/")) == 2:
                    variant_ids.append(variant_entry[0].split("|")[0])
                    transcript_id.append("transcript:" + variant_info[info_pos_Feature])
                    cons.append(variant_info[info_pos_consequence].split("&"))
                    oAA.append(
                        variant_info[info_pos_AA].split("/")[0]
                    )  # can also be "-" if there is an insertion
                    nAA.append(variant_info[info_pos_AA].split("/")[1])
                    if (
                        "-" in variant_info[info_pos_ProteinPosition].split("/")[0]
                    ):  # in case of frameshifts, vep only gives X as the new aa
                        protPosStart.append(
                            int(
                                variant_info[info_pos_ProteinPosition]
                                .split("/")[0]
                                .split("-")[0]
                            )
                        )
                        protPosEnd.append(
                            int(
                                variant_info[info_pos_ProteinPosition]
                                .split("/")[0]
                                .split("-")[1]
                            )
                        )
                    else:
                        protPosStart.append(
                            int(variant_info[info_pos_ProteinPosition].split("/")[0])
                        )
                        protPosEnd.append(
                            int(variant_info[info_pos_ProteinPosition].split("/")[0])
                        )
                    protPos_mod.append(False)

    # dissect file with all aa seqs to entries
    transcript_data = open(
        transcript_file, "r"
    )  # <Pfad zu "Homo_sapiens.GRCh38.pep.all.fa" >
    transcript_info_entries = transcript_data.read().split(
        ">"
    )  # evtl erstes > in file weglöschen
    transcript_data.close()
    transcript_info = []
    transcript_info_id = []

    # transcript info contains aa seqs, becomes processed later
    for i in range(0, len(transcript_info_entries), 1):
        if transcript_info_entries[i] != "":
            transcript_info.append(transcript_info_entries[i].split(" "))

    # transcript ids
    for i in range(0, len(transcript_info_entries), 1):
        if transcript_info_entries[i] != "":
            transcript_info_tmp = transcript_info_entries[i].split(" ")[4]
            pointAt = False
            # remove version of ENST ID vor comparison with vep annotation
            for p in range(0, len(transcript_info_tmp), 1):
                if transcript_info_tmp[p] == ".":
                    pointAt = p

            transcript_info_tmp = transcript_info_tmp[:pointAt]
            transcript_info_id.append(transcript_info_tmp)
    if (len(transcript_info_id)) != len(transcript_info):
        print("ERROR!!!!!!")

    # create list with aa_seq_refs of transcript_ids, mal gucken, ob man alle auf einmal uebergebenkann an esm model
    aa_seq_ref = []
    totalNumberOfStopCodons = []
    numberOfStopCodons = []
    numberOfStopCodonsInIndel = []
    for j in range(0, len(transcript_id), 1):
        transcript_found = False
        for i in range(
            1, len(transcript_info_id), 1
        ):  # start bei 1 statt 0 weil das inputfile mit ">" anfaengt und 0. element in aa_seq_ref einfach [] ist
            if transcript_info_id[i] == transcript_id[j]:
                transcript_found = True
                # prepare Seq remove remainings of header
                temp_seq = transcript_info[i][-1]
                for k in range(0, len(temp_seq), 1):
                    if temp_seq[k] != "\n":
                        k = k + 1
                    else:
                        k = k + 1
                        temp_seq = temp_seq[k:]
                        break

                # prepare seq (remove /n)
                forbidden_chars = "\n"
                for char in forbidden_chars:
                    temp_seq = temp_seq.replace(char, "")

                # count stop codons in seq before site of mutation
                numberOfStopCodons.append(0)
                if "*" in temp_seq:
                    for k in range(0, len(temp_seq), 1):
                        if temp_seq[k] == "*" and k < protPosStart[j]:
                            numberOfStopCodons[j] = numberOfStopCodons[j] + 1

                # count stop codons in Indel
                numberOfStopCodonsInIndel.append(0)
                if "*" in temp_seq:
                    for k in range(0, len(temp_seq), 1):
                        if (
                            temp_seq[k] == "*"
                            and k >= protPosStart[j]
                            and k < protPosEnd[j]
                        ):
                            numberOfStopCodonsInIndel[j] = (
                                numberOfStopCodonsInIndel[j] + 1
                            )

                # count stop codons in seq
                totalNumberOfStopCodons.append(0)
                if "*" in temp_seq:
                    for k in range(0, len(temp_seq), 1):
                        if temp_seq[k] == "*":
                            totalNumberOfStopCodons[j] = totalNumberOfStopCodons[j] + 1

                # remove additional stop codons (remove *)
                forbidden_chars = "*"
                for char in forbidden_chars:
                    temp_seq = temp_seq.replace(char, "")

                aa_seq_ref.append(temp_seq)
        if transcript_found == False:
            aa_seq_ref.append("NA")
            numberOfStopCodons.append(9999)
            totalNumberOfStopCodons.append(9999)
            numberOfStopCodonsInIndel.append(9999)

    conseq = []
    aa_seq_alt = []
    for j in range(0, len(aa_seq_ref), 1):
        # print(nAA[j])
        # print(oAA[j])
        # print("\n")

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
                + str(transcript_id[j])
            )

    # prepare data array for esm model

    window = 250
    data_ref = []
    for i in range(0, len(transcript_id), 1):
        if len(aa_seq_ref[i]) < window:
            data_ref.append((transcript_id[i], aa_seq_ref[i]))
            protPos_mod[i] = protPosStart[i] - numberOfStopCodons[i]

        elif (
            (len(aa_seq_ref[i]) >= window)
            and (
                protPosStart[i] - numberOfStopCodons[i] + 1 + window / 2
                <= len(aa_seq_ref[i])
            )
            and (protPosStart[i] - numberOfStopCodons[i] + 1 - window / 2 >= 1)
        ):
            data_ref.append(
                (
                    transcript_id[i],
                    aa_seq_ref[i][
                        protPosStart[i]
                        - numberOfStopCodons[i]
                        - int(window / 2) : protPosStart[i]
                        - numberOfStopCodons[i]
                        + int(window / 2)
                    ],
                )
            )  # esm model can only handle 1024 amino acids, so if the sequence is longer , just the sequece around the mutaion i
            protPos_mod[i] = int(
                len(
                    aa_seq_ref[i][
                        protPosStart[i]
                        - numberOfStopCodons[i]
                        - int(window / 2) : protPosStart[i]
                        - numberOfStopCodons[i]
                        + int(window / 2)
                    ]
                )
                / 2
            )

        elif (
            len(aa_seq_ref[i]) >= window
            and protPosStart[i] - numberOfStopCodons[i] + 1 - window / 2 < 1
        ):
            data_ref.append((transcript_id[i], aa_seq_ref[i][:window]))
            protPos_mod[i] = protPosStart[i] - numberOfStopCodons[i]

        else:
            data_ref.append((transcript_id[i], aa_seq_ref[i][-window:]))
            protPos_mod[i] = (
                protPosStart[i] - numberOfStopCodons[i] - (len(aa_seq_ref[i]) - window)
            )

    data_alt = []

    for i in range(0, len(transcript_id), 1):
        if len(aa_seq_alt[i]) < window:
            data_alt.append((transcript_id[i], aa_seq_alt[i]))

        elif (
            (len(aa_seq_alt[i]) >= window)
            and (
                protPosStart[i] - numberOfStopCodons[i] + 1 + window / 2
                <= len(aa_seq_alt[i])
            )
            and (protPosStart[i] - numberOfStopCodons[i] + 1 - window / 2 >= 1)
        ):
            data_alt.append(
                (
                    transcript_id[i],
                    aa_seq_alt[i][
                        protPosStart[i]
                        - numberOfStopCodons[i]
                        - int(window / 2) : protPosStart[i]
                        - numberOfStopCodons[i]
                        + int(window / 2)
                    ],
                )
            )  # esm model can only handle 1024 amino acids, so if the sequence is longer , just the sequece around the mutaion i

        elif (
            len(aa_seq_alt[i]) >= window
            and protPosStart[i] - numberOfStopCodons[i] + 1 - window / 2 < 1
        ):
            data_alt.append((transcript_id[i], aa_seq_alt[i][:window]))

        else:
            data_alt.append((transcript_id[i], aa_seq_alt[i][-window:]))

    ref_alt_scores = []
    # preloading reduces disk I/O time and bottlenecks
    preloaded_models = list()
    for k in range(0,len(modelsToUse),1):
        # print("Preloading model {}".format(k))
        model, alphabet = pretrained.load_model_and_alphabet(modelsToUse[k])
        preloaded_models.append([model,alphabet])
    # load esm model(s)
    for o in range(0, len([data_ref, data_alt]), 1):
        data = [data_ref, data_alt][o]
        modelScores = []  # scores of different models
        if len(data) >= 1:
            for k in range(0, len(modelsToUse), 1):
                torch.cuda.empty_cache()
                # print("fetch model {} from preloaded models".format(k))
                model, alphabet = preloaded_models[k] # pretrained.load_model_and_alphabet(modelsToUse[k])
                model.eval()  # disables dropout for deterministic results
                batch_converter = alphabet.get_batch_converter()

                if torch.cuda.is_available():
                    model = model.cuda()
                    # print("transferred to GPU")

                # apply es model to sequence, tokenProbs hat probs von allen aa an jeder pos basierend auf der seq in "data"
                seq_scores = []
                for t in range(0, len(data), batch_size):
                    # print("modelling allele {}/2 : model {}/{} : {}/{}".format(o+1,k+1,len(modelsToUse),t+1,len(data)))
                    if t + batch_size > len(data):
                        batch_data = data[t:]
                    else:
                        batch_data = data[t : t + batch_size]

                    batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
                    with torch.no_grad():  # setzt irgeineine flag auf false
                        if torch.cuda.is_available():
                            token_probs = torch.log_softmax(
                                model(batch_tokens.cuda())["logits"], dim=-1
                            )
                        else:
                            token_probs = torch.log_softmax(
                                model(batch_tokens)["logits"], dim=-1
                            )

                    # test and extract scores from tokenProbs
                    if o == 1:  # alt seqences
                        for i in range(0, len(batch_data), 1):
                            # print (str(t+i)+" of "+ str(len(data))+ "alt seqs")
                            if conseq[i + t] == "FS":
                                score = 0
                                for y in range(
                                    0, len(batch_data[i][1]), 1
                                ):  # iterating over single AA in sequence
                                    if y < protPos_mod[i + t]:
                                        score = (
                                            score
                                            + token_probs[
                                                i,
                                                y + 1,
                                                alphabet.get_idx(batch_data[i][1][y]),
                                            ]
                                        )
                                    else:
                                        # calc mean of all possible aa at this position
                                        aa_scores = []
                                        for l in range(4, 24, 1):
                                            aa_scores.append(
                                                token_probs[i, y + 1, l]
                                            )  # for all aa (except selenocystein)
                                        aa_scores.append(
                                            token_probs[i, y + 1, 26]
                                        )  # for selenocystein
                                        aa_scores.sort()
                                        mid = len(aa_scores) // 2
                                        median = (aa_scores[mid] + aa_scores[~mid]) / 2
                                        score = score + median

                                seq_scores.append(float(score))
                            elif conseq[i + t] == "NA":
                                score = 0
                                seq_scores.append(float(score))
                    elif o == 0:  # ref sequences
                        for i in range(0, len(batch_data), 1):
                            if conseq[i + t] == "FS":
                                score = 0
                                for y in range(
                                    0, len(batch_data[i][1]), 1
                                ):  # iterating over single AA in sequence
                                    score = (
                                        score
                                        + token_probs[
                                            i,
                                            y + 1,
                                            alphabet.get_idx(batch_data[i][1][y]),
                                        ]
                                    )
                                seq_scores.append(float(score))
                            elif conseq[i + t] == "NA":
                                score = 999  # sollte nacher rausgeschissen werden, kein score sollte -999 sein
                                seq_scores.append(float(score))

                modelScores.append(seq_scores)
        ref_alt_scores.append(modelScores)

    np_array_scores = np.array(ref_alt_scores)
    np_array_score_diff = np_array_scores[1] - np_array_scores[0]

    # write scores in cvf. file

    # get information from vcf file with SNVs and write them into lists (erstmal Bsp, später automatisch aus info zeile extrahieren)

    # identify positions of annotations important for esm score
    header_end = False
    for i in range(0, len(vcf_data), 1):
        if vcf_data[i][0:6] == "#CHROM":
            vcf_data[i - 1] = (
                vcf_data[i - 1]
                + "##INFO=<ID=EsmScoreFrameshift"
                + ',Number=.,Type=String,Description="esmScore for one submodels. Format: esmScore">\n'
            )
            header_end = i
            break

    for i in range(header_end + 1, len(vcf_data), 1):
        j = 0
        while j < len(variant_ids):
            if vcf_data[i].split("|")[0] == variant_ids[j]:
                # count number of vep entires per variant that result in an esm score (i.e. with consequence "missense")
                numberOfEsmScoresPerVariant = 0
                for l in range(j, len(variant_ids), 1):
                    if vcf_data[i].split("|")[0] == variant_ids[l]:
                        numberOfEsmScoresPerVariant = numberOfEsmScoresPerVariant + 1
                    else:
                        break

                # annotate vcf line with esm scores
                # for k in range (0, len(modelsToUse), 1):
                vcf_data[i] = (
                    vcf_data[i][:-1] + ";EsmScoreFrameshift" + "=" + vcf_data[i][-1:]
                )
                for h in range(0, numberOfEsmScoresPerVariant, 1):
                    if aa_seq_ref[j + h] != "NA":
                        average_score = 0
                        for k in range(0, len(modelsToUse), 1):
                            average_score = average_score + float(
                                np_array_score_diff[k][j + h]
                            )
                        average_score = average_score / len(modelsToUse)
                        vcf_data[i] = (
                            vcf_data[i][:-1]
                            + str(transcript_id[j + h][11:])
                            + "|"
                            + str(round(float(average_score), 3))
                            + vcf_data[i][-1:]
                        )
                    else:
                        vcf_data[i] = (
                            vcf_data[i][:-1]
                            + str(transcript_id[j + h][11:])
                            + "|"
                            + "NA"
                            + vcf_data[i][-1:]
                        )

                    if h != numberOfEsmScoresPerVariant - 1:
                        vcf_data[i] = vcf_data[i][:-1] + "," + vcf_data[i][-1:]

                j = j + numberOfEsmScoresPerVariant
            else:
                j = j + 1

    vcf_file_output = BgzfWriter(output_file, "w")
    for line in vcf_data:
        vcf_file_output.write(line)

    vcf_file_output.close()


if __name__ == "__main__":
    cli()
