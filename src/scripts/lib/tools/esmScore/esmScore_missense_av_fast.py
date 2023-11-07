from esm_lib import *


"""
Description: input is vep annotated vcf and a file containing all peptide sequences with Ensemble transcript IDs. 
The script adds a score for frameshifts and stop gains to the info column of the vcf file. In brief, Scores for missense variants were calculated for variants annotated with the Ensembl
VEP tools' missense consequence annotation. Stop gain, stop lost, and stop retained consequences were explicitly excluded. Multiple amino acid substitutions were treated as inframe 
InDel . The amino acid sequence was obtained by matching the Ensembl transcript ID with entries in a FASTA file containing amino acid sequences of all transcripts in the 
Ensemble data base . Optional stop codons were removed from the sequence and the amino acid substitution was centered within a 350 amino acid window. As final score, we used the
averaged log odds ratio between the alternative
and reference amino acid employing the five different ESM-1v models "esm1v_t33_650M_UR90S_1", "esm1v_t33_650M_UR90S_2", "esm1v_t33_650M_UR90S_3", "esm1v_t33_650M_UR90S_4", and "esm1v_t33_650M_UR90S_5". 
Author: thorben Maass, Max Schubach
Contact: tho.maass@uni-luebeck.de
Year:2023
"""

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
def cli(input_file, transcript_file, model_directory, modelsToUse, output_file, batch_size):

    # load vcf data
    vcf_data, info_pos_Feature, info_pos_ProteinPosition, info_pos_AA, info_pos_consequence = read_variants(input_file)
    # load transcripts
    transcript_info, transcript_info_ids = read_transcripts(transcript_file)
    # get missense variants
    variant_ids, transcript_ids, oAA, nAA, proteinPositions = get_missense_variants_for_esm(
        vcf_data, info_pos_Feature, info_pos_ProteinPosition, info_pos_AA, info_pos_consequence)

    print(len(variant_ids), " missense variants found")
    # create the amino acid sequence
    data, proteinPositions_mod, aa_seq = create_aa_seq(transcript_ids,
                                                       transcript_info,
                                                       transcript_info_ids,
                                                       proteinPositions,
                                                       variant_ids,
                                                       oAA,)

    # starting computing scores
    torch.hub.set_dir(model_directory)
    modelScores = []  # scores of different models

    if len(data) >= 1:
        # load esm model(s)
        for k in range(0, len(modelsToUse), 1):
            print("load and run model" + modelsToUse[k])
            torch.cuda.empty_cache()
            model, alphabet = pretrained.load_model_and_alphabet(modelsToUse[k])
            model.eval()  # disables dropout for deterministic results
            batch_converter = alphabet.get_batch_converter()
            if torch.cuda.is_available():
                model = model.cuda()

            # apply es model to sequence, tokenProbs hat probs von allen aa an jeder pos basierend auf der seq in "data"

            print("run model " + modelsToUse[k])
            all_scores = []  # scores for all variants of a particular model
            for t in range(0, len(data), batch_size):
                if t + batch_size >= len(data):
                    batch_data = data[t:]
                else:
                    batch_data = data[t: t + batch_size]

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
                for i in range(0, len(batch_data), 1):
                    if batch_data[i][1] == "NA":
                        score = 0
                        all_scores.append(float(score))
                        continue
                    else:
                        score = (
                            token_probs[
                                i,
                                proteinPositions_mod[i + t] + 1 - 1,
                                alphabet.get_idx(nAA[i + t]),
                            ]
                            - token_probs[
                                i,
                                proteinPositions_mod[i + t] + 1 - 1,
                                alphabet.get_idx(oAA[i + t]),
                            ]
                        )  # protPos+1 weil Sequenz bei 0 los geht, +1 weil esm model das so will (hebt sich eigentlich auf, aber der vollstaendigkeit halber...)
                        all_scores.append(float(score))

            modelScores.append(all_scores)

    # write out variants
    vcf_data = variants_average_score(vcf_data, variant_ids, transcript_ids, aa_seq, len(modelsToUse), modelScores)
    write_variants(output_file, vcf_data)


if __name__ == "__main__":
    cli()
