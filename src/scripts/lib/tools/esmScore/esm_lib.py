from Bio.bgzf import BgzfReader, BgzfWriter
import warnings


def read_variants(vcf_file):
    vcf_file_data = BgzfReader(vcf_file, "r")  # TM_example.vcf.gz
    vcf_data = []
    print("Read VCF data")
    for line in vcf_file_data:
        vcf_data.append(line)

    info_pos_Feature = False  # TranscriptID
    info_pos_ProteinPosition = False  # resdidue in protein that is mutated
    info_pos_AA = False  # mutation from aa (amino acid) x to y
    info_pos_consequence = False  # consequence of mutation

    print("identify positions of annotations importnat for esm score")
    for line in vcf_data:
        if line[0:14] == "##INFO=<ID=CSQ":
            info = line.split("|")
            for i in range(0, len(info), 1):
                if info[i] == "Feature":
                    info_pos_Feature = i
                elif info[i] == "Protein_position":
                    info_pos_ProteinPosition = i
                elif info[i] == "Amino_acids":
                    info_pos_AA = i
                elif info[i] == "Consequence":
                    info_pos_consequence = i
            break
    return (vcf_data, info_pos_Feature, info_pos_ProteinPosition, info_pos_AA, info_pos_consequence)


def get_vatiant_id(variant_line):
    return "-".join(variant_line.split("\t")[0:7])


def get_missense_variants_for_esm(vcf_data, info_pos_Feature, info_pos_ProteinPosition, info_pos_AA, info_pos_consequence):
    # extract annotations important for esm score, "NA" for non-coding variants
    print("extract annotations important for esm score, NA for non-coding variants")
    variant_ids = []
    transcript_ids = []
    oAA = []
    nAA = []
    proteinPositions = []

    for variant in vcf_data:
        if variant[0:1] != "#":

            variant_id = get_vatiant_id(variant)
            variant_entry = []
            esmMissense_entry = []
            info_fields = variant.split("\t")[7].split(";")

            for info_field in info_fields:
                if info_field[0:3] == "CSQ":
                    variant_entry = info_field[4:].split(",")
                elif info_field[0:16] == "EsmScoreMissense":
                    esmMissense_entry = info_field[17:].split(",")
            for i in range(0, len(variant_entry), 1):
                variant_info = variant_entry[i].split("|")
                consequences = variant_info[info_pos_consequence].split("&")
                if (
                    "missense_variant" in consequences
                    and len(variant_info[info_pos_AA].split("/")) == 2
                    and "stop_gained" not in consequences
                    and "stop_lost" not in consequences
                    and "stop_retained_variant" not in consequences
                    and "-" not in variant_info[info_pos_ProteinPosition].split("/")[0]
                ):
                    # check if ESM missense score is already annotated
                    is_already_annotated = False
                    for missense_info in esmMissense_entry:
                        if variant_info[info_pos_Feature] in missense_info.split("|"):
                            is_already_annotated = True
                            break
                    if is_already_annotated:
                        continue

                    # add variant info to annotation list
                    variant_ids.append(variant_id)
                    transcript_ids.append(
                        "transcript:" + variant_info[info_pos_Feature]
                    )
                    oAA.append(variant_info[info_pos_AA][0])
                    nAA.append(variant_info[info_pos_AA][2])
                    proteinPositions.append(
                        int(
                            variant_info[info_pos_ProteinPosition].split(
                                sep="/", maxsplit=1
                            )[0]
                        )
                    )
    return (variant_ids, transcript_ids, oAA, nAA, proteinPositions)


def read_transcripts(transcript_file):
    # dissect file with all aa seqs to entries
    print("read transcript file")
    transcript_data = open(
        transcript_file, "r"
    )  # <Pfad zu "Homo_sapiens.GRCh38.pep.all.fa" >
    transcript_info_entries = transcript_data.read().split(
        ">"
    )  # evtl erstes > in file weglöschen
    transcript_data.close()
    transcript_info = []
    transcript_info_ids = []

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
            transcript_info_ids.append(transcript_info_tmp)

    if (len(transcript_info_ids)) != len(transcript_info):
        raise ValueError("Error0: some transcripts may be missing in transcripts file ")

    return (transcript_info, transcript_info_ids)


def create_aa_seq(
    transcript_ids,
    transcript_info,
    transcript_info_ids,
    proteinPositions,
    variant_ids,
    oAA,
):
    # create list with aa_seqs of transcript_ids, mal gucken, ob man alle auf einmal uebergebenkann an esm model
    aa_seq = []
    numberOfStopCodons = (
        []
    )  # number of stop codons in fasta aa sequence before the site of mutation
    totalNumberOfStopCodons = []  # number of stop codons in fasta aa sequence
    for j in range(0, len(transcript_ids), 1):
        transcript_found = False
        for i in range(
            1, len(transcript_info), 1
        ):  # start bei 1 statt 0 weil das inputfile mit ">" anfaengt und 0. element in aa_seq einfach [] ist
            if (
                transcript_info_ids[i] == transcript_ids[j]
            ):  # -2 damit ".9" usw wegfaellt
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
                        if temp_seq[k] == "*" and k < proteinPositions[j]:
                            numberOfStopCodons[j] = numberOfStopCodons[j] + 1

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

                aa_seq.append(temp_seq)

        if transcript_found == False:
            aa_seq.append("NA")
            numberOfStopCodons.append("NA")
            totalNumberOfStopCodons.append("NA")

    # prepare data array for esm model, Problem: only give the coding sequences i
    print("prepare data array for esm model")
    window = 350
    data = []
    proteinPositions_mod = [False] * len(proteinPositions)  # if a protein is longer than 1024 aa
    for i in range(0, len(transcript_ids), 1):
        if aa_seq[i] == "NA":
            data.append((transcript_ids[i], aa_seq[i]))
            proteinPositions_mod[i] = "NA"

        elif len(aa_seq[i]) < window:
            data.append((transcript_ids[i], aa_seq[i]))
            proteinPositions_mod[i] = proteinPositions[i] - numberOfStopCodons[i]
            if (
                aa_seq[i][proteinPositions[i] - 1 - numberOfStopCodons[i]]
                != data[i][1][proteinPositions_mod[i] - 1]
            ):
                raise ValueError(
                    "Error1: amino acid sequence could not be processed for esm models. Error at entry:  "
                    + str(variant_ids[i])
                )
            if aa_seq[i][proteinPositions[i] - 1 - numberOfStopCodons[i]] != oAA[i]:
                # print(oAA[i])
                # print(aa_seq[i][protPos[i] - 1 - numberOfStopCodons[i]])
                # print(variant_ids[i])
                warnings.warn(
                    "Error1.1: ref amino acid according to vep does not match transcript data base. ESM scores for this variant are not correct. Error at entry:  "
                    + str(variant_ids[i])
                )

        elif (
            (len(aa_seq[i]) >= window)
            and (
                proteinPositions[i] - numberOfStopCodons[i] + 1 + window / 2
                <= len(aa_seq[i])
            )
            and (proteinPositions[i] - numberOfStopCodons[i] + 1 - window / 2 >= 1)
        ):
            data.append(
                (
                    transcript_ids[i],
                    aa_seq[i][
                        proteinPositions[i]
                        - numberOfStopCodons[i]
                        - int(window / 2): proteinPositions[i]
                        - numberOfStopCodons[i]
                        + int(window / 2)
                    ],
                )
            )  # esm model can only handle 1024 amino acids, so if the sequence is longer , just the sequece around the mutaion i
            proteinPositions_mod[i] = int(
                len(
                    aa_seq[i][
                        proteinPositions[i]
                        - numberOfStopCodons[i]
                        - int(window / 2): proteinPositions[i]
                        - numberOfStopCodons[i]
                        + int(window / 2)
                    ]
                )
                / 2
            )

            if (
                aa_seq[i][proteinPositions[i] - 1 - numberOfStopCodons[i]]
                != data[i][1][proteinPositions_mod[i] - 1]
            ):
                # print(aa_seq[i][protPos[i]-1-numberOfStopCodons[i]-2:protPos[i]-1-numberOfStopCodons[i]+2])
                # print(data[i][1][protPos_mod[i]-1-2:protPos_mod[i]-1+2])
                # print(numberOfStopCodons[i])
                raise ValueError(
                    "Error2: amino acid sequence could not be processed for esm models. Error at entry:  "
                    + str(variant_ids[i])
                )
            if aa_seq[i][proteinPositions[i] - 1 - numberOfStopCodons[i]] != oAA[i]:
                warnings.warn(
                    "Error1.1: ref amino acid according to vep does not match transcript data base. ESM scores for this variant are not correct. Error at entry:  "
                    + str(variant_ids[i])
                )

        elif (
            len(aa_seq[i]) >= window
            and proteinPositions[i] - numberOfStopCodons[i] + 1 - window / 2 < 1
        ):
            data.append((transcript_ids[i], aa_seq[i][:window]))
            proteinPositions_mod[i] = proteinPositions[i] - numberOfStopCodons[i]

            if (
                aa_seq[i][proteinPositions[i] - 1 - numberOfStopCodons[i]]
                != data[i][1][proteinPositions_mod[i] - 1]
            ):
                raise ValueError(
                    "Error3: amino acid sequence could not be processed for esm models. Error at entry:  "
                    + str(variant_ids[i])
                )
            if aa_seq[i][proteinPositions[i] - 1 - numberOfStopCodons[i]] != oAA[i]:
                warnings.warn(
                    "Error1.1: ref amino acid according to vep does not match transcript data base. ESM scores for this variant are not correct. Error at entry:  "
                    + str(variant_ids[i])
                )

        else:
            data.append((transcript_ids[i], aa_seq[i][-window:]))
            proteinPositions_mod[i] = (
                proteinPositions[i] - (len(aa_seq[i]) - window) - numberOfStopCodons[i]
            )

            if (
                aa_seq[i][proteinPositions[i] - 1 - numberOfStopCodons[i]]
                != data[i][1][proteinPositions_mod[i] - 1]
            ):
                raise ValueError(
                    "Error4: amino acid sequence could not be processed for esm models. Error at entry:  "
                    + str(variant_ids[i])
                )
            if aa_seq[i][proteinPositions[i] - 1 - numberOfStopCodons[i]] != oAA[i]:
                warnings.warn(
                    "Error1.1: ref amino acid according to vep does not match transcript data base. ESM scores for this variant are not correct. Error at entry:  "
                    + str(variant_ids[i])
                )
    return (data, proteinPositions_mod, aa_seq)



def modify_header(vcf_data):
        # identify positions of annotations important for esm score
    header_end = False
    header_contain_missense = False
    for i in range(0, len(vcf_data), 1):
        if vcf_data[i][0:27] == "##INFO=<ID=EsmScoreMissense":
            header_contain_missense = True
        elif vcf_data[i][0:6] == "#CHROM":
            header_end = i - 1
            if not header_contain_missense:
                vcf_data[i - 1] = (
                    vcf_data[i - 1]
                    + "##INFO=<ID=EsmScoreMissense"
                    + ',Number=.,Type=String,Description="esmScore av of 5 submodels. Format: SYMBOL|esmScore">\n'
                )
                header_end += 1
            break
    return(vcf_data, header_end)

def variants_average_score(
    vcf_data, variant_ids, transcript_ids, aa_seq, numModels, modelScores
):
    # identify positions of annotations important for esm score
    vcf_data, header_end = modify_header(vcf_data)

    for i in range(header_end + 1, len(vcf_data), 1):
        j = 0
        while j < len(variant_ids):
            if get_vatiant_id(vcf_data[i]) == variant_ids[j]:
                # count number of vep entires per variant that result in an esm score (i.e. with consequence "missense")
                numberOfEsmScoresPerVariant = 0
                for l in range(j, len(variant_ids), 1):
                    if get_vatiant_id(vcf_data[i]) == variant_ids[l]:
                        numberOfEsmScoresPerVariant = numberOfEsmScoresPerVariant + 1
                    else:
                        break

                # annotate vcf line with esm scores
                if "EsmScoreMissense" in vcf_data[i]:
                    vcf_data[i] = vcf_data[i][:-1] + "," + vcf_data[i][-1:]
                else:
                    vcf_data[i] = (
                        vcf_data[i][:-1] + ";EsmScoreMissense" + "=" + vcf_data[i][-1:]
                    )
                for h in range(0, numberOfEsmScoresPerVariant, 1):
                    if aa_seq[j + h] != "NA":
                        average_score = 0
                        for k in range(0, numModels, 1):
                            average_score = average_score + float(modelScores[k][j + h])
                        average_score = average_score / numModels
                        vcf_data[i] = (
                            vcf_data[i][:-1]
                            + str(transcript_ids[j + h][11:])
                            + "|"
                            + str(round(float(average_score), 3))
                            + vcf_data[i][-1:]
                        )
                    else:
                        vcf_data[i] = (
                            vcf_data[i][:-1]
                            + str(transcript_ids[j + h][11:])
                            + "|"
                            + "NA"
                            + vcf_data[i][-1:]
                        )
                    if h != numberOfEsmScoresPerVariant - 1:
                        vcf_data[i] = vcf_data[i][:-1] + "," + vcf_data[i][-1:]

                j = j + numberOfEsmScoresPerVariant
            else:
                j = j + 1

    return(vcf_data)

def write_variants(output_file, vcf_data):

    vcf_file_output = BgzfWriter(output_file, "w")
    for line in vcf_data:
        vcf_file_output.write(line)

    vcf_file_output.close()
