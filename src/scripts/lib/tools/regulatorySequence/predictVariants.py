import click


# options
@click.command()
@click.option(
    "--variants",
    "variants_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Variant file to predict in VCF format.",
)
@click.option(
    "--model",
    "model_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Tensorflow model in json format.",
)
@click.option(
    "--weights",
    "weights_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Model weights in hdf5 format.",
)
@click.option(
    "--reference",
    "reference_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Reference sequence in FASTA format (indexed).",
)
@click.option(
    "--genome",
    "genome_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Genome file of the reference with lengths of contigs.",
)
@click.option(
    "--output",
    "output_file",
    type=click.Path(writable=True),
    default="/dev/stdout",
    help="Output file with predictions in tsv.gz format.",
)
def cli(
    variants_file, model_file, weights_file, reference_file, genome_file, output_file
):
    import os
    os.environ["CUDA_VISIBLE_DEVICES"]="-1"    
    import numpy as np
    import math
    import vcfpy

    import copy

    import tensorflow as tf

    from seqiolib import Interval, Encoder, VariantType, Variant
    from seqiolib import utils

    from pyfaidx import Fasta
    import pybedtools

    strategy = tf.distribute.MirroredStrategy()

    def loadAndPredict(sequences, model, variants=None):
        X = []
        i = 0
        for sequence in sequences:
            if variants is not None:
                sequence.replace(variants[i])
            X.append(Encoder.one_hot_encode_along_channel_axis(sequence.getSequence()))
            i += 1
        prediction = model.predict(np.array(X))
        return prediction

    def extendIntervals(intervals, region_length, genome_file):
        left = math.ceil((region_length - 1) / 2)
        right = math.floor((region_length - 1) / 2)
        click.echo("Extending intervals left=%d, right=%d..." % (left, right))
        return list(
            map(
                pybedtoolsIntervalToInterval,
                intervals.slop(r=right, l=left, g=str(genome_file)),
            )
        )

    def getCorrectedChrom(chrom):
        if chrom.startswith("chr"):
            return chrom
        elif chrom == "MT":
            return "chrM"
        else:
            return "chr" + chrom

    def variantToPybedtoolsInterval(record):
        return pybedtools.Interval(getCorrectedChrom(record.CHROM), record.POS - 1, record.POS)

    def pybedtoolsIntervalToInterval(interval_pybed):
        return Interval(
            interval_pybed.chrom, interval_pybed.start + 1, interval_pybed.stop
        )

    # load variants
    click.echo("Loading variants...")
    records = []
    vcf_reader = vcfpy.Reader.from_path(variants_file)

    for record in vcf_reader:
        records.append(record)
    click.echo("Found %d variants" % len(records))

    if len(records) == 0:
        click.echo("No variants found. Writing file with header only and exiting...")
        vcf_writer = vcfpy.Writer.from_path(output_file, vcf_reader.header)
        vcf_writer.close()
        exit(0)
    # convert to intervals (pybedtools)
    click.echo("Convert to bed tools intervals...")
    intervals = pybedtools.BedTool(list(map(variantToPybedtoolsInterval, records)))

    with strategy.scope():
        click.echo("Load model...")
        model = utils.io.ModelIO.loadModel(model_file, weights_file)

        input_length = model.input_shape[1]
        click.echo("Detecting interval length of %d" % input_length)
        intervals = extendIntervals(intervals, input_length, genome_file)

        # load sequence for variants
        reference = Fasta(reference_file)
        sequences_ref = []
        sequences_alt = []
        predict_avail_idx = set()

        click.echo("Load reference and try to get ref and alt.")
        alt_idx = 0
        for i in range(len(records)):
            record = records[i]
            interval = intervals[i]

            # can be problematic if we are on the edges of a chromose.
            # Workaround. It is possible to extend the intreval left or right to get the correct length
            if interval.length != input_length:
                click.echo(
                    "Cannot use variant %s because of wrong size of interval %s "
                    % (str(record), str(interval))
                )
                alt_idx += len(record.ALT)
                continue

            sequence_ref = utils.io.SequenceIO.readSequence(reference, interval)

            for j in range(len(record.ALT)):
                alt_record = record.ALT[j]
                variant = Variant(
                    getCorrectedChrom(record.CHROM), record.POS, record.REF, alt_record.value
                )
                # INDEL
                if (
                    variant.type == VariantType.DELETION
                    or variant.type == VariantType.INSERTION
                ):
                    # DELETION
                    if variant.type == VariantType.DELETION:
                        extend = len(variant.ref) - len(variant.alt)
                        if interval.isReverse():
                            interval.position = interval.position + extend
                        else:
                            interval.position = interval.position - extend
                        interval.length = interval.length + extend
                    # INSERTION
                    elif variant.type == VariantType.INSERTION:
                        extend = len(variant.alt) - len(variant.ref)
                        if interval.isReverse():
                            interval.position = interval.position - extend
                        else:
                            interval.position = interval.position + extend
                        interval.length = interval.length - extend
                    if interval.length > 0:
                        sequence_alt = utils.io.SequenceIO.readSequence(
                            reference, interval
                        )
                        sequence_alt.replace(variant)
                        if len(sequence_alt.sequence) == input_length:
                            # FIXME: This is a hack. it seems that for longer indels the replacement does not work
                            sequences_alt.append(sequence_alt)
                            sequences_ref.append(sequence_ref)
                            predict_avail_idx.add(alt_idx)
                        else:
                            print(
                                "Cannot use variant %s because of wrong interval %s has wrong size after InDel Correction"
                                % (str(variant), str(interval))
                            )
                    else:
                        print(
                            "Cannot use variant %s because interval %s has negative size"
                            % (str(variant), str(interval))
                        )
                # SNV
                else:
                    sequence_alt = copy.copy(sequence_ref)
                    sequence_alt.replace(variant)
                    sequences_alt.append(sequence_alt)
                    sequences_ref.append(sequence_ref)
                    predict_avail_idx.add(alt_idx)
                alt_idx += 1
        click.echo("Predict reference...")
        results_ref = loadAndPredict(sequences_ref, model)
        click.echo("Predict alternative...")
        results_alt = loadAndPredict(sequences_alt, model)

    num_targets = results_ref.shape[1] if len(results_alt.shape) > 1 else 1

    for task_id in range(num_targets):
        vcf_reader.header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "RegSeq%d" % task_id),
                    ("Number", "A"),
                    ("Type", "Float"),
                    (
                        "Description",
                        "Regulatory sequence prediction of the alt minus reference, output task %d"
                        % task_id,
                    ),
                ]
            )
        )

    vcf_writer = vcfpy.Writer.from_path(output_file, vcf_reader.header)

    alt_idx = 0
    predict_idx = 0
    for i in range(len(records)):
        record = records[i]
        to_add = {}
        for j in range(len(record.ALT)):
            if alt_idx in predict_avail_idx:
                for task_id in range(num_targets):
                    to_add["RegSeq%d" % task_id] = to_add.get(
                        "RegSeq%d" % task_id, []
                    ) + [round(results_alt[predict_idx][task_id] - results_ref[predict_idx][task_id],6)]
                predict_idx += 1
            else:
                for task_id in range(num_targets):
                    to_add["RegSeq%d" % task_id] = to_add.get(
                        "RegSeq%d" % task_id, []
                    ) + [np.nan]
            alt_idx += 1

        for key, value in to_add.items():
            record.INFO[key] = value

        vcf_writer.write_record(record)
    vcf_writer.close()


if __name__ == "__main__":
    cli()
