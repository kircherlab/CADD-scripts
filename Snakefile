"""
Note that we are declaring many files temporary here that should be located in a temporary folder to begin with
"""


# Min version of snakemake
from snakemake.utils import min_version

min_version("7.32.3")

# Validation of config file
from snakemake.utils import validate

validate(config, schema="schemas/config_schema.yaml")

# CADD environment variable
import os


envvars:
    "CADD",


# START Rules


rule decompress:
    conda:
        "envs/environment_minimal.yml"
    input:
        "{file}.vcf.gz",
    output:
        temp("{file}.vcf"),
    log:
        "{file}.decompress.log",
    shell:
        """
        zcat {input} > {output} 2> {log}
        """


rule prepare:
    conda:
        "envs/environment_minimal.yml"
    input:
        "{file}.vcf",
    output:
        temp("{file}.prepared.vcf"),
    log:
        "{file}.prepare.log",
    params:
        cadd=os.environ["CADD"],
    shell:
        """
        cat {input} \
        | python {params.cadd}/src/scripts/VCF2vepVCF.py \
        | sort -k1,1 -k2,2n -k4,4 -k5,5 \
        | uniq > {output} 2> {log}
        """


rule prescore:
    conda:
        "envs/environment_minimal.yml"
    input:
        vcf="{file}.prepared.vcf",
        prescored="%s/%s" % (os.environ["CADD"], config["PrescoredFolder"]),
    output:
        novel=temp("{file}.novel.vcf"),
        prescored=temp("{file}.pre.tsv"),
    log:
        "{file}.prescore.log",
    params:
        cadd=os.environ["CADD"],
    shell:
        """
        # Prescoring
        echo '## Prescored variant file' > {output.prescored} 2> {log};
        PRESCORED_FILES=`find -L {input.prescored} -maxdepth 1 -type f -name \\*.tsv.gz | wc -l`
        if [ ${{PRESCORED_FILES}} -gt 0 ];
        then
            for PRESCORED in $(ls {input.prescored}/*.tsv.gz)
            do
                cat {input.vcf} \
                | python {params.cadd}/src/scripts/extract_scored.py --header \
                    -p $PRESCORED --found_out={output.prescored}.tmp \
                > {input.vcf}.tmp 2>> {log};
                cat {output.prescored}.tmp >> {output.prescored}
                mv {input.vcf}.tmp {input.vcf} &> {log};
            done;
            rm {output.prescored}.tmp &>> {log}
        fi
        mv {input.vcf} {output.novel} &>> {log}
        """


rule annotation_vep:
    conda:
        "envs/vep.yml"
    input:
        vcf="{file}.novel.vcf",
        veppath="%s/%s" % (os.environ["CADD"], config["VEPpath"]),
    output:
        temp("{file}.vep.vcf.gz"),
    log:
        "{file}.annotation_vep.log",
    params:
        cadd=os.environ["CADD"],
        genome_build=config["GenomeBuild"],
        ensembl_db=config["EnsemblDB"],
    shell:
        """
        cat {input.vcf} \
        | vep --quiet --cache --offline --dir {input.veppath} \
            --buffer 1000 --no_stats --species homo_sapiens \
            --db_version={params.ensembl_db} --assembly {params.genome_build} \
            --format vcf --regulatory --sift b --polyphen b --per_gene --ccds --domains \
            --numbers --canonical --total_length --vcf --force_overwrite --output_file STDOUT \
        | bgzip -c > {output} 2> {log}
        """


rule annotate_esm:
    conda:
        "envs/esm.yml"
    input:
        vcf="{file}.vep.vcf.gz",
        models=expand(
            "{path}/{model}.pt",
            path=config["ESMpath"],
            model=config["ESMmodels"],
        ),
        transcripts="%s/pep.%s.fa" % (config["ESMpath"], config["EnsemblDB"]),
    output:
        missens=temp("{file}.esm_missens.vcf.gz"),
        frameshift=temp("{file}.esm_frameshift.vcf.gz"),
        final=temp("{file}.esm.vcf.gz"),
    log:
        "{file}.annotate_esm.log",
    params:
        cadd=os.environ["CADD"],
        models=["--model %s " % model for model in config["ESMmodels"]],
        batch_size=config["ESMbatchsize"],
    shell:
        """
        model_directory=`dirname {input.models[0]}`;
        model_directory=`dirname $model_directory`;

        python {params.cadd}/src/scripts/lib/tools/esmScore/esmScore_missense_av_fast.py \
        --input {input.vcf} \
        --transcripts {input.transcripts} \
        --model-directory $model_directory {params.models} \
        --output {output.missens} --batch-size {params.batch_size} &> {log}

        python {params.cadd}/src/scripts/lib/tools/esmScore/esmScore_frameshift_av.py \
        --input {output.missens} \
        --transcripts {input.transcripts} \
        --model-directory $model_directory {params.models} \
        --output {output.frameshift} --batch-size {params.batch_size} &>> {log}

        python {params.cadd}/src/scripts/lib/tools/esmScore/esmScore_inFrame_av.py \
        --input {output.frameshift} \
        --transcripts {input.transcripts} \
        --model-directory $model_directory {params.models} \
        --output {output.final} --batch-size {params.batch_size} &>> {log}
        """


rule annotate_regseq:
    conda:
        "envs/regulatorySequence.yml"
    input:
        vcf="{file}.esm.vcf.gz",
        reference="%s/%s" % (config["REGSEQpath"], "reference.fa"),
        genome="%s/%s" % (config["REGSEQpath"], "reference.fa.genome"),
        model="%s/%s" % (config["REGSEQpath"], "Hyperopt400InclNegatives.json"),
        weights="%s/%s" % (config["REGSEQpath"], "Hyperopt400InclNegatives.h5"),
    output:
        temp("{file}.regseq.vcf.gz"),
    log:
        "{file}.annotate_regseq.log",
    params:
        cadd=os.environ["CADD"],
    shell:
        """
        python {params.cadd}/src/scripts/lib/tools/regulatorySequence/predictVariants.py \
        --variants {input.vcf} \
        --model {input.model} \
        --weights {input.weights} \
        --reference {input.reference} \
        --genome {input.genome} \
        --output {output} &> {log}
        """


rule annotate_mmsplice:
    conda:
        "envs/mmsplice.yml"
    input:
        vcf="{file}.regseq.vcf.gz",
        transcripts="%s/homo_sapiens.110.gtf" % config.get("MMSPLICEpath", ""),
        reference="%s/reference.fa" % config.get("REFERENCEpath", ""),
    output:
        mmsplice=temp("{file}.mmsplice.vcf.gz"),
        idx=temp("{file}.regseq.vcf.gz.tbi"),
    log:
        "{file}.annotate_mmsplice.log",
    params:
        cadd=os.environ["CADD"],
    shell:
        """
        tabix -p vcf {input.vcf} &> {log};
        python {params.cadd}/src/scripts/lib/tools/MMSplice.py -i {input.vcf} \
        -g {input.transcripts} \
        -f {input.reference} | \
        grep -v '^Variant(CHROM=' | \
        bgzip -c > {output.mmsplice} 2>> {log}
        """


rule annotation:
    conda:
        "envs/environment_minimal.yml"
    input:
        vcf=lambda wc: "{file}.%s.vcf.gz"
        % ("mmsplice" if config["GenomeBuild"] == "GRCh38" else "regseq"),
        reference_cfg="%s/%s" % (os.environ["CADD"], config["ReferenceConfig"]),
    output:
        temp("{file}.anno.tsv.gz"),
    log:
        "{file}.annotation.log",
    params:
        cadd=os.environ["CADD"],
    shell:
        """
        zcat {input.vcf} \
        | python {params.cadd}/src/scripts/annotateVEPvcf.py \
            -c {input.reference_cfg} \
        | gzip -c > {output} 2> {log}
        """


rule imputation:
    conda:
        "envs/environment_minimal.yml"
    input:
        tsv="{file}.anno.tsv.gz",
        impute_cfg="%s/%s" % (os.environ["CADD"], config["ImputeConfig"]),
    output:
        temp("{file}.csv.gz"),
    log:
        "{file}.imputation.log",
    params:
        cadd=os.environ["CADD"],
    shell:
        """
        zcat {input.tsv} \
        | python {params.cadd}/src/scripts/trackTransformation.py -b \
            -c {input.impute_cfg} -o {output} --noheader &>> {log};
        """


rule score:
    conda:
        "envs/environment_minimal.yml"
    input:
        impute="{file}.csv.gz",
        anno="{file}.anno.tsv.gz",
        conversion_table="%s/%s" % (os.environ["CADD"], config["ConversionTable"]),
        model_file="%s/%s" % (os.environ["CADD"], config["Model"]),
    output:
        temp("{file}.novel.tsv"),
    log:
        "{file}.score.log",
    params:
        cadd=os.environ["CADD"],
        use_anno=config["Annotation"],
        columns=config["Columns"],
    shell:
        """
        python {params.cadd}/src/scripts/predictSKmodel.py \
            -i {input.impute} -m {input.model_file} -a {input.anno} \
        | python {params.cadd}/src/scripts/max_line_hierarchy.py --all \
        | python {params.cadd}/src/scripts/appendPHREDscore.py \
            -t {input.conversion_table} > {output} 2>> {log};
    
        if [ "{params.use_anno}" = 'False' ]
        then
            cat {output} | cut -f {params.columns} | uniq > {output}.tmp 2>> {log};
            mv {output}.tmp {output} &>> {log}
        fi
        """


rule join:
    conda:
        "envs/environment_minimal.yml"
    input:
        pre="{file}.pre.tsv",
        novel="{file}.novel.tsv",
    output:
        "{file}.tsv.gz",
    log:
        "{file}.join.log",
    params:
        header=config["Header"],
    shell:
        """
        (
            echo "{params.header}";
            head -n 1 {input.novel};
            cat {input.pre} {input.novel} | \
            grep -v "^#" | \
            sort -k1,1 -k2,2n -k3,3 -k4,4 || true;
        ) | bgzip -c > {output} 2>> {log};
        """


# END Rules
