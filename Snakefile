"""
Note that we are declaring many files temporary here that should be located in a temporary folder to begin with
"""


# container with conda environments
containerized: "docker://visze/cadd-scripts-v1_7:0.1.0"

# need glob to get chunked files
import psutil 

# Min version of snakemake
from snakemake.utils import min_version

min_version("7.32.3")

# Validation of config file
from snakemake.utils import validate

validate(config, schema="schemas/config_schema.yaml")

# CADD environment variable
import os

###################################
## RULE SPECFIC THREADING LIMITS ##
###################################
system_memory = psutil.virtual_memory().total / (1024 ** 3)
gpu_memory = 0
try:
    lines = subprocess.check_output("nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null", shell=True).decode('utf-8').splitlines()
    # sum over gpu(s)
    gpu_memory = int(sum([int(x) for x in lines])/1024)
except Exception as e:
    pass
    
# esm takes ~16 Gb of system memory per tread & 4 Gb of GPU ram (if available). 
if gpu_memory > 0:
    config['esm_slots'] = int(min((0.9*system_memory)/16, 0.95*gpu_memory/4))
else:
    config['esm_slots'] = int(0.9*system_memory / 16)
if config['esm_slots'] < 1:
    config['esm_slots'] = 1
# then assign other resources
config['esm_load'] = 100 if config['esm_slots'] < 1 else int(100/config['esm_slots']) 
config['esm_threads'] = int(workflow.cores /  config['esm_slots'])
config['vep_load'] = 100 if system_memory < 4 else int(100/int(0.9*system_memory / 4 ))  # up to 4Gb/ram
config['regseq_load'] = 100 if system_memory < 2 else int(100/int(0.9*system_memory / 2 ))   # up to 2Gb of ram
config['mms_load'] = 100 if system_memory < 16 else int(100/int(0.9*system_memory / 16 ))   # up to 16Gb/ram
config['mms_threads'] = int(workflow.cores / (100/config['mms_load']))
config['anno_load'] = 1 # disk IO intensive
config['impute_load'] = 1 
config['prescore_load'] = 1 
config['score_load'] = 1

print("Threading Overview",flush=True)
print("##################",flush=True)
print("Assigned cores: {}".format(workflow.cores),flush=True)
print("Available system memory: {}GB".format(int(system_memory)),flush=True)
if gpu_memory > 0:
   print("Total GPU memory: {}GB".format(gpu_memory),flush=True)
else:
   print("No gpu found",flush=True)
print("Task Parallelization: ",flush=True)
print("  PreScore : {}x".format(min(workflow.cores,int(100/config['prescore_load']))))
print("  VEP : {}x".format(min(workflow.cores,int(100/config['vep_load']))))
print("  ESM : {}x with {} threads each (memory/gpu constraints)".format(min(workflow.cores,config['esm_slots']),config['esm_threads']))
print("  RegSeq : {}x".format(min(workflow.cores,int(100/config['regseq_load']))))
print("  MMsplice : {}x with {} threads each (memory constraints)".format(min(workflow.cores,int(100/config['mms_load'])),config['mms_threads']))
print("  Annotate : {}x".format(min(workflow.cores,int(100/config['anno_load']))))
print("  Impute : {}x".format(min(workflow.cores,int(100/config['impute_load']))))
          



## allowed scattering 
scattergather:
    split=workflow.cores,

envvars:
    "CADD",


#wildcard_constraints:
#    basefile="[^/]+"

# START Rules

rule decompress:
    conda:
        "envs/environment_minimal.yml"
    input:
        "{file}.vcf.gz",
    output:
        "{file}.vcf",
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
        vcf="{file}.vcf",
        
    output:
        #prep="{file}.prepared.vcf.tmp",
        #split=directory("{file}_splits"),
        splits=scatter.split("{{file}}_splits/chunk_{scatteritem}.prepared.vcf"),
    log:
        "{file}.prepare.log",
    params:
        cadd=os.environ["CADD"],
        threads=workflow.cores,
    resources:
        # < 1GB of memory
        load=1,
    shell:
        """
        mkdir -p {wildcards.file}_splits/ 2>> {log}
        cat {input.vcf} \
        | python {params.cadd}/src/scripts/VCF2vepVCF.py \
        | grep -v '^#' \
        | sed 's/^chr//' \
        | sort -k1,1 -k2,2n -k4,4 -k5,5 \
        | uniq > {wildcards.file}_splits/full.vcf 2> {log} 

        # split
        LC=$(wc -l {wildcards.file}_splits/full.vcf | cut -f1 -d' ')
        LC=$(((LC / {params.threads})+1))
        
        split -l $LC --numeric-suffixes=1 --additional-suffix="-of-{params.threads}.prepared.vcf" {wildcards.file}_splits/full.vcf {wildcards.file}_splits/chunk_ 2>> {log}

        rm -f {wildcards.file}_splits/full.vcf
        
        # strip padding zeros in the file names 
        for f in {wildcards.file}_splits/chunk_*.prepared.vcf
        do
            mv -n "$f" "$(echo "$f" | sed -E 's/(chunk_)0*([1-9][0-9]*)(-of-{params.threads}\\.prepared\\.vcf)/\\1\\2\\3/')"
        done
        """


checkpoint prescore:
    conda:
        "envs/environment_minimal.yml"
    input:
        vcf="{file}_splits/chunk_{chunk}.prepared.vcf",
        prescored="%s/%s" % (os.environ["CADD"], config["PrescoredFolder"]),
    output:
        novel="{file}_splits/chunk_{chunk}.novel.vcf",
        prescored="{file}_splits/chunk_{chunk}.pre.tsv",
    log:
        "{file}.chunk_{chunk}.prescore.log",
    params:
        cadd=os.environ["CADD"],
    resources:
        # < 1GB of memory
        load=int(config['prescore_load']),
    shell:
        """
        # Prescoring
        echo '## Prescored variant file' > {output.prescored} 2> {log};
        PRESCORED_FILES=`find -L {input.prescored} -maxdepth 1 -type f -name \\*.tsv.gz | wc -l`
        cp {input.vcf} {input.vcf}.new
        if [ ${{PRESCORED_FILES}} -gt 0 ];
        then
            for PRESCORED in $(ls {input.prescored}/*.tsv.gz)
            do
                cat {input.vcf}.new \
                | python {params.cadd}/src/scripts/extract_scored.py --header \
                    -p $PRESCORED --found_out={output.prescored}.tmp \
                > {input.vcf}.tmp 2>> {log};
                cat {output.prescored}.tmp >> {output.prescored}
                mv {input.vcf}.tmp {input.vcf}.new &> {log};
            done;
            rm {output.prescored}.tmp &>> {log}
        fi
        mv {input.vcf}.new {output.novel} &>> {log}
        """


rule annotation_vep:
    conda:
        "envs/vep.yml"
    input:
        vcf="{file}_splits/chunk_{chunk}.novel.vcf",
        veppath="%s/%s" % (os.environ["CADD"], config["VEPpath"]),
    output:
        "{file}_splits/chunk_{chunk}.vep.vcf.gz",
    log:
        "{file}.chunk_{chunk}.annotation_vep.log",
    params:
        cadd=os.environ["CADD"],
        genome_build=config["GenomeBuild"],
        ensembl_db=config["EnsemblDB"],
    resources:
        # < 1GB of memory
        load=int(config['vep_load']),
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
        #vcf="{file}_splits/chunk_{chunk}.vep.vcf.gz",
        vcf="{file}_splits/chunk_{chunk}.vep.vcf.gz",
        models=expand(
            "{path}/{model}.pt",
            path=config["ESMpath"],
            model=config["ESMmodels"],
        ),
        transcripts="%s/pep.%s.fa" % (config["ESMpath"], config["EnsemblDB"]),
    output:
        missens="{file}_splits/chunk_{chunk}.esm_missens.vcf.gz",
        frameshift="{file}_splits/chunk_{chunk}.esm_frameshift.vcf.gz",
        final="{file}_splits/chunk_{chunk}.esm.vcf.gz",
    log:
        "{file}.chunk_{chunk}.annotate_esm.log",
    resources:
        load=int(config['esm_load']),
    threads: 
        config['esm_threads'],
    params:
        cadd=os.environ["CADD"],
        models=["--model %s " % model for model in config["ESMmodels"]],
        batch_size=config["ESMbatchsize"],
        #header=config["Header"],
    
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

        #rm -f {wildcards.file}.esm_in.vcf.gz
        """


rule annotate_regseq:
    conda:
        "envs/regulatorySequence.yml"
    input:
        vcf="{file}_splits/chunk_{chunk}.esm.vcf.gz",
        reference="%s/%s" % (config["REGSEQpath"], "reference.fa"),
        genome="%s/%s" % (config["REGSEQpath"], "reference.fa.genome"),
        model="%s/%s" % (config["REGSEQpath"], "Hyperopt400InclNegatives.json"),
        weights="%s/%s" % (config["REGSEQpath"], "Hyperopt400InclNegatives.h5"),
    output:
        "{file}_splits/chunk_{chunk}.regseq.vcf.gz",
    log:
        "{file}.chunk_{chunk}.annotate_regseq.log",
    params:
        cadd=os.environ["CADD"],
    resources:
        # roughly 4GB of memory
        load=int(config['regseq_load']),
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
        vcf="{file}_splits/chunk_{chunk}.regseq.vcf.gz",
        transcripts="%s/homo_sapiens.110.gtf" % config.get("MMSPLICEpath", ""),
        reference="%s/reference.fa" % config.get("REFERENCEpath", ""),
    output:
        mmsplice="{file}_splits/chunk_{chunk}.mmsplice.vcf.gz",
        idx="{file}_splits/chunk_{chunk}.regseq.vcf.gz.tbi",
    log:
        "{file}.chunk_{chunk}.annotate_mmsplice.log",
    params:
        cadd=os.environ["CADD"],
    resources:
        load=int(config['mms_load']),
    threads: 
        config['mms_threads'],
    shell:
        """
        tabix -p vcf {input.vcf} &> {log};
        KERAS_BACKEND=tensorflow python {params.cadd}/src/scripts/lib/tools/MMSplice.py -i {input.vcf} \
        -g {input.transcripts} \
        -f {input.reference} | \
        grep -v '^Variant(CHROM=' | \
        bgzip -c > {output.mmsplice} 2>> {log}
        """


rule annotation:
    conda:
        "envs/environment_minimal.yml"
    input:
        vcf=lambda wc: "{file}_splits/chunk_{chunk}.%s.vcf.gz"
        % ("mmsplice" if config["GenomeBuild"] == "GRCh38" else "regseq"),
        reference_cfg="%s/%s" % (os.environ["CADD"], config["ReferenceConfig"]),
    output:
        "{file}_splits/chunk_{chunk}.anno.tsv.gz",
    log:
        "{file}.chunk_{chunk}.annotation.log",
    params:
        cadd=os.environ["CADD"],
    resources:
        load=int(config['anno_load']),
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
        tsv="{file}_splits/chunk_{chunk}.anno.tsv.gz",
        impute_cfg="%s/%s" % (os.environ["CADD"], config["ImputeConfig"]),
    output:
        "{file}_splits/chunk_{chunk}.csv.gz",
    log:
        "{file}.chunk_{chunk}.imputation.log",
    params:
        cadd=os.environ["CADD"],
    resources:
        load=int(config['impute_load']),
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
        impute="{file}_splits/chunk_{chunk}.csv.gz",
        anno="{file}_splits/chunk_{chunk}.anno.tsv.gz",
        conversion_table="%s/%s" % (os.environ["CADD"], config["ConversionTable"]),
        model_file="%s/%s" % (os.environ["CADD"], config["Model"]),
    output:
        "{file}_splits/chunk_{chunk}.novel.tsv",
    log:
        "{file}.chunk_{chunk}.score.log",
    params:
        cadd=os.environ["CADD"],
        use_anno=config["Annotation"],
        columns=config["Columns"],
    resources:
        load=config['score_load'],
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


# def aggregate_input(wildcards):
#     # Find all chunk files for the given wildcard
#     chunk_files = glob.glob(f"{wildcards.file}_splits/{wildcards.file}.chunk_*.novel.vcf")
#     pre_files = glob.glob(f"{wildcards.file}_splits/{wildcards.file}.chunk_*.pre.tsv")
#     
#     # Combine the novel and prescore chunk files if not empty
#     output = [f for f in chunk_files + pre_files if os.path.getsize(f) > 0]
#     if not output:
#         # no output : make empty file 
#         open(f"{wildcards.file}.empty", "w").close()
#         output = [f"{wildcards.file}.empty"]
# 
#     return output



#def aggregate_input(wildcards):
#    with checkpoints.prescore.get(file=wildcards.file).output["novel"].open() as f:
#        output = ["{file}.pre.tsv"]
#        for line in f:
#            if line.strip() != "":
#                output.append("{file}.novel.tsv")
#                break
#        return output


rule join:
    conda:
        "envs/environment_minimal.yml"
    input:
        #aggregate_input,
        pre=gather.split("{{file}}_splits/chunk_{scatteritem}.pre.tsv"),
        scored=gather.split("{{file}}_splits/chunk_{scatteritem}.novel.tsv"),
    output:
        "{file,.+(?<!\\.anno)}.tsv.gz",
    log:
        "{file}.join.log",
    params:
        header=config["Header"],
    shell:
        """
        (
            echo "{params.header}";
            cat {input.pre} {input.scored} | grep -v "^##" | grep "^#" | tail -n 1;
            cat {input.pre} {input.scored}| \
            grep -v "^#" | \
            sort -k1,1 -k2,2n -k3,3 -k4,4 || true;
        ) | bgzip -c > {output} 2>> {log};
        """


# END Rules
