
"""Note that we are declaring many files temporary here that should be located in a temporary folder to begin with"""

rule decompress:
    input: '{file}.vcf.gz'
    output: temp('{file}.vcf')
    shell:
        '''
        zcat {input} > {output}
        '''

rule prepare:
    input: '{file}.vcf'
    output: temp('{file}.prepared.vcf')
    conda: 'envs/environment.yml'
    shell:
        '''
        cat {input} \
        | python $CADD/src/scripts/VCF2vepVCF.py \
        | sort -k1,1 -k2,2n -k4,4 -k5,5 \
        | uniq > {output}
        '''        

rule prescore:
    input: '{file}.prepared.vcf'
    output:
        novel=temp('{file}.novel.vcf'),
        prescored=temp('{file}.pre.tsv')
    conda: 'envs/environment.yml'
    shell:
        '''
        # Prescoring
        echo '## Prescored variant file' > {output.prescored};
        if [ -d $CADD/{config[PrescoredFolder]} ]
        then
            for PRESCORED in $(ls $CADD/{config[PrescoredFolder]}/*.tsv.gz)
            do
                cat {input} \
                | python $CADD/src/scripts/extract_scored.py --header \
                    -p $PRESCORED --found_out={output.prescored}.tmp \
                > {input}.tmp;
                cat {output.prescored}.tmp >> {output.prescored}
                mv {input}.tmp {input};
            done;
            rm {output.prescored}.tmp
        fi
        mv {input} {output.novel}
        '''

rule annotation:
    input: '{file}.novel.vcf'
    output: temp('{file}.anno.tsv.gz')
    conda: 'envs/environment.yml'
    shell:
        '''
        cat {input} \
        | vep --quiet --cache --offline --dir $CADD/{config[VEPpath]} \
            --buffer 1000 --no_stats --species homo_sapiens \
            --db_version={config[EnsemblDB]} --assembly {config[GenomeBuild]} \
            --format vcf --regulatory --sift b --polyphen b --per_gene --ccds --domains \
            --numbers --canonical --total_length --vcf --force_overwrite --output_file STDOUT \
        | python $CADD/src/scripts/annotateVEPvcf.py \
            -c $CADD/{config[ReferenceConfig]} \
        | gzip -c > {output}
        '''

rule imputation:
    input: '{file}.anno.tsv.gz'
    output: temp('{file}.csv.gz')
    conda: 'envs/environment.yml'
    shell:
        '''
        zcat {input} \
        | python $CADD/src/scripts/trackTransformation.py -b \
            -c $CADD/{config[ImputeConfig]} -o {output} --noheader;
        '''

rule score:
    input:
        impute='{file}.csv.gz',
        anno='{file}.anno.tsv.gz'
    output: temp('{file}.novel.tsv')
    conda: 'envs/environment.yml'
    shell:
        '''
        python $CADD/src/scripts/predictSKmodel.py \
            -i {input.impute} -m $CADD/{config[Model]} -a {input.anno} \
        | python $CADD/src/scripts/max_line_hierarchy.py --all \
        | python $CADD/src/scripts/appendPHREDscore.py \
            -t $CADD/{config[ConversionTable]} > {output};
    
        if [ "{config[Annotation]}" = 'False' ]
        then
            cat {output} | cut -f {config[Columns]} | uniq > {output}.tmp
            mv {output}.tmp {output}
        fi
        '''

rule join:
    input:
        pre='{file}.pre.tsv',
        novel='{file}.novel.tsv'
    output: '{file}.tsv.gz'
    conda: 'envs/environment.yml'
    shell:
        '''
        (
        echo "{config[Header]}";
        head -n 1 {input.novel};
        cat {input.pre} {input.novel} \
        | grep -v "^#" \
        | sort -k1,1 -k2,2n -k3,3 -k4,4 || true;
        ) | bgzip -c > {output};
        '''
