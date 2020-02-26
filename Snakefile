# Variables and filepaths stored here
include: "scripts/snakemake_variables.py"

shell.executable("/usr/bin/bash")
shell.prefix("source ~/.bashrc; ")

# TODO use glob + expand to properly specify input of make_anc_bed and partition_anc_het

rule all:
    input:
        DATA_DIR + "hits/hits_ascertainment_Eur.txt",
        DATA_DIR + "hits/hits_estimation_Eur.txt",
        DATA_DIR + "hits/hits_estimation_het.txt",
        DATA_DIR + "hits/hits_estimation_Afr.txt"

########################### GENERATE GENE ANNOTATION ###########################

rule annotate_genes:
    input:
        DATA_DIR + GENCODE_ANNO 
    output:
        DATA_DIR + "gencode.v30.genes.gtf"
    shell:
        """
        conda activate gene-annotation-env
        python {ANNO_SCRIPT} {input} {output}
        conda deactivate
        """

################################ CHOOSE SAMPLES ################################

rule make_geno_list:
    input:
        DATA_DIR + VCF 
    output:
        DATA_DIR + "genotyped_individual_IDs.txt"
    shell:
        """
        zgrep -m 1 "^#CHROM" {input} | sed 's/   /\n/g' > {output}
        """

rule make_exp_list:
    input:
        DATA_DIR + READS 
    output:
        DATA_DIR + "rnaseqd_individual_IDs.txt"
    shell:
        """
        zgrep -m 1 "Name" {input} | cut -f 1,2 --complement | tr '\t' '\n' > {output}
        """

rule select_samples:
    input:
        sample_data=DATA_DIR + RNASEQ_METADATA,
        geno_ID=rules.make_geno_list.output,
        exp_ID=rules.make_exp_list.output,
        ind_data=DATA_DIR + IND_METADATA 
    output:
        DATA_DIR + "sample_participant_lookup.txt"
    shell:
        """
        conda activate py36
        python {SAMPLE_SELECTION_SCRIPT} {input.sample_data} {input.ind_data} {input.geno_ID} {input.exp_ID} {output}
        conda deactivate
        """

############################# NORMALIZE EXPRESSION #############################

rule make_chr_list:
    input:
        vcf=ancient(DATA_DIR + VCF)
    output:
        chr_list=DATA_DIR + "chr_list.txt"
    shell:
        """
        conda activate tabix-env
        tabix -l {input.vcf} > {output.chr_list}
        conda deactivate
        """

rule normalize_expression:
    input:
        tpm=DATA_DIR + TPM,
        counts=DATA_DIR + READS,
        anno=rules.annotate_genes.output,
        sample_map=rules.select_samples.output,
        chr_list=rules.make_chr_list.output
    output:
        bed=DATA_DIR + "mesa.expression.bed.gz",
        idx=DATA_DIR + "mesa.expression.bed.gz.tbi"
    params:
        prefix=DATA_DIR + "mesa"
    shell:
        """
        conda activate norm-exp-env
        {GTEX_NORM_SCRIPT} {input.tpm} {input.counts} {input.anno} \
            {input.sample_map} {input.chr_list} {params.prefix} \
            --tpm_threshold 0.1 \
            --count_threshold 6 \
            --sample_frac_threshold 0.2 
        conda deactivate
        """

rule make_gene_list:
    input:
        ancient(rules.normalize_expression.output.bed)
    output:
        DATA_DIR + "all_genes_TSS.txt"
    shell:
        """
        less {input} | cut -f 1,2,3,4 | tail -n +2 | sed 's/chr//' > {output}
        """

############################# GENERATE COVARIATES #############################

rule MAF_filter:
    input:
        ancient(DATA_DIR + VCF)
    output:
        vcf=DATA_DIR + "mesa.maf05.vcf.gz"
    shell:
        """
        conda activate bcftools-env
        bcftools view -i 'MAF[0] > 0.05' -Ou {input} | \
            bcftools view --genotype ^miss --phased -Oz -o {output.vcf}
        conda deactivate
        """

rule find_indep_snps:
    input:
        vcf=rules.MAF_filter.output
    output:
        keep=DATA_DIR + "mesa.maf05.prune.in",
        exclude=temp(expand(DATA_DIR + "mesa.maf05.{ext}", ext=["log", "nosex", "prune.out"]))
    params:
        prefix=lambda wildcards, output: output.keep[:-9]
    shell:
        """
        conda activate plink-env
        plink --vcf {input} --vcf-half-call m --allow-no-sex \
            --keep-allele-order --indep-pairwise 50 5 0.5 --out {params.prefix}
        conda deactivate
        """

rule prune:
    input:
        vcf=rules.MAF_filter.output,
        snps=rules.find_indep_snps.output.keep
    output:
        DATA_DIR + "mesa.maf05.prune.vcf.gz"
    shell:
        """
        conda activate bcftools-env
        bcftools view -i 'ID=@{input.snps}' -Oz -o {output} {input.vcf}
        conda deactivate
        """

rule prep_covariates:
    input:
        vcf=rules.prune.output,
        exp=rules.normalize_expression.output.bed,
        samples=rules.select_samples.output,
        indiv_data=DATA_DIR + IND_METADATA,
        sample_data=DATA_DIR + RNASEQ_METADATA 
    output:
        DATA_DIR + "mesa.sample_covariates.txt"
    shell:
        """
        conda activate py36
        python {COV_PREP_SCRIPT} --samples {input.samples} \
            --indiv_metadata {input.indiv_data} \
            --sample_metadata {input.sample_data} \
            --out {output}
        conda deactivate
        """

############################## PARTITION SAMPLES ##############################

rule make_anc_bed:
    output:
        DATA_DIR + "anc_tracts.bed"
    params:
        dir=DATA_DIR + "bed/"
    shell:
        """
        conda activate py36
        python {TRACT_SCRIPT} --tract_dir {params.dir} --out {output}
        conda deactivate
        """

rule make_genes_bed:
    input:
        genes=rules.make_gene_list.output,
        anno=rules.annotate_genes.output,
        chrom_lengths=DATA_DIR + CHR_LEN 
    output:
        DATA_DIR + "genes.bed"
    shell:
        """
        conda activate py36
        python {GENE_BED_SCRIPT} --genes {input.genes} \
            --anno {input.anno} --chrom_lengths {input.chrom_lengths} \
            --out {output}
        conda deactivate
        """

rule intersect_tracts_genes:
    input:
        tracts=rules.make_anc_bed.output,
        genes=rules.make_genes_bed.output
    output:
        DATA_DIR + "intersection_anc_genes.bed"
    shell:
        """
        conda activate bedtools
        sort -k 1,1 -k2,2n {input.tracts} | bedtools intersect -a stdin -b {input.genes} -wao > {output}
        conda deactivate
        """

checkpoint partition_samples:
    input:
        intersect=rules.intersect_tracts_genes.output,
        samples=rules.select_samples.output,
        metadata=DATA_DIR + IND_METADATA,
        genes=rules.make_gene_list.output
    output:
        directory(DATA_DIR + "fastqtl_sample_input")
    params:
        output_dir=DATA_DIR + "fastqtl_sample_input"
    shell:
        """
        mkdir -p {params.output_dir}/ascertainment/Eur
        mkdir -p {params.output_dir}/estimation/Afr
        mkdir -p {params.output_dir}/estimation/het
        mkdir -p {params.output_dir}/estimation/Eur
        conda activate py36
        python {PARTITION_SCRIPT} --intersect {input.intersect} \
            --samples {input.samples} --metadata {input.metadata} \
            --genes {input.genes} --out {params.output_dir}
        conda deactivate
        """

################################## CALL EQTLS ##################################

rule subset_geno:
    input:
        vcf=DATA_DIR + VCF,
        gene_TSS_map=rules.make_gene_list.output
    output:
        temp(DATA_DIR + "fastqtl_geno_input/{gene}.vcf")
    params:
        output_dir=DATA_DIR + "fastqtl_geno_input/"
    shell:
        """
        mkdir -p {params.output_dir}
        CHR=$(grep {wildcards.gene} {input.gene_TSS_map} | cut -f 1)
        REGION_START=$(grep {wildcards.gene} {input.gene_TSS_map} | cut -f 2)
        REGION_STOP=$(grep {wildcards.gene} {input.gene_TSS_map} | cut -f 3)
        
        conda activate bcftools-env
        bcftools view -r chr$CHR:$REGION_START-$REGION_STOP -Ov -o {output} {input.vcf}
        conda deactivate
        """

rule call_eQTLs:
    input:
        geno_input=rules.subset_geno.output,
        pheno_input=rules.normalize_expression.output.bed,
        sample_input=DATA_DIR + "fastqtl_sample_input/ascertainment/{anc}/{gene}.txt",
        covariates=rules.combine_covariates.output
    output:
        DATA_DIR + "fastqtl_output/{type}/{anc}/{gene}.txt"
    params:
        output_dir=DATA_DIR + "fastqtl_output/{type}/{anc}"
    shell:
        """
        mkdir -p {params.output_dir}

        conda activate pystats
        cat <(zcat {input.pheno_input} | head -n 1) <(zgrep -m 1 {wildcards.gene} {input.pheno_input}) \
            | python {EQTL_SCRIPT}
            --genotypes {input.geno_input} \
            --samples {input.sample_input} \
            --covariates {input.covariates} \
            --out {output}
        conda deactivate
        """

def merge_input(wildcards):
    checkpoint_output = checkpoints.partition_samples.get(**wildcards).output[0]
    return expand(DATA_DIR + "fastqtl_output/{type}/{anc}/{gene}.txt",
           type=wildcards.type, anc=wildcards.anc,
           gene=glob_wildcards(os.path.join(checkpoint_output, wildcards.type, wildcards.anc, "{gene}.txt")).gene)

rule merge_output:
    input:
        merge_input
    output:
        DATA_DIR + "fastqtl_output/merged_{type}_{anc}.txt"
    params:
        dir=DATA_DIR + "fastqtl_output/{type}/{anc}/*"
    shell:
        """
        cat {params.dir} > {output}
        """

rule identify_hits:
    input:
        asc=DATA_DIR + "fastqtl_output/merged_ascertainment_Eur.txt", 
        eur=DATA_DIR + "fastqtl_output/merged_estimation_Eur.txt",
        het=DATA_DIR + "fastqtl_output/merged_estimation_het.txt",
        afr=DATA_DIR + "fastqtl_output/merged_estimation_Afr.txt"
    output:
        DATA_DIR + "hits/hits_ascertainment_Eur.txt",
        DATA_DIR + "hits/hits_estimation_Eur.txt",
        DATA_DIR + "hits/hits_estimation_het.txt",
        DATA_DIR + "hits/hits_estimation_Afr.txt"
    params:
        DATA_DIR + "hits/"
    shell:
        """
        mkdir -p {params}
        module reset
        module load R
        Rscript {HITS_SCRIPT} --asc {input.asc} --eur {input.eur} \
            --het {input.het} --afr {input.afr} --out {params}
        conda deactivate
        """

################################# SIMULATIONS ################################

rule perform_simulation:
    input:
        sd=DATA_DIR + "expression_std_dev.txt",
        afr="results/v2/hits/hits_estimation_Afr.txt",
        eur="results/v2/hits/hits_estimation_Eur.txt" 
    output:
        DATA_DIR + "sims/scaling_{afr}_{eur}/sim_{sim}.txt"
    shell:
        """
        mkdir -p DATA_DIR/sims
        module reset
        module load R
        Rscript --vanilla {SIM_SCRIPT} --sd_file {input.sd} \
		--afr_effect_file {input.afr} \
		--eur_effect_file {input.eur} \
                --afr_scaling {wildcards.afr} \
                --eur_scaling {wildcards.eur} \
		--out {output}
        """

rule merge_sims:
    input:
        expand(DATA_DIR + "sims/scaling_{{afr}}_{{eur}}/sim_{sim}.txt", sim=range(SIM_ITER))
    output:
        DATA_DIR + "sims/scaling_{afr}_{eur}/merged_results.txt"
    shell:
        """
        cat {input} > {output}
        """
