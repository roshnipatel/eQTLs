# Variables and filepaths stored here
include: "scripts/snakemake_variables.py"

shell.executable("/usr/bin/bash")
shell.prefix("source ~/.bashrc; ")

# TODO use glob + expand to properly specify input of make_anc_bed and partition_anc_het

rule all:
    input:
        DATA_DIR + "fastqtl_output/merged_estimation_het.txt",
        DATA_DIR + "fastqtl_output/merged_estimation_Afr.txt",
        DATA_DIR + "fastqtl_output/merged_estimation_Eur.txt",
        DATA_DIR + "fastqtl_output/merged_ascertainment_Eur.txt"

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

rule make_rnaseq_list:
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
        rnaseq_ID=rules.make_rnaseq_list.output,
        ind_data=DATA_DIR + IND_METADATA 
    output:
        DATA_DIR + "sample_participant_lookup.txt"
    shell:
        """
        conda activate py36
        python {SAMPLE_SELECTION_SCRIPT} {input.sample_data} {input.ind_data} {input.geno_ID} {input.rnaseq_ID} {output}
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
            --sample_frac_threshold 0.2 \
            --normalization_method tmm
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

rule calculate_PEER:
    input:
        bed=ancient(rules.normalize_expression.output.bed)
    output:
        res=DATA_DIR + "mesa.PEER_residuals.txt",
        alpha=DATA_DIR + "mesa.PEER_alpha.txt",
        cov=DATA_DIR + "mesa.PEER_covariates.txt"
    params:
        prefix=DATA_DIR + "mesa",
        num_PEER=30
    shell:
        """
        Rscript {PEER_SCRIPT} {input.bed} {params.prefix} {params.num_PEER}
        """

rule MAF_filter:
    input:
        ancient(DATA_DIR + VCF)
    output:
        vcf=DATA_DIR + "mesa.maf05.vcf.gz"
    shell:
        """
        conda activate bcftools-env
        bcftools view -i 'MAF[0] > 0.05' -Ob {input} | \
            bcftools view --genotype ^miss --phased -Oz -o {output.vcf}
        conda deactivate
        """

rule indiv_filter:
    input:
        vcf=rules.MAF_filter.output.vcf,
        filter=rules.select_samples.output
    output:
        DATA_DIR + "mesa.maf05.filtered.vcf.gz"
    shell:
        """
        conda activate bcftools-env
        bcftools view -Oz -s $(tail -n +2 {input.filter} | cut -f2 | tr '\n' ',' | sed 's/,$//') \
            -o {output} {input.vcf}
        conda deactivate
        """

rule find_indep_snps:
    input:
        vcf=rules.indiv_filter.output
    output:
        keep=DATA_DIR + "mesa.maf05.filtered.prune.in",
        exclude=temp(expand(DATA_DIR + "mesa.maf05.filtered.{ext}", ext=["log", "nosex", "prune.out"]))
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
        vcf=rules.indiv_filter.output,
        snps=rules.find_indep_snps.output.keep
    output:
        bed=DATA_DIR + "mesa.maf05.filtered.prune05.bed",
        bim=temp(DATA_DIR + "mesa.maf05.filtered.prune05.bim"),
        fam=temp(DATA_DIR + "mesa.maf05.filtered.prune05.fam")
    params:
        prefix=lambda wildcards, output: output.bed[:-4]
    shell:
        """
        conda activate plink-env
        plink --vcf {input.vcf} --make-bed --extract {input.snps} --out {params.prefix}
        conda deactivate
        """

rule genotype_PC:
    input:
        rules.prune.output.bed,
        rules.prune.output.bim,
        rules.prune.output.fam
    output:
        raw=DATA_DIR + "mesa.maf05.filtered.prune05.pc.eigenvec",
        exclude=temp(expand(DATA_DIR + "mesa.maf05.filtered.prune05.pc.{ext}", ext=["eigenval", "log"])),
        trimmed=DATA_DIR + "mesa.maf05.filtered.prune05.pc.txt"
    params:
        out_prefix=lambda wildcards, output: output.trimmed[:-4],
        in_prefix=lambda wildcards, output: output.trimmed[:-7]
    shell:
        """
        conda activate plink-env
        plink --bfile {params.in_prefix} --pca --out {params.out_prefix}
        cut -f 2-22 -d$' ' {output.raw} > {output.trimmed}
        conda deactivate
        """

rule prep_sample_covariates:
    input:
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

rule combine_covariates:
    input:
        PEER=rules.calculate_PEER.output.cov,
        PC=rules.genotype_PC.output.trimmed,
        addl_cov=rules.prep_sample_covariates.output
    output:
        DATA_DIR + "mesa.combined_covariates.txt"
    params:
        prefix=DATA_DIR + "mesa",
    shell:
        """
        conda activate py36
        {COV_COMB_SCRIPT} {input.PEER} {params.prefix} \
            {input.PC} {input.addl_cov}
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

################################# RUN FASTQTL #################################

rule prep_pheno_file:
    output:
        temp(DATA_DIR + "fastqtl_gene_input/{gene}.txt")
    params:
        output_dir=DATA_DIR + "fastqtl_gene_input/"
    shell:
        """
        mkdir -p {params.output_dir}
        echo "{wildcards.gene}" > {output}
        """

rule ascertain_eQTLs:
    input:
        vcf=DATA_DIR + VCF,
        pheno=ancient(rules.normalize_expression.output.bed),
        sample_input=DATA_DIR + "fastqtl_sample_input/ascertainment/{anc}/{gene}.txt",
        gene_input=rules.prep_pheno_file.output,
        covariates=rules.combine_covariates.output,
        gene_TSS_map=rules.make_gene_list.output
    output:
        DATA_DIR + "fastqtl_output/ascertainment/{anc}/{gene}.txt"
    singularity:
        FASTQTL_DOCKER 
    params:
        output_dir=DATA_DIR + "fastqtl_output/"
    shell:
        """
        mkdir -p {params.output_dir}/ascertainment/{wildcards.anc}
        CHR=$(grep {wildcards.gene} {input.gene_TSS_map} | cut -f 1)
        REGION_START=$(grep {wildcards.gene} {input.gene_TSS_map} | cut -f 2)
        REGION_STOP=$(grep {wildcards.gene} {input.gene_TSS_map} | cut -f 3)
        fastQTL --vcf {input.vcf} --bed {input.pheno} \
	       --include-samples {input.sample_input} \
           --include-phenotypes {input.gene_input} \
	       --include-covariates {input.covariates} \
           --out {output} \
           --window {WINDOW} \
           --region chr$CHR:$REGION_START-$REGION_STOP \
           --permute {PERM}
        """

rule estimate_eQTLs:
    input:
        vcf=DATA_DIR + VCF,
        pheno=ancient(rules.normalize_expression.output.bed),
        sample_input=DATA_DIR + "fastqtl_sample_input/estimation/{anc}/{gene}.txt",
        gene_input=rules.prep_pheno_file.output,
        covariates=rules.combine_covariates.output,
        gene_TSS_map=rules.make_gene_list.output
    output:
        DATA_DIR + "fastqtl_output/estimation/{anc}/{gene}.txt"
    singularity:
        FASTQTL_DOCKER 
    params:
        output_dir=DATA_DIR + "fastqtl_output/"
    shell:
        """
        mkdir -p {params.output_dir}/estimation/{wildcards.anc}
        CHR=$(grep {wildcards.gene} {input.gene_TSS_map} | cut -f 1)
        REGION_START=$(grep {wildcards.gene} {input.gene_TSS_map} | cut -f 2)
        REGION_STOP=$(grep {wildcards.gene} {input.gene_TSS_map} | cut -f 3)
        fastQTL --vcf {input.vcf} --bed {input.pheno} \
	       --include-samples {input.sample_input} \
           --include-phenotypes {input.gene_input} \
	       --include-covariates {input.covariates} \
           --out {output} \
           --window {WINDOW} \
           --region chr$CHR:$REGION_START-$REGION_STOP
        """

##################################### RUN #####################################

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
        dir=DATA_DIR + "fastqtl_output/{type}/{anc}"
    shell:
        """
        cat {params.dir} > {output}
        """
