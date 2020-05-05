# Variables and filepaths stored here
include: "scripts/snakemake_variables.py"

shell.executable("/usr/bin/bash")
shell.prefix("source ~/.bashrc; ")

# TODO use glob + expand to properly specify input of make_anc_bed
 
rule all:
    input:
        DATA_DIR + "fastqtl_output/randsamp_merged_ascertainment_Eur.txt",
        DATA_DIR + "fastqtl_output/randsamp_merged_estimation_Eur.txt",
        DATA_DIR + "fastqtl_output/randsamp_merged_validation_Eur.txt",
        DATA_DIR + "fastqtl_output/randsamp_merged_validation_Afr.txt",
        DATA_DIR + "fastqtl_output/randsamp_merged_estimation_het.txt",
        DATA_DIR + "fastqtl_output/randsamp_merged_estimation_Afr.txt",
        DATA_DIR + "fastqtl_output/hits_ascertainment_Eur.txt",
        DATA_DIR + "fastqtl_output/hits_estimation_Eur.txt",
        DATA_DIR + "fastqtl_output/hits_validation_Eur.txt",
        DATA_DIR + "fastqtl_output/hits_validation_Afr.txt",
        DATA_DIR + "fastqtl_output/hits_estimation_het.txt",
        DATA_DIR + "fastqtl_output/hits_estimation_Afr.txt"
#         DATA_DIR + "fastqtl_output/merged_ascertainment_Eur.txt",
#         DATA_DIR + "fastqtl_output/merged_estimation_Afr.txt",
#         DATA_DIR + "fastqtl_output/merged_estimation_het.txt",
#         DATA_DIR + "fastqtl_output/merged_estimation_Eur.txt",
#         DATA_DIR + "fastqtl_output/merged_validation_Eur.txt"
#         DATA_DIR + "fastqtl_output/ascertainment/Eur/ENSG00000000419.12.txt"
#         DATA_DIR + "hits/hits_estimation_Eur.txt",
#         DATA_DIR + "hits/hits_estimation_het.txt",
#         DATA_DIR + "hits/hits_estimation_Afr.txt"
#         expand(DATA_DIR + "mesa.residual_expression.{anc}.txt", anc=["Afr", "Eur"])
#         DATA_DIR + "fastqtl_output/estimation/het/ENSG00000000419.12.txt"

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
        DATA_DIR + "sample_participant_lookup_{anc}.txt"
    shell:
        """
        conda activate py36
        python {SAMPLE_SELECTION_SCRIPT} \
            --sample_data {input.sample_data} \
            --ind_data {input.ind_data} \
            --ancestry {wildcards.anc} \
            --genotyped_individuals {input.geno_ID} \
            --rnaseq_individuals {input.exp_ID} \
            --output {output}
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
        bed=DATA_DIR + "mesa.{anc}.expression.bed.gz",
        idx=DATA_DIR + "mesa.{anc}.expression.bed.gz.tbi"
    params:
        prefix=DATA_DIR + "mesa.{anc}"
    shell:
        """
        conda activate norm-exp-env
        {NORM_SCRIPT} {input.tpm} {input.counts} {input.anno} \
            {input.sample_map} {input.chr_list} {params.prefix} \
            --tpm_threshold 0.1 \
            --count_threshold 6 \
            --sample_frac_threshold 0.2 
        conda deactivate
        """

rule make_gene_list:
    input:
        expand(rules.normalize_expression.output.bed, anc=["Afr", "Eur"])
    output:
        DATA_DIR + "all_genes_TSS.txt"
    shell:
        """
        conda activate py36
        python {MERGE_GENES} --bed {input} --out {output}
        conda deactivate
        """

############################# GENERATE COVARIATES #############################

rule indiv_filter:
    input:
        vcf=ancient(DATA_DIR + VCF),
        samples=rules.select_samples.output
    output:
        vcf=DATA_DIR + "mesa.{anc}.vcf.gz",
        idx=DATA_DIR + "mesa.{anc}.vcf.gz.tbi"
    shell:
        """
        conda activate bcftools-env
        SAMP=$(tail -n +2 {input.samples} | cut -f2 | tr '\n' ',' | sed 's/,$//')
        bcftools view -s $SAMP --force-samples -Ou {input.vcf} | \
            bcftools view --genotype ^miss --phased -Oz -o {output.vcf}
        conda deactivate
        conda activate tabix-env
        tabix -p vcf {output.vcf}
        conda deactivate
        """

rule MAF_filter:
    input:
        rules.indiv_filter.output.vcf
    output:
        DATA_DIR + "maf05.mesa.{anc}.vcf.gz"
    shell:
        """
        conda activate bcftools-env
        bcftools view --min-af .05 -Oz -o {output} {input}
        conda deactivate
        """

rule find_indep_snps:
    input:
        vcf=rules.MAF_filter.output
    output:
        keep=DATA_DIR + "maf05.mesa.{anc}.prune.in",
        exclude=temp(expand(DATA_DIR + "maf05.mesa.{{anc}}.{ext}", ext=["log", "nosex", "prune.out"]))
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
        DATA_DIR + "prune.maf05.mesa.{anc}.vcf.gz"
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
        DATA_DIR + "mesa.{anc}.sample_covariates.txt"
    shell:
        """
        conda activate pystats
        python {COV_PREP_SCRIPT} \
            --samples {input.samples} \
            --indiv_metadata {input.indiv_data} \
            --sample_metadata {input.sample_data} \
            --vcf {input.vcf} \
            --exp {input.exp} \
            --out {output}
        conda deactivate
        """

rule generate_group_PC:
    input:
        afr=expand(rules.select_samples.output, anc="Afr"),
        eur=expand(rules.select_samples.output, anc="Eur"),
        exp=expand(rules.normalize_expression.output.bed, anc=["Afr", "Eur"])
    output:
        expand(DATA_DIR + "mesa.residual_expression.{anc}.bed.gz", anc=["Afr", "Eur"])
    params:
        DATA_DIR + "mesa.residual_expression"
    shell:
        """
        conda activate pystats
        python {GRP_PREP_SCRIPT} \
            --afr {input.afr} \
            --eur {input.eur} \
            --exp {input.exp} \
            --out {params}
        conda deactivate
        """

rule exclude_all_PC:
    input:
        rules.prep_covariates.output
    output:
        DATA_DIR + "mesa.{anc}.sample_covariates.no_PC.txt"
    shell:
        """
        cut -f 1-5 {input} > {output}
        """

rule exclude_exp_PC:
    input:
        rules.prep_covariates.output
    output:
        DATA_DIR + "mesa.{anc}.sample_covariates.no_exp_PC.txt"
    shell:
        """
        cut -f 1-5,16-30 {input} > {output}
        """

rule extract_exp_PC:
    input:
        rules.prep_covariates.output
    output:
        DATA_DIR + "mesa.{anc}.sample_covariates.exp_PC.txt"
    shell:
        """
        cut -f 1,6-15 {input} > {output}
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
	afr_samples=expand(rules.select_samples.output, anc=["Afr"]),
	eur_samples=expand(rules.select_samples.output, anc=["Eur"]),
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
        mkdir -p {params.output_dir}/validation/Eur
        mkdir -p {params.output_dir}/validation/Afr
        conda activate py36
        python {PARTITION_SCRIPT} \
            --intersect {input.intersect} \
            --afr_samples {input.afr_samples} \
            --eur_samples {input.eur_samples} \
            --genes {input.genes} \
            --afr_validation_run \
            --out {params.output_dir}
        conda deactivate
        """

# checkpoint add_validation_samp:
#     input:
#         intersect=rules.intersect_tracts_genes.output,
# 	afr_samples=expand(rules.select_samples.output, anc=["Afr"]),
# 	eur_samples=expand(rules.select_samples.output, anc=["Eur"]),
#         genes=rules.make_gene_list.output
#     output:
#         directory(DATA_DIR + "secondary_validation")
#     params:
#         output_dir=DATA_DIR + "secondary_validation"
#     shell:
#         """
#         mkdir -p {params.output_dir}/validation/Afr
#         conda activate py36
#         python {PARTITION_SCRIPT} \
#             --intersect {input.intersect} \
#             --afr_samples {input.afr_samples} \
#             --eur_samples {input.eur_samples} \
#             --genes {input.genes} \
#             --afr_validation_run \
#             --out {params.output_dir}
#         conda deactivate
#         """
# 
# rule move_validation_files:
#     input:
#         DATA_DIR + "secondary_validation/validation/Afr/{gene}.txt"
#     output:
#         DATA_DIR + "fastqtl_sample_input/validation/Afr/{gene}.txt"
#     params:
#         output_dir=DATA_DIR + "fastqtl_sample_input/validation/Afr"
#     shell:
#         """
#         mkdir -p {params.output_dir}
#         cp {input} {output}
#         """

################################## CALL EQTLS ##################################

rule subset_geno:
    input:
        vcf=rules.indiv_filter.output,
        gene_window=rules.make_genes_bed.output
    output:
        DATA_DIR + "fastqtl_geno_input/{anc}/{gene}.vcf.gz"
    params:
        out_dir=DATA_DIR + "fastqtl_geno_input/{anc}"
    shell:
        """
        mkdir -p {params.out_dir}
        CHR=$(grep {wildcards.gene} {input.gene_window} | cut -f 1)
        REGION_START=$(grep {wildcards.gene} {input.gene_window} | cut -f 2)
        REGION_STOP=$(grep {wildcards.gene} {input.gene_window} | cut -f 3)
        
        conda activate bcftools-env
        bcftools view -r chr$CHR:$REGION_START-$REGION_STOP -Oz -o {output} {input.vcf}
        conda deactivate
        """

rule subset_pheno:
    input:
        bed=DATA_DIR + "mesa.{anc}.expression.bed.gz"
    output:
        DATA_DIR + "fastqtl_pheno_input/{anc}/{gene}.txt"
    params:
        out_dir=DATA_DIR + "fastqtl_pheno_input/{anc}"
    shell:
        """
        mkdir -p {params.out_dir}
        cat <(zcat {input.bed} | head -n 1) <(zgrep -m 1 {wildcards.gene} {input.bed}) > {output}
        """

rule call_eQTLs:
    input:
        geno_input=lambda wildcards: expand(rules.subset_geno.output, anc="Eur", \
            gene=wildcards.gene) if wildcards.anc == "Eur" else \
            expand(rules.subset_geno.output, anc="Afr", gene=wildcards.gene),
        pheno_input=lambda wildcards: expand(rules.subset_pheno.output, anc="Eur", \
            gene=wildcards.gene) if wildcards.anc == "Eur" else \
            expand(rules.subset_pheno.output, anc="Afr", gene=wildcards.gene),
        sample_input=DATA_DIR + "fastqtl_sample_input/{type}/{anc}/{gene}.txt"
    output:
        DATA_DIR + "fastqtl_output/{type}/{anc}/{gene}.txt"
    params:
        output_dir=DATA_DIR + "fastqtl_output/{type}/{anc}"
    shell:
        """
        mkdir -p {params.output_dir}
        conda activate pystats
        python {EQTL_SCRIPT} \
            --phenotypes {input.pheno_input} \
            --genotypes {input.geno_input} \
            --samples {input.sample_input} \
            --type {wildcards.type} \
            --out {output}
        conda deactivate
        """

def merge_input(wildcards):
    checkpoint_output = checkpoints.partition_samples.get(**wildcards).output[0]
    # checkpoint_output = checkpoints.add_validation_samp.get(**wildcards).output[0]
    return expand(DATA_DIR + "fastqtl_output/{type}/{anc}/{gene}.txt",
                  type=wildcards.type, anc=wildcards.anc,
                  gene=glob_wildcards(os.path.join(checkpoint_output, wildcards.type, 
                                      wildcards.anc, "{gene}.txt")).gene)

rule merge_output:
    input:
        merge_input
    output:
        DATA_DIR + "fastqtl_output/merged_{type}_{anc}.txt"
    params:
        dir=DATA_DIR + "fastqtl_output/{type}/{anc}/"
    shell:
        """
        conda activate py36
        python {MERGE_RES} --dir {params.dir} --out {output}
        conda deactivate 
        """

rule sample_merged:
    input:
        DATA_DIR + "fastqtl_output/merged_{type}_{anc}.txt"
    output:
        DATA_DIR + "fastqtl_output/randsamp_merged_{type}_{anc}.txt"
    shell:
        """
        cat <(head -n 1 {input}) <(tail -n +2 {input} | shuf -n 100000) > {output}
        """

rule identify_hits:
    input:
        asc=DATA_DIR + "fastqtl_output/merged_ascertainment_Eur.txt", 
        eur=DATA_DIR + "fastqtl_output/merged_estimation_Eur.txt",
        val=DATA_DIR + "fastqtl_output/merged_validation_Eur.txt",
        val_afr=DATA_DIR + "fastqtl_output/merged_validation_Afr.txt",
        het=DATA_DIR + "fastqtl_output/merged_estimation_het.txt",
        afr=DATA_DIR + "fastqtl_output/merged_estimation_Afr.txt"
    output:
        DATA_DIR + "fastqtl_output/hits_ascertainment_Eur.txt",
        DATA_DIR + "fastqtl_output/hits_estimation_Eur.txt",
        DATA_DIR + "fastqtl_output/hits_validation_Eur.txt",
        DATA_DIR + "fastqtl_output/hits_validation_Afr.txt",
        DATA_DIR + "fastqtl_output/hits_estimation_het.txt",
        DATA_DIR + "fastqtl_output/hits_estimation_Afr.txt"
    params:
        DATA_DIR + "fastqtl_output/"
    shell:
        """
        mkdir -p {params}
        module reset
        module load R
        Rscript --vanilla {HITS_SCRIPT} --asc {input.asc} --eur {input.eur} \
            --het {input.het} --afr {input.afr} --val {input.val} --val_afr {input.val_afr} \
            --out {params}
        conda deactivate
        """

############################## SIMULATIONS #################################

rule analyze_expression_sd:
    input:
        afr=DATA_DIR + "mesa.Afr.expression.bed.gz",
        eur=DATA_DIR + "mesa.Eur.expression.bed.gz"
    params:
        samp=DATA_DIR + "fastqtl_sample_input/estimation/"
    output:
        DATA_DIR + "expression_sd.txt"
    shell:
        """
        module reset
        module load R
        Rscript --vanilla {SD_SCRIPT} --afr {input.afr} --eur {input.eur} \
                --samp_prefix {params.samp} --out {output}
        """

rule perform_simulation:
    input:
        sd=rules.analyze_expression_sd.output,
        afr=DATA_DIR + "fastqtl_output/hits_estimation_Afr.txt",
        eur=DATA_DIR + "fastqtl_output/hits_estimation_Eur.txt"
    params:
        out_dir=DATA_DIR + "sims"
    output:
        DATA_DIR + "sims/scaling_{afr}_{eur}/sim_{sim}.txt"
    shell:
        """
        mkdir -p {params.out_dir}
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
        expand(DATA_DIR + "sims/scaling_{{afr}}_{{eur}}/sim_{num}.txt", num=range(SIM_ITER))
    output:
        DATA_DIR + "sims/scaling_{afr}_{eur}/merged_results.txt"
    params:
        dir=DATA_DIR + "sims/scaling_{afr}_{eur}/"
    shell:
        """
        module reset
        module load R
        Rscript --vanilla {MERGE_SIMS} --dir {params.dir} \
                --out {output}
        """
