# Variables and filepaths stored here
include: "scripts/snakemake_variables.py"

shell.executable("/usr/bin/bash")
shell.prefix("source ~/.bashrc; ")

rule all:
    input:
        ### For eQTL calling
        # DATA_DIR + "fastqtl_output/randsamp_merged_ascertainment_Eur.txt",
        # DATA_DIR + "fastqtl_output/randsamp_merged_reestimation_primary_Eur.txt",
        # DATA_DIR + "fastqtl_output/randsamp_merged_reestimation_validation_Eur.txt",
        # DATA_DIR + "fastqtl_output/randsamp_merged_reestimation_validation_Afr.txt",
        # DATA_DIR + "fastqtl_output/randsamp_merged_reestimation_primary_het.txt",
        # DATA_DIR + "fastqtl_output/randsamp_merged_reestimation_primary_Afr.txt",
        # DATA_DIR + "fastqtl_output/hits_ascertainment_Eur.txt",
        # DATA_DIR + "fastqtl_output/hits_reestimation_validation_Eur.txt",
        # DATA_DIR + "fastqtl_output/hits_reestimation_validation_Afr.txt",
        # DATA_DIR + "fastqtl_output/hits_reestimation_primary_Eur.txt",
        # DATA_DIR + "fastqtl_output/hits_reestimation_primary_het.txt",
        # DATA_DIR + "fastqtl_output/hits_reestimation_primary_Afr.txt"
        #
        ### For likelihood model stuff:
        # DATA_DIR + "MLE/merged_data.txt",
        # expand(DATA_DIR + "MLE/{cov}/optimize_{option}_{param}.txt", cov=["global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1"], option=["neg_control", "drop_asc"], param=["delta", "betas"]),
        # expand(DATA_DIR + "MLE/{cov}/results_bootstrap_{option}_delta.summary.txt", cov=["global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1"], option=['neg_control', 'drop_asc']),
        # expand(DATA_DIR + "MLE/{cov}/likelihood_vs_delta_{option}.txt", option=["neg_control", "drop_asc"], cov=["global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1"])
        # expand(DATA_DIR + "MLE/simulated_optimization/delta{delta}.afr_error{afr_error}.eur_error{eur_error}.corr{corr}.n_snp{n_snp}.n_idv{n_idv}/bias_variance_report.txt", delta=[0, 0.5], eur_error=[0.12], afr_error=[0.13, 0.17, 0.30], corr=[0.8], n_snp=[3200], n_idv=["200_150"]),
        # expand(DATA_DIR + "MLE/simulated_optimization/delta{delta}.afr_error{afr_error}.eur_error{eur_error}.corr{corr}.n_snp{n_snp}.n_idv{n_idv}/bias_variance_report.txt", delta=[0, 0.5], eur_error=[0.12], afr_error=[0.17], corr=[0.5, 0.8, 0.9], n_snp=[3200], n_idv=["200_150"]),
        # expand(DATA_DIR + "MLE/simulated_optimization/delta{delta}.afr_error{afr_error}.eur_error{eur_error}.corr{corr}.n_snp{n_snp}.n_idv{n_idv}/bias_variance_report.txt", delta=[0, 0.5], eur_error=[0.12], afr_error=[0.17], corr=[0.8], n_snp=[1000, 3200, 6000], n_idv=["200_150"]),
        # expand(DATA_DIR + "MLE/simulated_optimization/delta{delta}.afr_error{afr_error}.eur_error{eur_error}.corr{corr}.n_snp{n_snp}.n_idv{n_idv}/bias_variance_report.txt", delta=[0, 0.5], eur_error=[0.12], afr_error=[0.17], corr=[0.8], n_snp=[3200], n_idv=["200_150", "350_150", "500_300"])
        #
        ### For covariate/PC stuff:
        DATA_DIR + "correlation.expression_regressed.global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1.txt"

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

############################# NORMALIZE EXPRESSION #############################

rule select_samples:
    input:
        sample_data=DATA_DIR + RNASEQ_METADATA,
        exclusion_list=EXCLUSION_FILE
    output:
        expand(DATA_DIR + "select_samples_{anc}.txt", anc=["Afr", "Eur"])
    params:
        prefix=DATA_DIR + "select_samples"
    shell:
        """
        conda activate py36
        python {SAMPLE_SELECTION_SCRIPT} \
            --sample_data {input.sample_data} \
            --exclusion_list {input.exclusion_list} \
            --output {params.prefix}
        conda deactivate
        """

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
        sample_map=DATA_DIR + "select_samples_{anc}.txt",
        chr_list=rules.make_chr_list.output,
        flag=DATA_DIR + "vcf_filter_{anc}.done"
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

rule snp_filter:
    input:
        vcf=ancient(DATA_DIR + VCF)
    output:
        vcf=DATA_DIR + "mesa.filt.bcf.gz",
        idx=DATA_DIR + "mesa.filt.bcf.gz.csi"
    shell:
        """
        conda activate bcftools-env
        bcftools view -m2 -M2 -v snps -Ob --genotype ^miss --phased -o {output.vcf} {input.vcf}
        bcftools index -c {output.vcf}
        conda deactivate
        """
    
rule indiv_filter:
    input:
        vcf=rules.snp_filter.output.vcf,
        samples=DATA_DIR + "select_samples_{anc}.txt"
    output:
        vcf=DATA_DIR + "mesa.{anc,[A-Za-z]+}.bcf.gz",
        idx=DATA_DIR + "mesa.{anc,[A-Za-z]+}.bcf.gz.csi",
        flag=DATA_DIR + "vcf_filter_{anc}.done"
    shell:
        """
        conda activate bcftools-env
        SAMP=$(tail -n +2 {input.samples} | cut -f2 | tr '\n' ',' | sed 's/,$//')
        bcftools view -m2 -M2 -v snps -Ou --genotype ^miss --phased {input.vcf} | \
            bcftools view -s $SAMP --force-samples -Ob -o {output.vcf}
        bcftools index -c {output.vcf}
        # Hodge-podge way of updating sample list based on individuals actually present in VCF
        cat <(head -n 1 {input.samples}) <(grep -f <(bcftools query -l {input.vcf}) {input.samples}) > tmp_{wildcards.anc}.txt
        mv tmp_{wildcards.anc}.txt {input.samples}
        touch {output.flag}
        conda deactivate
        """

rule MAF_filter:
    input:
        rules.indiv_filter.output.vcf
    output:
        vcf=DATA_DIR + "mesa.{anc}.maf05.vcf.gz",
        idx=DATA_DIR + "mesa.{anc}.maf05.vcf.gz.csi"
    shell:
        """
        conda activate bcftools-env
        bcftools view --min-af .05 -Oz -o {output.vcf} {input}
        bcftools index -c {output.vcf}
        conda deactivate
        """

rule find_indep_snps:
    input:
        vcf=rules.MAF_filter.output.vcf
    output:
        keep=temp(DATA_DIR + "mesa.{anc}.maf05.prune.in"),
        exclude=temp(expand(DATA_DIR + "mesa.{{anc}}.maf05.{ext}", ext=["log", "nosex", "prune.out"]))
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
        vcf=rules.MAF_filter.output.vcf,
        snps=rules.find_indep_snps.output.keep
    output:
        DATA_DIR + "mesa.{anc}.prune.maf05.vcf.gz"
    shell:
        """
        conda activate bcftools-env
        bcftools view -i 'ID=@{input.snps}' -Oz -o {output} {input.vcf}
        conda deactivate
        """

rule PCA:
    input:
        expand(DATA_DIR + "mesa.{anc}.{{filetype}}.gz", anc=["Afr", "Eur"])
    output:
        pc=DATA_DIR + "mesa.{filetype}.principal_components.txt",
        eig=DATA_DIR + "mesa.{filetype}.eigenvalues.txt"
    params:
        out=DATA_DIR + "mesa.{filetype}"
    shell:
        """
        conda activate pystats
        python {PCA_SCRIPT} \
            --data {input} \
            --n 10 \
            --filetype {wildcards.filetype} \
            --out {params.out}
        conda deactivate
        """

rule prep_covariates:
    input:
        samples=expand(DATA_DIR + "select_samples_{anc}.txt", anc=["Afr", "Eur"]),
        ganc=DATA_DIR + "combined_global_anc_frac.txt",
        pc=expand(rules.PCA.output.pc, filetype=["prune.maf05.vcf", "expression.bed"]),
        sample_data=DATA_DIR + RNASEQ_METADATA,
        flag=expand(DATA_DIR + "vcf_filter_{anc}.done", anc=["Afr", "Eur"])
    output:
        DATA_DIR + "mesa.sample_covariates.txt"
    shell:
        """
        conda activate pystats
        python {COV_PREP_SCRIPT} \
            --samples {input.samples} \
            --sample_metadata {input.sample_data} \
            --global_anc {input.ganc} \
            --pca {input.pc} \
            --out {output}
        conda deactivate
        """

rule extract_linreg_covariates:
    input:
        DATA_DIR + "mesa.sample_covariates.txt"
    output:
        DATA_DIR + "mesa.sample_covariates_for_regression.txt"
    shell:
        """
        cut -f 2,3,4,10 {input} > {output}
        """

############################## PARTITION SAMPLES ##############################

rule make_anc_bed:
    input:
        expand(BED_DIR + "chr{chr}/job_complete.txt", chr=CHROMS)
    output:
        DATA_DIR + "anc_tracts.bed"
    params:
        dir=BED_DIR
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
	afr_samples=expand(DATA_DIR + "select_samples_{anc}.txt", anc=["Afr"]),
	eur_samples=expand(DATA_DIR + "select_samples_{anc}.txt", anc=["Eur"]),
        genes=rules.make_gene_list.output,
        flag=expand(DATA_DIR + "vcf_filter_{anc}.done", anc=["Afr", "Eur"])
    output:
        directory(DATA_DIR + "fastqtl_sample_input")
    params:
        output_dir=DATA_DIR + "fastqtl_sample_input"
    shell:
        """
        mkdir -p {params.output_dir}/ascertainment/Eur
        mkdir -p {params.output_dir}/reestimation_primary/Afr
        mkdir -p {params.output_dir}/reestimation_primary/het
        mkdir -p {params.output_dir}/reestimation_primary/Eur
        mkdir -p {params.output_dir}/reestimation_validation/Eur
        mkdir -p {params.output_dir}/reestimation_validation/Afr
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

############################ LIN REG EQTL CALLING ############################

rule subset_geno:
    input:
        vcf=DATA_DIR + "mesa.{anc}.maf05.vcf.gz",
        idx=DATA_DIR + "mesa.{anc}.maf05.vcf.gz.csi",
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
        sample_input=DATA_DIR + "fastqtl_sample_input/{type}/{anc}/{gene}.txt",
        covariates=DATA_DIR + "mesa.sample_covariates_for_regression.txt"
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
            --covariates {input.covariates} \
            --samples {input.sample_input} \
            --type {wildcards.type} \
            --out {output}
        conda deactivate
        """

def merge_asc_input(wildcards):
    checkpoint_output = checkpoints.partition_samples.get(**wildcards).output[0]
    return expand(DATA_DIR + "fastqtl_output/ascertainment/Eur/{gene}.txt",
                  gene=glob_wildcards(os.path.join(checkpoint_output, "ascertainment", 
                                      "Eur", "{gene}.txt")).gene)

def merge_est_input(wildcards):
    checkpoint_output = checkpoints.identify_hits.get(**wildcards).output[0]
    sig_genes = list(pd.read_csv(checkpoint_output, sep='\t')["gene"])
    return expand(DATA_DIR + "fastqtl_output/reestimation_{type}/{anc}/{gene}.txt",
                  type=wildcards.type, anc=wildcards.anc,
                  gene=sig_genes)

rule merge_ascertainment:
    input:
        merge_asc_input
    output:
        DATA_DIR + "fastqtl_output/merged_ascertainment_Eur.txt"
    params:
        dir=DATA_DIR + "fastqtl_output/ascertainment/Eur/"
    shell:
        """
        conda activate py36
        python {MERGE_RES} --dir {params.dir} --out {output}
        conda deactivate 
        """

rule merge_reestimation:
    input:
        merge_est_input
    output:
        DATA_DIR + "fastqtl_output/merged_reestimation_{type}_{anc}.txt"
    params:
        dir=DATA_DIR + "fastqtl_output/reestimation_{type}/{anc}/"
    shell:
        """
        conda activate py36
        python {MERGE_RES} --dir {params.dir} --out {output}
        conda deactivate 
        """

checkpoint identify_hits:
    input:
        asc=DATA_DIR + "fastqtl_output/merged_ascertainment_Eur.txt" 
    output:
        DATA_DIR + "fastqtl_output/hits_ascertainment_Eur.txt"
    params:
        DATA_DIR + "fastqtl_output/"
    shell:
        """
        mkdir -p {params}
        conda activate py36
        python {ID_HITS_SCRIPT} --ascertainment \
            --merged {input.asc} \
            --out {output} \
            --fdr {FDR}
        conda deactivate
        """

rule extract_hits:
    input:
        merged=DATA_DIR + "fastqtl_output/merged_reestimation_{type}_{anc}.txt",
        hits=DATA_DIR + "fastqtl_output/hits_ascertainment_Eur.txt"
    output:
        DATA_DIR + "fastqtl_output/hits_reestimation_{type}_{anc}.txt"
    params:
        DATA_DIR + "fastqtl_output/"
    shell:
        """
        mkdir -p {params}
        conda activate py36
        python {ID_HITS_SCRIPT} --reestimation \
            --merged {input.merged} \
            --hits {input.hits} \
            --out {output}
        conda deactivate
        """

############################# LIKELIHOOD MODEL ################################

def hit_genes():
    afr_df = pd.read_csv("data/fastqtl_output/hits_reestimation_primary_Afr.txt", sep='\t')
    eur_df = pd.read_csv("data/fastqtl_output/hits_reestimation_primary_Eur.txt", sep='\t')
    afr_genes = set(afr_df["gene"])
    eur_genes = set(eur_df["gene"])
    genes = list(afr_genes.intersection(eur_genes))
    return(genes)

rule merge_data:
    input:
        afr_geno=expand(rules.subset_geno.output, gene=hit_genes(), anc="Afr"),
        afr_pheno=expand(rules.subset_pheno.output, gene=hit_genes(), anc="Afr"),
        eur_geno=expand(rules.subset_geno.output, gene=hit_genes(), anc="Eur"),
        eur_pheno=expand(rules.subset_pheno.output, gene=hit_genes(), anc="Eur"),
        tracts=rules.make_anc_bed.output,
        cov=rules.prep_covariates.output,
        afr_hits="data/fastqtl_output/hits_reestimation_primary_Afr.txt",
        eur_hits="data/fastqtl_output/hits_reestimation_primary_Eur.txt"
    output:
        merged=DATA_DIR + "MLE/merged_data.txt"
    params:
        out_dir=DATA_DIR + "MLE/"
    shell:
        """
        mkdir -p {params.out_dir}
        conda activate pystats
        python scripts/merge_eQTL_data.py --tracts {input.tracts} \
            --afr_hits {input.afr_hits} \
            --eur_hits {input.eur_hits} \
            --covariates {input.cov} \
            --out {output.merged}
        conda deactivate
        """

rule merge_swapped_data:
    input:
        afr_geno=expand(rules.subset_geno.output, gene=hit_genes(), anc="Afr"),
        afr_pheno=expand(rules.subset_pheno.output, gene=hit_genes(), anc="Afr"),
        eur_geno=expand(rules.subset_geno.output, gene=hit_genes(), anc="Eur"),
        eur_pheno=expand(rules.subset_pheno.output, gene=hit_genes(), anc="Eur"),
        tracts=rules.make_anc_bed.output,
        afr_hits="results/QTL_calling/v3.1/hits/hits_estimation_Afr.txt",
        eur_hits="results/QTL_calling/v3.1/hits/hits_estimation_Eur.txt"
    output:
        merged=DATA_DIR + "MLE/merged_data_n5_swapped.txt"
    params:
        out_dir=DATA_DIR + "MLE/"
    shell:
        """
        mkdir -p {params.out_dir}
        conda activate pystats
        python scripts/merge_eQTL_data.py --tracts {input.tracts} \
            --afr_hits {input.afr_hits} \
            --eur_hits {input.eur_hits} \
            --n_genes 5 \
            --swap_ref_alt \
            --out {output.merged}
        conda deactivate
        """

rule optimize:
    input:
        DATA_DIR + "MLE/merged_data.txt"
    output:
        delta=DATA_DIR + "MLE/{cov}/optimize_{option}_delta.txt",
        betas=DATA_DIR + "MLE/{cov}/optimize_{option}_betas.txt"
    params:
        cov=lambda wildcards: wildcards.cov.split('.'),
        dir=DATA_DIR + "MLE/{cov}"
    shell:
        """
        mkdir -p {params.dir}
        conda activate pystats
        python scripts/iterative_parameter_optimization.py --merged {input} \
            --max_iter {MAX_ITER} \
            --covariates {params.cov} \
            --option {wildcards.option} \
            --betas_out {output.betas} \
            --delta_out {output.delta} 
        conda deactivate
        """

rule perform_simulation:
    input:
        afr="data/fastqtl_output/hits_reestimation_primary_Afr.txt",
        eur="data/fastqtl_output/hits_reestimation_primary_Eur.txt"
    output:
        DATA_DIR + "MLE/simulated_data/delta{delta}.afr_error{afr_error}.eur_error{eur_error}.corr{corr}.n_snp{n_snp}.n_idv{n_idv}/{idx}.txt"
    params:
        dir=DATA_DIR + "MLE/simulated_data/delta{delta}.afr_error{afr_error}.eur_error{eur_error}.corr{corr}.n_snp{n_snp}.n_idv{n_idv}"
    shell:
        """
        mkdir -p {params.dir}
        conda activate pystats
        python {SIM_SCRIPT} \
            --correlation {wildcards.corr} \
            --afr_error {wildcards.afr_error} \
            --eur_error {wildcards.eur_error} \
            --n_snp {wildcards.n_snp} \
            --n_idv {wildcards.n_idv} \
            --afr_hits {input.afr} \
            --eur_hits {input.eur} \
            --delta {wildcards.delta} \
            --out {output}
        conda deactivate
        """

rule optimize_simulated_data:
    input:
        DATA_DIR + "MLE/simulated_data/delta{delta}.afr_error{afr_error}.eur_error{eur_error}.corr{corr}.n_snp{n_snp}.n_idv{n_idv}/{idx}.txt"
    output:
        delta=temp(DATA_DIR + "MLE/simulated_optimization/delta{delta}.afr_error{afr_error}.eur_error{eur_error}.corr{corr}.n_snp{n_snp}.n_idv{n_idv}/delta_{idx}.txt"),
        betas=temp(DATA_DIR + "MLE/simulated_optimization/delta{delta}.afr_error{afr_error}.eur_error{eur_error}.corr{corr}.n_snp{n_snp}.n_idv{n_idv}/betas_{idx}.txt")
    params:
        dir=DATA_DIR + "MLE/simulated_optimization/delta{delta}.afr_error{afr_error}.eur_error{eur_error}.corr{corr}.n_snp{n_snp}.n_idv{n_idv}"
    shell:
        """
        mkdir -p {params.dir}
        conda activate pystats
        python scripts/iterative_parameter_optimization.py --merged {input} \
            --max_iter {MAX_ITER} \
            --covariates \
            --betas_out {output.betas} \
            --delta_out {output.delta} 
        conda deactivate
        """

rule parse_simulation_results:
    input:
        expand(DATA_DIR + "MLE/simulated_optimization/delta{{delta}}.afr_error{{afr_error}}.eur_error{{eur_error}}.corr{{corr}}.n_snp{{n_snp}}.n_idv{{n_idv}}/delta_{idx}.txt", idx=[str(i).zfill(3) for i in range(100)])
    output:
        DATA_DIR + "MLE/simulated_optimization/delta{delta}.afr_error{afr_error}.eur_error{eur_error}.corr{corr}.n_snp{n_snp}.n_idv{n_idv}/bias_variance_report.txt"
    shell:
        """
        conda activate pystats
        python scripts/compute_bias_variance.py --delta_files {input} \
            --simulated_delta {wildcards.delta} \
            --out {output}
        conda deactivate
        """

rule bootstrap:
    input:
        DATA_DIR + "MLE/merged_data.txt"
    output:
        delta=temp(DATA_DIR + "MLE/{cov}/bootstrap_{option}_delta_{idx}.txt"),
        betas=temp(DATA_DIR + "MLE/{cov}/bootstrap_{option}_betas_{idx}.txt")
    params:
        cov=lambda wildcards: wildcards.cov.split('.'),
        dir=DATA_DIR + "MLE/{cov}"
    shell:
        """
        mkdir -p {params.dir}
        conda activate pystats
        python scripts/iterative_parameter_optimization.py --merged {input} \
            --max_iter {MAX_ITER} \
            --bootstrap \
            --covariates {params.cov} \
            --option {wildcards.option} \
            --betas_out {output.betas} \
            --delta_out {output.delta} 
        conda deactivate
        """

rule parse_bootstrap:
    input:
        expand(DATA_DIR + "MLE/{{cov}}/bootstrap_{{option}}_delta_{idx}.txt", idx=[str(i).zfill(3) for i in range(1000)])
    output:
        DATA_DIR + "MLE/{cov}/results_bootstrap_{option}_delta.summary.txt",
        DATA_DIR + "MLE/{cov}/results_bootstrap_{option}_delta.all_values.txt"
    params:
        prefix=DATA_DIR + "MLE/{cov}/results_bootstrap_{option}_delta."
    shell:
        """
        conda activate py36
        python scripts/parse_bootstrap.py --bootstrap_files {input} --out {params.prefix}
        conda deactivate
        """

rule likelihood:
    input:
        merged=DATA_DIR + "MLE/merged_data.txt"
    output:
        DATA_DIR + "MLE/{cov}/likelihood_vs_delta_{option}.txt"
    params:
        cov=lambda wildcards: wildcards.cov.split('.'),
        dir=DATA_DIR + "MLE/{cov}"
    shell:
        """
        mkdir -p {params.dir}
        conda activate pystats
        python scripts/compute_likelihood.py --merged {input.merged} \
            --option {wildcards.option} \
            --covariates {params.cov} \
            --out {output}
        conda deactivate
        """

############################## MISC. ANALYSES #################################

rule sample_merged:
    input:
        DATA_DIR + "fastqtl_output/merged_{type}_{anc}.txt"
    output:
        DATA_DIR + "fastqtl_output/randsamp_merged_{type}_{anc}.txt"
    shell:
        """
        cat <(head -n 1 {input}) <(tail -n +2 {input} | shuf -n 100000) > {output}
        """

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

rule correlate_expression:
    input:
        merged=DATA_DIR + "MLE/merged_data.txt"
    output:
        DATA_DIR + "correlation.expression_regressed.{cov_regress}.txt"
    params:
        cov_regress=lambda wildcards: wildcards.cov_regress.split('.')
    shell:
        """
        conda activate pystats
        python scripts/correlate_expression_genotype_PC.py \
            --merged {input.merged} \
            --covariates_to_regress {params.cov_regress} \
            --out {output}
        conda deactivate
        """
