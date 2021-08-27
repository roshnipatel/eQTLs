# Variables and filepaths stored here
include: "scripts/snakemake_variables.py"

shell.executable("/usr/bin/bash")
shell.prefix("source ~/.bashrc; ")

# Define population for ascertainment of significant SNPs
ASC_POP = "Eur"

# Define lists for expanding input for rule all
NO_BETA_COV_LIST = ["seq_center.exam", "race_Afr.race_Eur.seq_center.exam", 
                    "global_ancestry.genotype_PC1.race_Afr.race_Eur.seq_center.exam",
                    "global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1"]
SINGLE_BETA_COV_LIST = ["global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1"]
NO_INT_COV_LIST = ["global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1"]
COV_LIST = ["global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1"]
GROUPS = ["test", "control"]
IND_SUBSETS = ["AA", "EA"]

rule all:
    input:
        ### For eQTL calling
        # DATA_DIR + "QTL_output/hits_ascertainment_" + ASC_POP + ".txt",
        # DATA_DIR + "QTL_output/hits_reestimation_primary_EA.txt",
        # DATA_DIR + "QTL_output/hits_reestimation_primary_AA.txt",
        # DATA_DIR + "QTL_output/hits_reestimation_primary_Eur.txt",
        # DATA_DIR + "QTL_output/hits_reestimation_primary_Afr.txt",
        # 
        ### For likelihood model stuff:
        # PROTECTED_DATA_DIR + "merged_data.txt",
        # expand(DATA_DIR + "model_fitting/{cov}/{group}.mode_fit_delta.optimize_{param}.txt", 
        #        param=["delta", "betas"], cov=COV_LIST, group=GROUPS),
        # expand(DATA_DIR + "model_fitting/{cov}/test.mode_fit_delta.idv_bin_{bin}.optimize_{param}.txt", 
        #        param=["delta", "betas"], cov=COV_LIST, bin=[1, 2, 3, 4, "low", "high"]),
        expand(DATA_DIR + "model_fitting/{cov}/test.jackknife_summary.txt", 
               cov=COV_LIST),
        # expand(DATA_DIR + "model_fitting/{cov}/test.idv_bin_{bin}.bootstrap_summary.txt", 
        #        cov=COV_LIST, bin=[1, 2, 3, 4, "low", "high"]),
        # expand(DATA_DIR + "model_fitting/{cov}/{group}.likelihood_vs_delta.txt", 
        #        cov=COV_LIST, group=GROUPS),
        # expand(DATA_DIR + "model_fitting/{cov}/{group}.mode_fit_delta.ind_{ind}.optimize.anova_analysis.txt", 
        #        cov=COV_LIST, group=["test"], ind=IND_SUBSETS),
        # expand(DATA_DIR + "model_fitting/{cov}/{group}.mode_fit_no_beta.ind_all.optimize.anova_analysis.txt", 
        #        cov=NO_BETA_COV_LIST, group=["test"]),
        # expand(DATA_DIR + "model_fitting/{cov}/{group}.mode_fit_single_beta.ind_all.optimize.anova_analysis.txt", 
        #        cov=SINGLE_BETA_COV_LIST, group=["test"]),
        # expand(DATA_DIR + "model_fitting/{cov}/{group}.mode_no_interaction.ind_{ind}.optimize.anova_analysis.txt", 
        #        cov=NO_INT_COV_LIST, group=["test"], ind=IND_SUBSETS),
        # expand(DATA_DIR + "model_fitting/simulated_data/vA{var_afr}.vE{var_eur}.errA{afr_error}.errE{eur_error}/betas.txt", 
        #        var_afr=[0.08], var_eur=[0.15], eur_error=[0.12], afr_error=[0.17]),
        # expand(DATA_DIR + "model_fitting/simulated_optimization/vA{var_afr}.vE{var_eur}.errA{afr_error}.errE{eur_error}/combined_optimization_outputs.txt", 
        #        eur_error=[0.12], afr_error=[0.17], 
        #        var_afr=[0.08], var_eur=[0.15])
        expand(DATA_DIR + "model_fitting/{cov}/{group}.mode_fit_delta.ind_{ind}.optimize.PGS_PVE_analysis.txt",
               cov=COV_LIST, group="test", ind=["EA", "AA"])
        ### For covariate/PC stuff:
        # DATA_DIR + "correlation.expression.txt"

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
        sample_data=PROTECTED_DATA_DIR + RNASEQ_METADATA,
        exclusion_list=PROTECTED_DATA_DIR + EXCLUSION_FILE
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
            --out {params.prefix}
        conda deactivate
        """

rule make_chr_list:
    input:
        vcf=PROTECTED_DATA_DIR + VCF
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
        tpm=PROTECTED_DATA_DIR + TPM,
        counts=PROTECTED_DATA_DIR + READS,
        anno=rules.annotate_genes.output,
        sample_map=DATA_DIR + "select_samples_{anc}.txt",
        chr_list=rules.make_chr_list.output,
        flag=DATA_DIR + "vcf_filter_{anc}.done"
    output:
        bed=PROTECTED_DATA_DIR + "mesa.{anc}.expression.bed.gz",
        idx=PROTECTED_DATA_DIR + "mesa.{anc}.expression.bed.gz.tbi"
    params:
        prefix=PROTECTED_DATA_DIR + "mesa.{anc}"
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
        python {GENE_LIST_SCRIPT} --bed {input} --out {output}
        conda deactivate
        """

############################# GENERATE COVARIATES ##############################

rule snp_filter:
    input:
        vcf=PROTECTED_DATA_DIR + VCF
    output:
        vcf=PROTECTED_DATA_DIR + "mesa.filt.bcf.gz",
        idx=PROTECTED_DATA_DIR + "mesa.filt.bcf.gz.csi"
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
        vcf=PROTECTED_DATA_DIR + "mesa.{anc,[A-Za-z]+}.bcf.gz",
        idx=PROTECTED_DATA_DIR + "mesa.{anc,[A-Za-z]+}.bcf.gz.csi",
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
        vcf=PROTECTED_DATA_DIR + "mesa.{anc}.maf05.vcf.gz",
        idx=PROTECTED_DATA_DIR + "mesa.{anc}.maf05.vcf.gz.csi"
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
        keep=temp(PROTECTED_DATA_DIR + "mesa.{anc}.maf05.prune.in"),
        exclude=temp(expand(PROTECTED_DATA_DIR + "mesa.{{anc}}.maf05.{ext}", 
                            ext=["log", "nosex", "prune.out"]))
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
        PROTECTED_DATA_DIR + "mesa.{anc}.prune.maf05.vcf.gz"
    shell:
        """
        conda activate bcftools-env
        bcftools view -i 'ID=@{input.snps}' -Oz -o {output} {input.vcf}
        conda deactivate
        """

rule PCA:
    input:
        expand(PROTECTED_DATA_DIR + "mesa.{anc}.{{filetype}}.gz", anc=["Afr", "Eur"])
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
        ganc=PROTECTED_DATA_DIR + "combined_global_anc_frac.txt",
        pc=expand(rules.PCA.output.pc, filetype=["prune.maf05.vcf", "expression.bed"]),
        sample_data=PROTECTED_DATA_DIR + RNASEQ_METADATA,
        flag=expand(DATA_DIR + "vcf_filter_{anc}.done", anc=["Afr", "Eur"])
    output:
        PROTECTED_DATA_DIR + "mesa.sample_covariates.txt"
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

rule extract_linreg_covariates_afr:
    input:
        PROTECTED_DATA_DIR + "mesa.sample_covariates.txt"
    output:
        PROTECTED_DATA_DIR + "mesa.Afr.sample_covariates_for_regression.txt"
    shell:
        """
        cut -f 2,3,4,6,8 {input} > {output}
        """

rule extract_linreg_covariates_eur:
    input:
        PROTECTED_DATA_DIR + "mesa.sample_covariates.txt"
    output:
        PROTECTED_DATA_DIR + "mesa.Eur.sample_covariates_for_regression.txt"
    shell:
        """
        cut -f 2,3,4,6,10 {input} > {output}
        """

############################## PARTITION SAMPLES ###############################

rule make_anc_bed:
    input:
        expand(PROTECTED_DATA_DIR + BED_DIR + "chr{chr}/job_complete.txt", chr=CHROMS)
    output:
        PROTECTED_DATA_DIR + "anc_tracts.bed"
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
        chrom_lengths=DATA_DIR + CHR_LEN 
    output:
        DATA_DIR + "genes.bed"
    shell:
        """
        conda activate py36
        python {GENE_BED_SCRIPT} --genes {input.genes} \
            --window {WINDOW} \
            --chrom_lengths {input.chrom_lengths} \
            --out {output}
        conda deactivate
        """

rule intersect_tracts_genes:
    input:
        tracts=rules.make_anc_bed.output,
        genes=rules.make_genes_bed.output
    output:
        PROTECTED_DATA_DIR + "intersection_anc_genes.bed"
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
        directory(DATA_DIR + "QTL_sample_input"),
        DATA_DIR + "ascertainment.txt",
        DATA_DIR + "validation.txt"
    params:
        output_dir=DATA_DIR + "QTL_sample_input"
    shell:
        """
        mkdir -p {params.output_dir}/ascertainment/Eur
        mkdir -p {params.output_dir}/ascertainment/Afr
        mkdir -p {params.output_dir}/reestimation_primary/Afr
        mkdir -p {params.output_dir}/reestimation_primary/Eur
        mkdir -p {params.output_dir}/reestimation_primary/AA
        mkdir -p {params.output_dir}/reestimation_primary/EA
        conda activate py36
        python {PARTITION_SCRIPT} \
            --intersect {input.intersect} \
            --afr_samples {input.afr_samples} \
            --eur_samples {input.eur_samples} \
            --genes {input.genes} \
            --ascertainment {output[1]} \
            --validation {output[2]} \
            --ascertainment_pop {ASC_POP} \
            --out_dir {params.output_dir}
        conda deactivate
        """

############################ LIN REG EQTL CALLING ##############################

rule subset_geno:
    input:
        vcf=PROTECTED_DATA_DIR + "mesa.{anc}.maf05.vcf.gz",
        idx=PROTECTED_DATA_DIR + "mesa.{anc}.maf05.vcf.gz.csi",
        gene_window=rules.make_genes_bed.output
    output:
        PROTECTED_DATA_DIR + "QTL_geno_input/{anc}/{gene}.vcf.gz"
    params:
        out_dir=PROTECTED_DATA_DIR + "QTL_geno_input/{anc}"
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
        bed=PROTECTED_DATA_DIR + "mesa.{anc}.expression.bed.gz"
    output:
        PROTECTED_DATA_DIR + "QTL_pheno_input/{anc}/{gene}.txt"
    params:
        out_dir=PROTECTED_DATA_DIR + "QTL_pheno_input/{anc}"
    shell:
        """
        mkdir -p {params.out_dir}
        cat <(zcat {input.bed} | head -n 1) <(zgrep -m 1 {wildcards.gene} {input.bed}) > {output}
        """

# We have a separate rule for ascertaining SNPs vs estimating effect sizes in case
# we want to use permutations for ascertainment.
rule ascertain_effect_sizes:
    input:
        geno_input=lambda wildcards: expand(rules.subset_geno.output, anc="Eur", \
            gene=wildcards.gene) if wildcards.anc[0] == "E" else \
            expand(rules.subset_geno.output, anc="Afr", gene=wildcards.gene),
        pheno_input=lambda wildcards: expand(rules.subset_pheno.output, anc="Eur", \
            gene=wildcards.gene) if wildcards.anc[0] == "E" else \
            expand(rules.subset_pheno.output, anc="Afr", gene=wildcards.gene),
        sample_input=DATA_DIR + "QTL_sample_input/ascertainment/{anc}/{gene}.txt",
        covariates=lambda wildcards: PROTECTED_DATA_DIR + \
            "mesa.Eur.sample_covariates_for_regression.txt" if \
            wildcards.anc[0] == "E" else PROTECTED_DATA_DIR + \
            "mesa.Afr.sample_covariates_for_regression.txt"
    output:
        DATA_DIR + "QTL_output/ascertainment/{anc}/{gene}.txt"
    params:
        output_dir=DATA_DIR + "QTL_output/ascertainment/{anc}"
    shell:
        """
        mkdir -p {params.output_dir}
        conda activate pystats
        python {ESTIMATION_SCRIPT} \
            --genotypes {input.geno_input} \
            --phenotypes {input.pheno_input} \
            --samples {input.sample_input} \
            --covariates {input.covariates} \
            --out {output} \
        conda deactivate
        """

rule estimate_effect_sizes:
    input:
        geno_input=lambda wildcards: expand(rules.subset_geno.output, anc="Eur", \
            gene=wildcards.gene) if wildcards.anc[0] == "E" else \
            expand(rules.subset_geno.output, anc="Afr", gene=wildcards.gene),
        pheno_input=lambda wildcards: expand(rules.subset_pheno.output, anc="Eur", \
            gene=wildcards.gene) if wildcards.anc[0] == "E" else \
            expand(rules.subset_pheno.output, anc="Afr", gene=wildcards.gene),
        sample_input=DATA_DIR + "QTL_sample_input/reestimation_{type}/{anc}/{gene}.txt",
        covariates=lambda wildcards: PROTECTED_DATA_DIR + \
            "mesa.Eur.sample_covariates_for_regression.txt" if \
            wildcards.anc[0] == "E" else PROTECTED_DATA_DIR + \
            "mesa.Afr.sample_covariates_for_regression.txt"
    output:
        DATA_DIR + "QTL_output/{type}/reestimation_{anc}/{gene}.txt"
    params:
        output_dir=DATA_DIR + "QTL_output/reestimation_{type}/{anc}"
    shell:
        """
        mkdir -p {params.output_dir}
        conda activate pystats
        python {ESTIMATION_SCRIPT} \
            --genotypes {input.geno_input} \
            --phenotypes {input.pheno_input} \
            --samples {input.sample_input} \
            --covariates {input.covariates} \
            --out {output}
        conda deactivate
        """

def concat_asc_input(wildcards):
    checkpoint_output = checkpoints.partition_samples.get(**wildcards).output[0]
    return expand(DATA_DIR + "QTL_output/ascertainment/{anc}/{gene}.txt",
                  anc=wildcards.anc,
                  gene=glob_wildcards(os.path.join(checkpoint_output, "ascertainment", 
                                      wildcards.anc, "{gene}.txt")).gene)

def concat_est_input(wildcards):
    checkpoint_output = checkpoints.identify_hits.get(**{"anc": ASC_POP}).output[0]
    sig_genes = list(pd.read_csv(checkpoint_output, sep='\t')["gene"])
    return expand(DATA_DIR + "QTL_output/reestimation_{type}/{anc}/{gene}.txt",
                  type=wildcards.type, anc=wildcards.anc,
                  gene=sig_genes)

rule concat_ascertainment:
    input:
        concat_asc_input
    output:
        DATA_DIR + "QTL_output/merged_ascertainment_{anc}.txt"
    params:
        dir=DATA_DIR + "QTL_output/ascertainment/{anc}/"
    shell:
        """
        conda activate py36
        python {CONCAT_QTL_SCRIPT} --dir {params.dir} --out {output}
        conda deactivate 
        """

rule concat_estimation:
    input:
        concat_est_input
    output:
        DATA_DIR + "QTL_output/merged_reestimation_{type}_{anc}.txt"
    params:
        dir=DATA_DIR + "QTL_output/reestimation_{type}/{anc}/"
    shell:
        """
        conda activate py36
        python {CONCAT_QTL_SCRIPT} --dir {params.dir} --out {output}
        conda deactivate 
        """

checkpoint identify_hits:
    input:
        asc=DATA_DIR + "QTL_output/merged_ascertainment_{anc}.txt" 
    output:
        DATA_DIR + "QTL_output/hits_ascertainment_{anc}.txt"
    params:
        DATA_DIR + "QTL_output/"
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

checkpoint extract_hits:
    input:
        merged=DATA_DIR + "QTL_output/merged_reestimation_{type}_{anc}.txt",
        hits=DATA_DIR + "QTL_output/hits_ascertainment_" + ASC_POP + ".txt"
    output:
        DATA_DIR + "QTL_output/hits_reestimation_{type}_{anc}.txt"
    params:
        DATA_DIR + "QTL_output/"
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

############################# LIKELIHOOD MODEL #################################

def fetch_genes(f):
    genes, gene_idx = [], None
    for line in f.readlines():
        if not gene_idx:
            gene_idx = line.strip().split('\t').index("gene")
            continue
        curr_gene = line.strip().split('\t')[gene_idx]
        genes.append(curr_gene)
    return(genes)

def expand_geno_pheno_data(wildcards):
    with open(checkpoints.extract_hits.get(type="primary", anc="Afr").output[0], 'r') as f:
        afr_genes = fetch_genes(f)
    with open(checkpoints.extract_hits.get(type="primary", anc="Eur").output[0], 'r') as f:
        eur_genes = fetch_genes(f)
    common_genes = list(set(afr_genes) & set(eur_genes))
    return(expand(PROTECTED_DATA_DIR + "QTL_geno_input/{anc}/{gene}.vcf.gz", 
                  anc=["Afr", "Eur"], gene=common_genes) + 
           expand(PROTECTED_DATA_DIR + "QTL_pheno_input/{anc}/{gene}.txt", 
                  anc=["Afr", "Eur"], gene=common_genes))

rule merge_data:
    input:
        tracts=rules.make_anc_bed.output,
        cov=rules.prep_covariates.output,
        afr_hits=ancient(DATA_DIR + "QTL_output/hits_reestimation_primary_Afr.txt"),
        eur_hits=ancient(DATA_DIR + "QTL_output/hits_reestimation_primary_Eur.txt"),
        data=expand_geno_pheno_data
    output:
        merged=PROTECTED_DATA_DIR + "merged_data.txt"
    params:
        out_dir=PROTECTED_DATA_DIR
    shell:
        """
        mkdir -p {params.out_dir}
        conda activate pystats
        python {MERGE_SCRIPT} --tracts {input.tracts} \
            --afr_hits {input.afr_hits} \
            --eur_hits {input.eur_hits} \
            --covariates {input.cov} \
            --out {output.merged}
        conda deactivate
        """

rule optimize:
    input:
        merged=PROTECTED_DATA_DIR + "merged_data.txt",
        ascertainment=DATA_DIR + "ascertainment.txt"
    output:
        delta=DATA_DIR + "model_fitting/{cov}/{group}.mode_{mode, ^(.)}.optimize_delta.txt",
        betas=DATA_DIR + "model_fitting/{cov}/{group}.mode_{mode, ^(.)}.optimize_betas.txt"
    params:
        cov=lambda wildcards: wildcards.cov.split('.'),
        dir=DATA_DIR + "model_fitting/{cov}"
    shell:
        """
        mkdir -p {params.dir}
        conda activate pystats
        python {OPTIMIZE_SCRIPT} --merged {input.merged} \
            --ascertainment {input.ascertainment} \
            --max_iter {MAX_ITER} \
            --covariates {params.cov} \
            --group {wildcards.group} \
            --mode {wildcards.mode} \
            --betas_out {output.betas} \
            --delta_out {output.delta} 
        conda deactivate
        """

rule optimize_idv_bin:
    input:
        merged=PROTECTED_DATA_DIR + "merged_data.txt",
        bin=DATA_DIR + "bins/idv_bin_{bin}.txt",
        ascertainment=DATA_DIR + "ascertainment.txt"
    output:
        delta=DATA_DIR + "model_fitting/{cov}/{group}.mode_{mode}.idv_bin_{bin}.optimize_delta.txt",
        betas=DATA_DIR + "model_fitting/{cov}/{group}.mode_{mode}.idv_bin_{bin}.optimize_betas.txt"
    params:
        cov=lambda wildcards: wildcards.cov.split('.'),
        dir=DATA_DIR + "model_fitting/{cov}"
    shell:
        """
        mkdir -p {params.dir}
        conda activate pystats
        python {OPTIMIZE_SCRIPT} --merged {input.merged} \
            --ascertainment {input.ascertainment} \
            --max_iter {MAX_ITER} \
            --covariates {params.cov} \
            --group {wildcards.group} \
            --mode {wildcards.mode} \
            --idv_bin {input.bin} \
            --betas_out {output.betas} \
            --delta_out {output.delta} 
        conda deactivate
        """

rule optimize_remove_val:
    input:
        merged=PROTECTED_DATA_DIR + "merged_data.txt",
        validation=DATA_DIR + "validation.txt",
        ascertainment=DATA_DIR + "ascertainment.txt"
    output:
        delta=DATA_DIR + "model_fitting/{cov}/{group}.mode_{mode}.optimize_wo_val_delta.txt",
        betas=DATA_DIR + "model_fitting/{cov}/{group}.mode_{mode}.optimize_wo_val_betas.txt",
    params:
        cov=lambda wildcards: wildcards.cov.split('.'),
        dir=DATA_DIR + "model_fitting/{cov}"
    shell:
        """
        mkdir -p {params.dir}
        conda activate pystats
        python {OPTIMIZE_SCRIPT} --merged {input.merged} \
            --validation {input.validation} \
            --ascertainment {input.ascertainment} \
            --max_iter {MAX_ITER} \
            --covariates {params.cov} \
            --group {wildcards.group} \
            --mode {wildcards.mode} \
            --betas_out {output.betas} \
            --delta_out {output.delta} 
        conda deactivate
        """

rule simulate_betas:
    input:
        afr=ancient("data/QTL_output/hits_reestimation_primary_Afr.txt"),
        eur=ancient("data/QTL_output/hits_reestimation_primary_Eur.txt")
    output:
        DATA_DIR + "model_fitting/simulated_data/vA{var_afr}.vE{var_eur}.errA{afr_error}.errE{eur_error}/betas.txt"
    params:
        dir=DATA_DIR + "model_fitting/simulated_data/vA{var_afr}.vE{var_eur}.errA{afr_error}.errE{eur_error}"
    shell:
        """
        mkdir -p {params.dir}
        conda activate pystats
        python scripts/simulate_betas.py \
            --var_afr {wildcards.var_afr} \
            --var_eur {wildcards.var_eur} \
            --error_afr {wildcards.afr_error} \
            --error_eur {wildcards.eur_error} \
            --afr_hits {input.afr} \
            --eur_hits {input.eur} \
            --out {output}
        conda deactivate
        """

rule simulate_data:
    input:
        rules.simulate_betas.output
    output:
        DATA_DIR + "model_fitting/simulated_data/vA{var_afr}.vE{var_eur}.errA{afr_error}.errE{eur_error}/delta{delta}.data.{idx}.txt"
    params:
        dir=DATA_DIR + "model_fitting/simulated_data/vA{var_afr}.vE{var_eur}.errA{afr_error}.errE{eur_error}"
    shell:
        """
        conda activate pystats
        python scripts/simulate_data.py \
            --sim_betas {input} \
            --error_afr {wildcards.afr_error} \
            --error_eur {wildcards.eur_error} \
            --delta {wildcards.delta} \
            --out {output}
        conda deactivate
        """

rule optimize_simulated_data:
    input:
        rules.simulate_data.output
    output:
        delta=temp(DATA_DIR + "model_fitting/simulated_optimization/vA{var_afr}.vE{var_eur}.errA{afr_error}.errE{eur_error}/delta{delta}.optimized_delta.{idx}.txt"),
        betas=temp(DATA_DIR + "model_fitting/simulated_optimization/vA{var_afr}.vE{var_eur}.errA{afr_error}.errE{eur_error}/delta{delta}.optimized_betas.{idx}.txt")
    params:
        dir=DATA_DIR + "model_fitting/simulated_optimization/vA{var_afr}.vE{var_eur}.errA{afr_error}.errE{eur_error}"
    shell:
        """
        mkdir -p {params.dir}
        conda activate pystats
        python {OPTIMIZE_SCRIPT} --merged {input} \
            --max_iter {MAX_ITER} \
            --covariates \
            --mode fit_delta \
            --betas_out {output.betas} \
            --delta_out {output.delta} 
        conda deactivate
        """

rule parse_simulation_results:
    input:
        expand(DATA_DIR + "model_fitting/simulated_optimization/vA{{var_afr}}.vE{{var_eur}}.errA{{afr_error}}.errE{{eur_error}}/delta{delta}.optimized_delta.{idx}.txt", 
               idx=[str(i).zfill(2) for i in range(10)], delta=[i/10 for i in range(11)])
    output:
        DATA_DIR + "model_fitting/simulated_optimization/vA{var_afr}.vE{var_eur}.errA{afr_error}.errE{eur_error}/combined_optimization_outputs.txt"
    params:
        dir=DATA_DIR + "model_fitting/simulated_optimization/vA{var_afr}.vE{var_eur}.errA{afr_error}.errE{eur_error}/"
    shell:
        """
        conda activate pystats
        python {SIM_PARSER} --delta_dir {params.dir} \
            --out {output}
        conda deactivate
        """

rule bootstrap:
    input:
        merged=rules.merge_data.output.merged,
        bin=ancient(DATA_DIR + "bins/idv_bin_{bin}.txt"),
        ascertainment=DATA_DIR + "ascertainment.txt"
    output:
        delta=temp(DATA_DIR + "model_fitting/{cov}/{group}.idv_bin_{bin}.bootstrap.delta_{idx}.txt"),
        betas=temp(DATA_DIR + "model_fitting/{cov}/{group}.idv_bin_{bin}.bootstrap.betas_{idx}.txt")
    params:
        cov=lambda wildcards: wildcards.cov.split('.'),
        dir=DATA_DIR + "model_fitting/{cov}"
    shell:
        """
        mkdir -p {params.dir}
        conda activate pystats
        python {OPTIMIZE_SCRIPT} --merged {input.merged} \
            --ascertainment {input.ascertainment} \
            --max_iter {MAX_ITER} \
            --bootstrap \
            --mode fit_delta \
            --idv_bin {input.bin} \
            --group {wildcards.group} \
            --covariates {params.cov} \
            --betas_out {output.betas} \
            --delta_out {output.delta} 
        conda deactivate
        """

rule parse_bootstrap:
    input:
        expand(DATA_DIR + "model_fitting/{{cov}}/{{group}}.bootstrap.delta_{idx}.txt", 
               idx=[str(i).zfill(3) for i in range(1000)])
    output:
        DATA_DIR + "model_fitting/{cov}/{group}.bootstrap_summary.txt",
        DATA_DIR + "model_fitting/{cov}/{group}.bootstrap_all.txt"
    params:
        prefix=DATA_DIR + "model_fitting/{cov}/{group}.bootstrap_"
    shell:
        """
        conda activate py36
        python {BOOTSTRAP_SCRIPT} --bootstrap_files {input} --out {params.prefix}
        conda deactivate
        """

rule jackknife:
    input:
        merged=rules.merge_data.output.merged,
        ascertainment=DATA_DIR + "ascertainment.txt"
    output:
        delta=temp(DATA_DIR + "model_fitting/{cov}/{group}.jackknife.delta_{gene}.txt"),
        betas=temp(DATA_DIR + "model_fitting/{cov}/{group}.jackknife.betas_{gene}.txt")
    params:
        cov=lambda wildcards: wildcards.cov.split('.'),
        dir=DATA_DIR + "model_fitting/{cov}"
    shell:
        """
        mkdir -p {params.dir}
        conda activate pystats
        python {OPTIMIZE_SCRIPT} --merged {input.merged} \
            --ascertainment {input.ascertainment} \
            --max_iter {MAX_ITER} \
            --jackknife {wildcards.gene} \
            --mode fit_delta \
            --group {wildcards.group} \
            --covariates {params.cov} \
            --betas_out {output.betas} \
            --delta_out {output.delta} 
        conda deactivate
        """

def expand_jackknife(wildcards):
    with open(checkpoints.extract_hits.get(type="primary", anc="Afr").output[0], 'r') as f:
        afr_genes = fetch_genes(f)
    with open(checkpoints.extract_hits.get(type="primary", anc="Eur").output[0], 'r') as f:
        eur_genes = fetch_genes(f)
    common_genes = list(set(afr_genes) & set(eur_genes))
    return(expand(DATA_DIR + "model_fitting/{cov}/{group}.jackknife.delta_{gene}.txt",
                  cov=wildcards.cov, group=wildcards.group, gene=common_genes)) 

rule parse_jackknife:
    input:
        expand_jackknife
    output:
        DATA_DIR + "model_fitting/{cov}/{group}.jackknife_summary.txt",
        DATA_DIR + "model_fitting/{cov}/{group}.jackknife_all.txt"
    params:
        prefix=DATA_DIR + "model_fitting/{cov}/{group}.jackknife_"
    shell:
        """
        conda activate py36
        python {BOOTSTRAP_SCRIPT} --bootstrap_files {input} --out {params.prefix}
        conda deactivate
        """

rule likelihood:
    input:
        merged=rules.merge_data.output.merged,
        ascertainment=DATA_DIR + "ascertainment.txt"
    output:
        DATA_DIR + "model_fitting/{cov}/{group}.likelihood_vs_delta.txt"
    params:
        cov=lambda wildcards: wildcards.cov.split('.'),
        dir=DATA_DIR + "model_fitting/{cov}"
    shell:
        """
        mkdir -p {params.dir}
        conda activate pystats
        python {LIKELIHOOD_SCRIPT} --merged {input.merged} \
            --ascertainment {input.ascertainment} \
            --group {wildcards.group} \
            --covariates {params.cov} \
            --out {output}
        conda deactivate
        """

############################## MISC. ANALYSES ##################################

rule sample_merged:
    input:
        DATA_DIR + "QTL_output/merged_{type}_{anc}.txt"
    output:
        DATA_DIR + "QTL_output/randsamp_merged_{type}_{anc}.txt"
    shell:
        """
        cat <(head -n 1 {input}) <(tail -n +2 {input} | shuf -n 100000) > {output}
        """

rule correlate_expression:
    input:
        merged=rules.merge_data.output.merged
    output:
        DATA_DIR + "correlation.expression.txt"
    shell:
        """
        conda activate pystats
        python {EXP_CORRELATION_SCRIPT} \
            --merged {input.merged} \
            --out {output}
        conda deactivate
        """

rule std_dev_analysis:
    input:
        rules.select_samples.output,
        PROTECTED_DATA_DIR + READS,
        expand(rules.normalize_expression.output.bed, anc=["Afr", "Eur"])
    output:
        DATA_DIR + "variances.txt"
    shell:
        """
        conda activate pystats
        python {SD_SCRIPT}
        conda deactivate
        """

rule anova_analysis:
    input:
        betas=DATA_DIR + "model_fitting/{cov}/{group}.mode_{mode}.optimize_wo_val_betas.txt",
        delta=DATA_DIR + "model_fitting/{cov}/{group}.mode_{mode}.optimize_wo_val_delta.txt",
        merged=PROTECTED_DATA_DIR + "merged_data.txt",
        validation=DATA_DIR + "validation.txt"
    output:
        DATA_DIR + "model_fitting/{cov}/{group}.mode_{mode}.ind_{ind}.optimize.anova_analysis.txt"
    params:
        terms=lambda wildcards: wildcards.cov.split('.')
    shell:
        """
        conda activate pystats
        python {ANOVA_SCRIPT} \
            --merged {input.merged} \
            --betas {input.betas} \
            --delta {input.delta} \
            --ind {wildcards.ind} \
            --validation {input.validation} \
            --model_terms {params.terms} \
            --residualized_terms \
            --mode {wildcards.mode} \
            --out {output}
        conda deactivate
        """

rule PGS_PVE_analysis:
    input:
        betas=DATA_DIR + "model_fitting/{cov}/{group}.mode_fit_delta.optimize_wo_val_betas.txt",
        delta=DATA_DIR + "model_fitting/{cov}/{group}.mode_fit_delta.optimize_wo_val_delta.txt",
        merged=PROTECTED_DATA_DIR + "merged_data.txt",
        validation=DATA_DIR + "validation.txt"
    output:
        DATA_DIR + "model_fitting/{cov}/{group}.mode_fit_delta.ind_{ind}.optimize.PGS_PVE_analysis.txt"
    params:
        terms=lambda wildcards: wildcards.cov.split('.')
    shell:
        """
        conda activate pystats
        python scripts/PGS_PVE.py \
            --merged {input.merged} \
            --betas {input.betas} \
            --delta {input.delta} \
            --ind {wildcards.ind} \
            --validation {input.validation} \
            --residualized_terms {params.terms} \
            --out {output}
        conda deactivate
        """
