# Variables and filepaths stored here
include: "scripts/snakemake_variables.py"

shell.executable("/usr/bin/bash")
shell.prefix("source ~/.bashrc; ")

rule all:
    input:
        DATA_DIR + "MESA.combined_covariates.txt"

rule annotate_genes:
    input:
        DATA_DIR + "gencode.v30.annotation.gtf"
    output:
        DATA_DIR + "gencode.v30.genes.gtf"
    shell:
        """
        conda activate gene-annotation-env
        python scripts/collapse_annotation.py {input} {output}
        conda deactivate
        """

rule make_chr_list:
    input:
        vcf=ancient(DATA_DIR + "genotype_freeze.6a.pass_only.phased.mesa_1331samples.maf01.biallelic.vcf.gz")
    output:
        chr_list=DATA_DIR + "chr_list.txt"
    shell:
        """
        conda activate tabix-env
        tabix -l {input.vcf} > {output.chr_list}
        conda deactivate
        """

rule select_samples:
    input:
        sample_data=DATA_DIR + "MESA_TOPMed_RNASeqSamples_11022018.txt"
    output:
        DATA_DIR + "sample_participant_lookup.txt"
    shell:
        """
        conda activate py36
        python scripts/select_samples.py {input.sample_data} {output}
        conda deactivate
        """

rule normalize_expression:
    input:
        tpm=DATA_DIR + "TOPMed_MESA_RNAseq_Pilot_expression_data/TOPMed_MESA_RNAseq_Pilot_RNASeQCv1.1.9.gene_tpm.gct.gz",
        counts=DATA_DIR + "TOPMed_MESA_RNAseq_Pilot_expression_data/TOPMed_MESA_RNAseq_Pilot_RNASeQCv1.1.9.gene_reads.gct.gz",
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
        scripts/src/eqtl_prepare_expression.py {input.tpm} {input.counts} {input.anno} \
            {input.sample_map} {input.chr_list} {params.prefix} \
            --tpm_threshold 0.1 \
            --count_threshold 6 \
            --sample_frac_threshold 0.2 \
            --normalization_method tmm
        conda deactivate
        """

rule calculate_PEER:
    input:
        bed=rules.normalize_expression.output.bed
    output:
        res=DATA_DIR + "MESA.PEER_residuals.txt",
        alpha=DATA_DIR + "MESA.PEER_alpha.txt",
        cov=DATA_DIR + "MESA.PEER_covariates.txt"
    params:
        prefix=DATA_DIR + "MESA",
        num_PEER=30
    shell:
        """
        Rscript scripts/src/run_PEER.R {input.bed} {params.prefix} {params.num_PEER}
        """

rule MAF_filter:
    input:
        ancient(DATA_DIR + "genotype_freeze.6a.pass_only.phased.mesa_1331samples.maf01.biallelic.vcf.gz")
    output:
        vcf=DATA_DIR + "mesa.maf05.vcf.gz"
    shell:
        """
        conda activate bcftools-env
        bcftools view -i 'MAF[0] > 0.05' -Ob {input} | \
            bcftools view --genotype ^miss --phased -Oz -o {output.vcf}
        conda deactivate
        """

rule find_indep_snps:
    input:
        vcf=rules.MAF_filter.output.vcf
    output:
        keep=DATA_DIR + "mesa.prune.in",
        exclude=temp(expand(DATA_DIR + "mesa.{ext}", ext=["log", "nosex", "prune.out"]))
    params:
        out_file=DATA_DIR + "mesa"
    shell:
        """
        conda activate plink-env
        plink --vcf {input} --vcf-half-call m --allow-no-sex \
            --keep-allele-order --indep-pairwise 50 5 0.5 --out {params.out_file}
        conda deactivate
        """

rule prune:
    input:
        vcf=rules.MAF_filter.output.vcf,
        snps=rules.find_indep_snps.output.keep
    output:
        DATA_DIR + "mesa.maf05.prune05.bed"
    shell:
        """
        conda activate plink-env
        plink --vcf {input.vcf} --make-bed --extract {input.snps} --out {output}
        conda deactivate
        """

rule genotype_PC:
    input:
        rules.prune.output
    output:
        DATA_DIR + "mesa.maf05.prune05.eigenvec"
    params:
        prefix=DATA_DIR + "mesa.maf05.prune05"
    shell:
        """
        conda activate plink-env
        plink --bfile {params.prefix} --pca --out {output}
        conda deactivate
        """

rule prep_covariates:
    input:
        samples=rules.select_samples.output,
        indiv_data=DATA_DIR + "MESA_sample_info.csv",
        sample_data=DATA_DIR + "MESA_TOPMed_RNASeqSamples_11022018.txt"
    output:
        DATA_DIR + "mesa.sample_covariates.txt"
    shell:
        """
        conda activate py36
        python scripts/prep_covariates.py --samples {input.samples} \
            --indiv_metadata {input.indiv_data} \
            --sample_metadata {input.sample_data} \
            --out {output}
        conda deactivate
        """

rule combine_covariates:
    input:
        PEER=rules.calculate_PEER.output.cov,
        PC=rules.genotype_PC.output,
        addl_cov=rules.prep_covariates.output
    output:
        DATA_DIR + "MESA.combined_covariates.txt"
    params:
        prefix=DATA_DIR + "MESA",
        PC_trim=DATA_DIR + "mesa.pc.txt"
    shell:
        """
        conda activate py36
        cut -f 2-22 -d$' ' {input.PC} > {params.PC_trim}
        scripts/src/combine_covariates.py {input.PEER} {params.prefix} \
            {params.PC_trim} {input.addl_cov}
        conda deactivate
        """
