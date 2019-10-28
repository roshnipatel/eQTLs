# Variables and filepaths stored here
include: "scripts/snakemake_variables.py"

shell.executable("/usr/bin/bash")
shell.prefix("source ~/.bashrc; ")

rule all:
    input:
        DATA_DIR + "gencode.v30.genes.gtf"

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
        vcf=DATA_DIR + "mesa.1331samples.genotypes.pass.phased.maf01.vcf.gz"
    output:
        chr_list=DATA_DIR + "chr_list.txt"
    shell:
        """
        conda activate tabix-env
        tabix --list-chroms {input.vcf} > {output.chr_list}
        conda deactivate
        """

rule select_samples:
    input:
        DATA_DIR + "MESA_TOPMed_RNASeqSamples_11022018.txt"
    output:
        DATA_DIR + "sample_participant_lookup.txt"
    shell:
        """
        conda activate py36
        python scripts/select_samples.py {input} {output}
        conda deactivate
        """

rule normalize_expression:
    input:
        tpm=DATA_DIR + "",
        counts=DATA_DIR + "",
        anno=rules.annotate_genes.output,
        sample_map=DATA_DIR + "sample_participant_lookup.txt",
        chr_list=rules.make_chr_list.output
    output:
        bed=DATA_DIR + "MESA.expression.bed.gz",
        idx=DATA_DIR + "MESA.expression.bed.gz.tbi"
    singularity:
        "docker://broadinstitute/gtex_eqtl:V8"
    params:
        prefix="MESA"
    shell:
        """
        scripts/src/eqtl_prepare_expression.py {input.tpm} {input.counts} {input.anno} \
            ${input.sample_map} ${input.chr_list} {params.prefix} \
            --tpm_threshold 0.1 \
            --count_threshold 6 \
            --sample_frac_threshold 0.2 \
            --normalization_method tmm
        """

rule calculate_PEER:
    input:
        bed=rules.normalize_expression.output.bed
    output:
        res=DATA_DIR + "MESA.PEER_residuals.txt",
        alpha=DATA_DIR + "MESA.PEER_alpha.txt",
        cov=DATA_DIR + "MESA.PEER_covariates.txt"
    params:
        prefix="MESA",
        num_PEER=30
    singularity:
        "docker://broadinstitute/gtex_eqtl:V8"
    shell:
        """
        Rscript scripts/src/run_PEER.R {input.bed} {params.prefix} {params.num_PEER}
        """

rule prep_covariates:
    input:
        genotypes=DATA_DIR + "",
        sample_data=DATA_DIR + ""
    output:
        PCs = DATA_DIR + "MESA.PC.txt",
        addl_cov = DATA_DIR + "MESA.sample_covariates.txt"
    shell:
        """

        """

rule combine_covariates:
    input:
        PEER=rules.calculate_PEER.output.cov,
        PCs=rules.prep_covariates.output.PCs,
        addl_cov=rules.prep_covariates.output.addl_cov
    output:
        DATA_DIR + "MESA.combined_covariates.txt"
    params:
        prefix="MESA"
    singularity:
        "docker://broadinstitute/gtex_eqtl:V8"
    shell:
        """
        scripts/src/combine_covariates.py {input.PEER} {params.prefix} \
            --genotype_pcs {input.PCs} --add_covariates {input.addl_cov}
        """
