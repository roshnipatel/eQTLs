# Variables and filepaths stored here
include: "scripts/snakemake_variables.py"

shell.executable("/usr/bin/bash")
shell.prefix("source ~/.bashrc; ")

# TODO use glob + expand to properly specify input of make_anc_bed and partition_anc_het

rule all:
    input:
        DATA_DIR + "fastqtl_anc_het_output/merged_anc_het.txt"

########################### GENERATE GENE ANNOTATION ###########################

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

################################ CHOOSE SAMPLES ################################

rule make_geno_list:
    input:
        DATA_DIR + "genotype_freeze.6a.pass_only.phased.mesa_1331samples.maf01.biallelic.vcf.gz"
    output:
        DATA_DIR + "genotyped_individual_IDs.txt"
    shell:
        """
        zgrep -m 1 "^#CHROM" {input} | sed 's/   /\n/g' > {output}
        """

rule select_samples:
    input:
        sample_data=DATA_DIR + "MESA_TOPMed_RNASeqSamples_11022018.txt",
        geno_ID=rules.make_geno_list.output
    output:
        DATA_DIR + "sample_participant_lookup.txt"
    shell:
        """
        conda activate py36
        python scripts/select_samples.py {input.sample_data} {input.geno_ID} {output}
        conda deactivate
        """

############################# NORMALIZE EXPRESSION #############################

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
        PC=rules.genotype_PC.output.trimmed,
        addl_cov=rules.prep_sample_covariates.output
    output:
        DATA_DIR + "mesa.combined_covariates.txt"
    params:
        prefix=DATA_DIR + "mesa",
    shell:
        """
        conda activate py36
        scripts/src/combine_covariates.py {input.PEER} {params.prefix} \
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
        python scripts/combine_tracts.py --tract_dir {params.dir} --out {output}
        conda deactivate
        """

rule make_genes_bed:
    input:
        genes=rules.make_gene_list.output,
        anno=rules.annotate_genes.output,
        chrom_lengths=DATA_DIR + "chrom_lengths.tsv"
    output:
        DATA_DIR + "genes.bed"
    shell:
        """
        conda activate py36
        python scripts/make_genes_bed.py --genes {input.genes} \
            --anno {input.anno} --chrom_lengths {input.chrom_lengths} \
            --out {output}
        conda deactivate
        """

rule intersect_tracts_genes:
    input:
        tracts=DATA_DIR + "anc_tracts.bed",
        genes=DATA_DIR + "genes.bed"
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
        metadata=DATA_DIR + "MESA_sample_info.csv",
        genes=rules.make_gene_list.output
    output:
        directory(DATA_DIR + "fastqtl_sample_input")
    params:
        output_dir=DATA_DIR + "fastqtl_sample_input"
    shell:
        """
        mkdir -p {params.output_dir}/ascertainment/Eur
        mkdir -p {params.output_dir}/estimation/Afr
        mkdir -p {params.output_dir}/estimation/Eur
        conda activate py36
        python scripts/partition_samples.py --intersect {input.intersect} \
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
        vcf=DATA_DIR + "genotype_freeze.6a.pass_only.phased.mesa_1331samples.maf01.biallelic.vcf.gz",
        pheno=ancient(rules.normalize_expression.output.bed),
        sample_input=DATA_DIR + "fastqtl_sample_input/ascertainment/{anc}/{gene}.txt",
        gene_input=rules.prep_pheno_file.output,
        covariates=rules.combine_covariates.output,
        gene_TSS_map=rules.make_gene_list.output
    output:
        DATA_DIR + "fastqtl_output/ascertainment/{anc}/{gene}.txt"
    singularity:
        "gtex_eqtl_V8.sif"
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
           --region chr$CHR:$REGION_START-$REGION_STOP \
           --permute 1000 10000
        """

rule estimate_eQTLs:
    input:
        vcf=DATA_DIR + "genotype_freeze.6a.pass_only.phased.mesa_1331samples.maf01.biallelic.vcf.gz",
        pheno=ancient(rules.normalize_expression.output.bed),
        sample_input=DATA_DIR + "fastqtl_sample_input/estimation/{anc}/{gene}.txt",
        gene_input=rules.prep_pheno_file.output,
        covariates=rules.combine_covariates.output,
        gene_TSS_map=rules.make_gene_list.output
    output:
        DATA_DIR + "fastqtl_output/estimation/{anc}/{gene}.txt"
    singularity:
        "gtex_eqtl_V8.sif"
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
        ls -l {params.dir} | wc -l > {output}
        """

rule identify_hits:
    input:
        DATA_DIR + "fastqtl_output/merged_estimation_Afr.txt",
        DATA_DIR + "fastqtl_output/merged_estimation_Eur.txt",
        DATA_DIR + "fastqtl_output/merged_ascertainment_Eur.txt"
    output:
        DATA_DIR + "hits/hits_estimation_Eur.txt",
        DATA_DIR + "hits/hits_estimation_Afr.txt",
        DATA_DIR + "hits/hits_ascertainment_Eur.txt"
    shell:
        """
        Rscript scripts/identify_hits.R
        """

###################### ANCESTRY-HETEROZYGOUS INDIVIDUALS ######################
# eQTLs in ancestry-heterozygous individuals  have to be called with a separate 
# set of rules because I partition individuals by SNPs rather than by genes.
# (Technically, you could also partition the Afr and Eur ancestry homozygous
# groups by SNPs rather than genes, but I already went to the trouble of calling
# eQTLs in those groups and I don't see any reason to duplicate the work.)

checkpoint partition_anc_het:
    input:
        hits=DATA_DIR + "hits/hits_ascertainment_Eur.txt",
        chrom_lengths=DATA_DIR + "chrom_lengths.tsv"
    output:
        samples=directory(DATA_DIR + "fastqtl_anc_het_sample_input"),
        regions=directory(DATA_DIR + "fastqtl_region_input"),
    params:
        sample_dir=DATA_DIR + "fastqtl_anc_het_sample_input/",
        region_dir=DATA_DIR + "fastqtl_region_input/",
        tract_dir=DATA_DIR + "bed/"
    shell:
        """
        mkdir -p {params.sample_dir}
        mkdir -p {params.region_dir}
        conda activate py36
        python partition_anc_het.py --hits {input.hits} \
            --tracts {params.tract_dir} \
            --chrom_lengths {input.chrom_lengths} \
            --sample_out {params.sample_dir} \
            --region_out {params.region_dir}
        conda deactivate
        """

rule prep_SNP_file:
    output:
        temp(DATA_DIR + "fastqtl_SNP_input/{SNP}.txt")
    params:
        output_dir=DATA_DIR + "fastqtl_SNP_input/"
    shell:
        """
        mkdir -p {params.output_dir}
        echo "{wildcards.SNP}" > {output}
        """

rule anc_het_eQTLs:
    input:
        vcf=DATA_DIR + "genotype_freeze.6a.pass_only.phased.mesa_1331samples.maf01.biallelic.vcf.gz",
        pheno=rules.normalize_expression.output.bed,
        sample_input=DATA_DIR + "fastqtl_anc_het_sample_input/{SNP}.txt",
        SNP_input=rules.prep_SNP_file.output,
        region_input=DATA_DIR + "fastqtl_region_input/{SNP}.txt",
        covariates=rules.combine_covariates.output
    output:
        DATA_DIR + "fastqtl_anc_het_output/{SNP}.txt"
    singularity:
        "gtex_eqtl_V8.sif"
    params:
        output_dir=DATA_DIR + "fastqtl_anc_het_output/"
    shell:
        """
        mkdir -p {params.output_dir}
        CHR=$(cut -f1 {input.region_input})
        REGION_START=$(cut -f2 input.region_input)
        REGION_STOP=$(cut -f3 input.region_input)
        fastQTL --vcf {input.vcf} --bed {input.pheno} \
            --include-samples {input.sample_input} \
            --include-sites {input.SNP_input} \
            --include-covariates {input.covariates} \
            --out {output} \
            --region chr$CHR:$REGION_START-$REGION_STOP
        """

def merge_anc_het_input(wildcards):
    checkpoint_output = checkpoints.partition_anc_het.get(**wildcards).output[0]
    return expand(DATA_DIR + "fastqtl_anc_het_output/{SNP}.txt",
           SNP=glob_wildcards(os.path.join(checkpoint_output, "{SNP}.txt")).SNP)

rule merge_anc_het_output:
    input:
        merge_anc_het_input
    output:
        DATA_DIR + "fastqtl_anc_het_output/merged_anc_het.txt"
    params:
        dir=DATA_DIR + "fastqtl_anc_het_output"
    shell:
        """
        ls -l {params.dir} | wc -l > {output}
        """

