library(tidyverse)
library(optparse)

# Parse arguments
option_list = list(make_option("--sd_file", type="character", default=NULL, 
			help="std dev file name", metavar="character"),
		   make_option("--afr_effect_file", type="character", default=NULL, 
			help="afr effect sizes file name", metavar="character"),
		   make_option("--eur_effect_file", type="character", default=NULL, 
			help="eur effect sizes file name", metavar="character"),
		   make_option("--afr_scaling", type="double", default=NULL, 
			help="afr scaling value", metavar="character"),
		   make_option("--eur_scaling", type="double", default=NULL, 
			help="eur scaling value", metavar="character"),
		   make_option("--out", type="character", default=NULL, 
              		help="output file name", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Define cis-heritability of gene expression
h2 = .2

# Define scaling used for alternate hypothesis
scaling_Afr = opt$afr_scaling
scaling_Eur = opt$eur_scaling

# Read files into tibbles
sd = read_tsv(opt$sd_file) %>% drop_na(Afr_sd, het_sd, Eur_sd)
Afr_effects = read_tsv(opt$afr_effect_file)
Eur_effects = read_tsv(opt$eur_effect_file)

effects = inner_join(Afr_effects, Eur_effects, by=c("gene", "variant"), suffix=c("_Afr", "_Eur")) %>% inner_join(sd, by=c("gene" = "gene_id"))

# Filter out variants with MAF < .05 in either ancestral population; we can't simulate genotypes in these situations
effects = effects %>% filter(maf_Afr > 0.05, maf_Eur > 0.05) 

genotype_variance = function(p) {
    var = p * (1 - p)
    return(var)
}

# Returns expected effect size for ancestry heterozygous individual under null or alternate 
# hypothesis.
generate_beta_het = function(beta_Afr, beta_Eur, p_Afr, p_Eur, mode) {
    var_Afr = genotype_variance(p_Afr)
    var_Eur = genotype_variance(p_Eur)
    if (mode == 'null') {
        beta = (beta_Afr * var_Afr + beta_Eur * var_Eur) / (var_Afr + var_Eur) 
    } else if (mode == 'alt') {
        beta = (scaling_Afr * beta_Afr * var_Afr + scaling_Eur * beta_Eur * var_Eur) / (scaling_Afr * var_Afr + scaling_Eur * var_Eur)
    } else {
        print('Incorrect mode specified for heterozygous effect size.')
        return(NA)
    }
    return(beta)
}
generate_beta_het_vect = Vectorize(generate_beta_het)

# Accepts p, a vector of length 2 containing the allele frequency of the alternate allele 
# on each haplotype; and n, the number of individuals. Simulates genotypes for all 
# individuals via binomial random sampling.
simulate_geno = function(p, n) {
    hap1 = tibble(hap1 = rbinom(n, 1, p[1]))
    hap2 = tibble(hap2 = rbinom(n, 1, p[2]))
    geno = bind_cols(hap1, hap2)
    return(geno)   
}

# Accepts genotypes, a nx2 table of genotypes; effect, a vector of length 2 containing 
# the effect of the alternate allele on each haplotype; and sd, a vector of length 2
# containing the standard deviation of expression for the given gene on each haplotype. 
# Simulates gene expression for all individuals as the sum of allelic effects, plus a 
# normally distributed error that depends on sd and h2.
simulate_pheno = function(genotypes, effect, sd) {
    effect_sum = rowSums(genotypes * effect)
    n = genotypes %>% nrow
    noise = rnorm(n, 0, sd * sqrt(1 - h2))
    pheno = tibble(pheno = effect_sum + noise)
    return(pheno)
}

normalize = function(phenotypes) {
    return(phenotypes / sd(phenotypes %>% pull(pheno)))
}

group_normalize = function(pheno1, pheno2, pheno3, pheno4) {
    tib = bind_rows(pheno1, pheno2, pheno3, pheno4)
    group_sd = sd(tib %>% pull(pheno))
    return(c(pheno1/group_sd, pheno2/group_sd, pheno3/group_sd, pheno4/group_sd))
}

estimate_effects = function(geno, pheno) {
    tib = bind_cols(geno, pheno)
    tib = tib %>% transmute(X = hap1 + hap2, Y = pheno)
    if (tib %>% pull(X) %>% sum > 0) {
        beta = summary(lm(Y~X, data=tib))$coeff[2,1]
        return(beta)
    }
    return(NA)
}

# Determine ancestry-heterozygous effect sizes under null and alternative hypotheses
effects = effects %>% mutate(beta_het_null = generate_beta_het_vect(slope_Afr, slope_Eur, maf_Afr, maf_Eur, "null"),
                             beta_het_alt = generate_beta_het_vect(slope_Afr, slope_Eur, maf_Afr, maf_Eur, "alt"))

simulate_eqtl_calling = function(p_Afr, p_Eur, beta_Afr, beta_het_null, beta_het_alt, beta_Eur, sd_Afr, sd_het, sd_Eur) {
    # Simulate genotypes for all three ancestries
    Afr_geno = simulate_geno(c(p_Afr, p_Afr), 100)
    het_geno = simulate_geno(c(p_Afr, p_Eur), 50)
    Eur_geno = simulate_geno(c(p_Eur, p_Eur), 100)

    # Simulate phenotypes for all ancestries, under null and alternative hypotheses
    Afr_pheno = simulate_pheno(Afr_geno, c(beta_Afr, beta_Afr), sd_Afr)
    het_null_pheno = simulate_pheno(het_geno, c(beta_het_null, beta_het_null), sd_het)
    het_alt_pheno = simulate_pheno(het_geno, c(beta_het_alt, beta_het_alt), sd_het)
    Eur_pheno = simulate_pheno(Eur_geno, c(beta_Eur, beta_Eur), sd_Eur)
    
    # Normalize all phenotypes individually
    Afr_pheno_indnorm = normalize(Afr_pheno)
    het_null_pheno_indnorm = normalize(het_null_pheno)
    het_alt_pheno_indnorm = normalize(het_alt_pheno)
    Eur_pheno_indnorm = normalize(Eur_pheno)
    
    # Normalize phenotypes together
    group = group_normalize(Afr_pheno, het_null_pheno, het_alt_pheno, Eur_pheno)
    Afr_pheno_grpnorm = group[1]
    het_null_pheno_grpnorm = group[2]
    het_alt_pheno_grpnorm = group[3]
    Eur_pheno_grpnorm = group[4]

    # Estimate effect sizes for all ancestries from simulated data
    ## African effect size
    est_beta_Afr_indnorm = estimate_effects(Afr_geno, Afr_pheno_indnorm)
    est_beta_Afr_grpnorm = estimate_effects(Afr_geno, Afr_pheno_grpnorm)
    
    ## Ancestry-heterozygous effect size under null model
    est_beta_het_null_indnorm = estimate_effects(het_geno, het_null_pheno_indnorm)
    est_beta_het_null_grpnorm = estimate_effects(het_geno, het_null_pheno_grpnorm)

    ## Ancestry-heterozygous effect size under alternative model
    est_beta_het_alt_indnorm = estimate_effects(het_geno, het_alt_pheno_indnorm)
    est_beta_het_alt_grpnorm = estimate_effects(het_geno, het_alt_pheno_grpnorm)

    ## European effect size
    est_beta_Eur_indnorm = estimate_effects(Eur_geno, Eur_pheno_indnorm)
    est_beta_Eur_grpnorm = estimate_effects(Eur_geno, Eur_pheno_grpnorm)

    # Calculate allele frequencies from simulations
    sim_Afr_maf = sum(Afr_geno) / (nrow(Afr_geno) * ncol(Afr_geno))
    sim_Eur_maf = sum(Eur_geno) / (nrow(Eur_geno) * ncol(Eur_geno))
    
    # Calculate expected ancestry-heterozygous effect size under null and
    # alternative hypotheses, using individual and group normalization
    exp_beta_het_indnorm = generate_beta_het(est_beta_Afr_indnorm, est_beta_Eur_indnorm, sim_Afr_maf, sim_Eur_maf, "null")
    exp_beta_het_grpnorm = generate_beta_het(est_beta_Afr_grpnorm, est_beta_Eur_grpnorm, sim_Afr_maf, sim_Eur_maf, "null")

    # Return simulated data
    return(c(exp_indnorm = exp_beta_het_indnorm, exp_grpnorm = exp_beta_het_grpnorm,
             est_null_indnorm = est_beta_het_null_indnorm, est_null_grpnorm = est_beta_het_null_grpnorm,
             est_alt_indnorm = est_beta_het_alt_indnorm, est_alt_grpnorm = est_beta_het_alt_grpnorm))
}

simulation_results = effects %>% group_by(gene, variant) %>% do(data.frame(t(simulate_eqtl_calling(.$maf_Afr, .$maf_Eur, .$slope_Afr, .$beta_het_null, .$beta_het_alt, .$slope_Eur, .$Afr_sd, .$het_sd, .$Eur_sd)))) %>% drop_na()
pval_null_indnorm = t.test(simulation_results$exp_indnorm, simulation_results$est_null_indnorm, paired=TRUE)$p.value
pval_alt_indnorm = t.test(simulation_results$exp_indnorm, simulation_results$est_alt_indnorm, paired=TRUE)$p.value
pval_null_grpnorm = t.test(simulation_results$exp_grpnorm, simulation_results$est_null_grpnorm, paired=TRUE)$p.value
pval_alt_grpnorm = t.test(simulation_results$exp_grpnorm, simulation_results$est_alt_grpnorm, paired=TRUE)$p.value

null_results = tibble(indnorm = pval_null_indnorm, grpnorm = pval_null_grpnorm, hyp="null")
alt_results = tibble(indnorm = pval_alt_indnorm, grpnorm = pval_alt_grpnorm, hyp="alt")

bind_rows(null_results, alt_results) %>% write.table(opt$out, row.names=FALSE, sep='\t')
