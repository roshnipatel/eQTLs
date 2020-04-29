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

# Define scaling used for alternate hypothesis
scaling_Afr = opt$afr_scaling
scaling_Eur = opt$eur_scaling

# Define probability that variant is causal
prob_causal = .7

# Read files into tibbles
sd = read_tsv(opt$sd_file) %>% drop_na(Afr_sd, het_sd, Eur_sd)
Afr_effects = read_tsv(opt$afr_effect_file)
Eur_effects = read_tsv(opt$eur_effect_file)

effects = inner_join(Afr_effects, Eur_effects, by=c("Gene", "ID"), suffix=c("_Afr", "_Eur")) %>% inner_join(sd, by=c("Gene" = "gene_id"))

# Filter out variants with MAF < .05 in either ancestral population; we can't simulate genotypes in these situations
effects = effects %>% filter(maf_Afr > 0.05, maf_Eur > 0.05) 

genotype_variance = function(p) {
    var = p * (1 - p)
    return(var)
}

generate_beta = function(n, p, sd) {
    prior = rbinom(n, 1, p) 
    betas = prior * rnorm(n, 0, sd)
    return(betas)
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
simulate_pheno = function(genotypes, effect, exp_sd, maf) {
    effect_sum = rowSums(genotypes * c(effect, effect))
    n = genotypes %>% nrow
    gen_var = effect ^ 2 * sum(genotype_variance(maf[1]), genotype_variance(maf[2]))
    env_var = (exp_sd ^ 2) - gen_var
    noise = rnorm(n, 0, sqrt(env_var))
    pheno = tibble(pheno = effect_sum + noise)
    return(pheno)
}

normalize = function(phenotypes) {
    mean_pheno = phenotypes %>% summarize(mean(pheno)) %>% first
    sd_pheno = phenotypes %>% summarize(sd(pheno)) %>% first
    return((phenotypes - mean_pheno) / sd_pheno)
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
effects = effects %>% mutate(beta_het_null = generate_beta_het_vect(effect_Afr, effect_Eur, maf_Afr, maf_Eur, "null"),
                             beta_het_alt = generate_beta_het_vect(effect_Afr, effect_Eur, maf_Afr, maf_Eur, "alt"))

simulate_eqtl_calling = function(p_Afr, p_Eur, beta_Afr, beta_het_null, beta_het_alt, beta_Eur, sd_Afr, sd_het, sd_Eur) {
    # Simulate genotypes for all three ancestries
    Afr_geno = simulate_geno(c(p_Afr, p_Afr), 100)
    het_geno = simulate_geno(c(p_Afr, p_Eur), 50)
    Eur_geno = simulate_geno(c(p_Eur, p_Eur), 100)

    # Simulate phenotypes for all ancestries, under null and alternative hypotheses
    Afr_pheno = simulate_pheno(Afr_geno, beta_Afr, sd_Afr, c(p_Afr, p_Afr))
    het_null_pheno = simulate_pheno(het_geno, beta_het_null, sd_het, c(p_Afr, p_Eur))
    het_alt_pheno = simulate_pheno(het_geno, beta_het_alt, sd_het, c(p_Afr, p_Eur))
    Eur_pheno = simulate_pheno(Eur_geno, beta_Eur, sd_Eur, c(p_Eur, p_Eur))
    
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
    est_beta_Afr = estimate_effects(Afr_geno, Afr_pheno)
    
    ## Ancestry-heterozygous effect size under null model
    est_beta_het_null_indnorm = estimate_effects(het_geno, het_null_pheno_indnorm)
    est_beta_het_null = estimate_effects(het_geno, het_null_pheno)

    ## Ancestry-heterozygous effect size under alternative model
    est_beta_het_alt_indnorm = estimate_effects(het_geno, het_alt_pheno_indnorm)
    est_beta_het_alt = estimate_effects(het_geno, het_alt_pheno)

    ## European effect size
    est_beta_Eur_indnorm = estimate_effects(Eur_geno, Eur_pheno_indnorm)
    est_beta_Eur = estimate_effects(Eur_geno, Eur_pheno)

    # Calculate allele frequencies from simulations
    sim_Afr_maf = sum(Afr_geno) / (nrow(Afr_geno) * ncol(Afr_geno))
    sim_Eur_maf = sum(Eur_geno) / (nrow(Eur_geno) * ncol(Eur_geno))
    
    # Calculate expected ancestry-heterozygous effect size under null and
    # alternative hypotheses, using individual and group normalization
    exp_beta_het_indnorm = generate_beta_het(est_beta_Afr_indnorm, est_beta_Eur_indnorm, sim_Afr_maf, sim_Eur_maf, "null")
    exp_beta_het = generate_beta_het(est_beta_Afr, est_beta_Eur, sim_Afr_maf, sim_Eur_maf, "null")

    # Return simulated data
    return(c(est_Afr_ind = est_beta_Afr_indnorm,
             est_Afr = est_beta_Afr,
             est_Eur_ind = est_beta_Eur_indnorm,
             est_Eur = est_beta_Eur,
             est_het_null_ind = est_beta_het_null_indnorm,
             est_het_null = est_beta_het_null,
             est_het_alt_ind = est_beta_het_alt_indnorm,
             est_het_alt = est_beta_het_alt,
             exp_indnorm_dif = exp_beta_het_indnorm - est_beta_Eur_indnorm, 
             exp_dif = exp_beta_het - est_beta_Eur,
             est_null_indnorm_dif = est_beta_het_null_indnorm - est_beta_Eur_indnorm,
             est_null_dif = est_beta_het_null - est_beta_Eur,
             est_alt_indnorm_dif = est_beta_het_alt_indnorm - est_beta_Eur_indnorm, 
             est_alt_dif = est_beta_het_alt - est_beta_Eur))
}

one_dim_PCA = function(x, y) {
    mean_x = mean(x)
    mean_y = mean(y)
    sd_x = sd(x)
    sd_y = sd(y)
    scaled = tibble(x = (x - mean_x) / sd_x, y = (y - mean_y) / sd_y)
    scaled_pca = prcomp(scaled)
    scaled_loading = scaled_pca$rotation[,1]
    unscaled_loading = scaled_loading * c(sd_x, sd_y)
    slope = unscaled_loading[2] / unscaled_loading[1]
    return(slope)
}

regress_betas = function(x, y) {
    model = lm(y ~ x)
    res = coef(summary(model))
    coeff = res["x", "Estimate"]
    se = res["x", "Std. Error"]
    return(c(coeff, se))
}

simulation_results = effects %>% group_by(Gene, ID) %>% do(data.frame(t(simulate_eqtl_calling(.$maf_Afr, .$maf_Eur, .$effect_Afr, .$beta_het_null, .$beta_het_alt, .$effect_Eur, .$Afr_sd, .$het_sd, .$Eur_sd)))) %>% drop_na()
res_null_indnorm_dif_OLS = regress_betas(simulation_results$exp_indnorm_dif, simulation_results$est_null_indnorm_dif)
res_null_dif_OLS = regress_betas(simulation_results$exp_dif, simulation_results$est_null_dif)
res_alt_indnorm_dif_OLS = regress_betas(simulation_results$exp_indnorm_dif, simulation_results$est_alt_indnorm_dif)
res_alt_dif_OLS = regress_betas(simulation_results$exp_dif, simulation_results$est_alt_dif)
res_null_indnorm_dif_1dPCA = one_dim_PCA(simulation_results$exp_indnorm_dif, simulation_results$est_null_indnorm_dif)
res_null_dif_1dPCA = one_dim_PCA(simulation_results$exp_dif, simulation_results$est_null_dif)
res_alt_indnorm_dif_1dPCA = one_dim_PCA(simulation_results$exp_indnorm_dif, simulation_results$est_alt_indnorm_dif)
res_alt_dif_1dPCA = one_dim_PCA(simulation_results$exp_dif, simulation_results$est_alt_dif)

ind_results = tibble(null_coeff_OLS = res_null_indnorm_dif_OLS[1], null_se = res_null_indnorm_dif_OLS[2], alt_coeff_OLS = res_alt_indnorm_dif_OLS[1], alt_se = res_alt_indnorm_dif_OLS[2], norm="ind", null_coeff_1dPCA = res_null_indnorm_dif_1dPCA, alt_coeff_1dPCA = res_alt_indnorm_dif_1dPCA)
results = tibble(null_coeff_OLS = res_null_dif_OLS[1], null_se = res_null_dif_OLS[2], alt_coeff_OLS = res_alt_dif_OLS[1], alt_se = res_alt_dif_OLS[2], norm="none", null_coeff_1dPCA = res_null_dif_1dPCA, alt_coeff_1dPCA = res_alt_dif_1dPCA)

# simulation_results %>% write.table(opt$out, row.names=FALSE, sep='\t')
bind_rows(ind_results, results) %>% write.table(opt$out, row.names=FALSE, sep='\t')
