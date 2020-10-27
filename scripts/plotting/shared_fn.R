library(tidyverse)
library(ggpubr)

merge_hits = function(res_dir) {
  # Load estimates of significant eQTLs
  Eur = read_tsv(paste0(res_dir, "hits_reestimation_primary_Eur.txt"))
  Afr = read_tsv(paste0(res_dir, "hits_reestimation_primary_Afr.txt"))
  het = read_tsv(paste0(res_dir, "hits_reestimation_primary_het.txt"))
  val = read_tsv(paste0(res_dir, "hits_reestimation_validation_Eur.txt"))
  
  # Merge significant eQTLs together
  merge_hetAfr = inner_join(het, Afr, by=c("gene", "ID"), suffix=c("_het", "_Afr")) 
  merge_valEur = inner_join(val, Eur, by=c("gene", "ID"), suffix=c("_val", "_Eur")) 
  merged = inner_join(merge_hetAfr, merge_valEur, by=c("gene", "ID"))
  return(merged)
}

compute_plot_cols = function(merged) {
  # Compute desired row quantities
  merged = merged %>% mutate(AfrEur_avg = (effect_Afr + effect_Eur) / 2)
  merged = merged %>% mutate(varAfr = maf_Afr * (1 - maf_Afr), varEur = maf_Eur * (1 - maf_Eur))
  merged = merged %>% mutate(alpha = varAfr / (varAfr + varEur))
  merged = merged %>% mutate(AfrEur_wtavg = alpha * effect_Afr + (1 - alpha) * effect_Eur)
  merged = merged %>% mutate(obs_dif_het = effect_het - effect_Eur,
                             exp_dif_avg = AfrEur_avg - effect_Eur, 
                             exp_dif_wtavg = AfrEur_wtavg - effect_Eur,
                             obs_dif_val = effect_val - effect_Eur)
  return(merged)
}

scale_col = function(x) {
  return((x - mean(x)) / sd(x))
}

one_dim_PCA = function(tib) {
  mean_tib = tib %>% summarize_all(mean)
  sd_tib = tib %>% summarize_all(sd)
  scaled = tib %>% mutate_all(scale_col)
  scaled_pca = prcomp(scaled)
  scaled_loading = scaled_pca$rotation[,1]
  unscaled_loading = scaled_loading * sd_tib
  slope = unscaled_loading[2] / unscaled_loading[1]
  intercept = mean_tib[2] - slope * mean_tib[1]
  return(c(slope = slope, intercept = intercept))
}

set_up = function(dir_path) {
  tib = merge_hits(dir_path)
  tib = compute_plot_cols(tib)
  return(tib)
}
