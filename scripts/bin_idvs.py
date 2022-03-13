# Run in py36 environment
import pandas as pd

# Create 2 bins for global ancestry
ind_table = pd.read_csv("protected_data/mesa.sample_covariates.txt", sep='\t')
EA = ind_table.loc[ind_table.race == 0].nwd_id
ganc_quantiles = ind_table.loc[ind_table.race == 1].global_ancestry.describe()
AA_low = ind_table.loc[(ind_table.race == 1) & (ind_table.global_ancestry < ganc_quantiles["50%"])].nwd_id
AA_high = ind_table.loc[ind_table.global_ancestry >= ganc_quantiles["50%"]].nwd_id
pd.concat([AA_low, EA]).to_csv("data/bins/idv_bin_low.txt", index=False)
pd.concat([AA_high, EA]).to_csv("data/bins/idv_bin_high.txt", index=False)

AA_1 = ind_table.loc[(ind_table.race == 1) & (ind_table.global_ancestry < ganc_quantiles["25%"])].nwd_id
AA_2 = ind_table.loc[(ind_table.global_ancestry >= ganc_quantiles["25%"]) & (ind_table.global_ancestry < ganc_quantiles["50%"])].nwd_id
AA_3 = ind_table.loc[(ind_table.global_ancestry >= ganc_quantiles["50%"]) & (ind_table.global_ancestry < ganc_quantiles["75%"])].nwd_id
AA_4 = ind_table.loc[ind_table.global_ancestry >= ganc_quantiles["75%"]].nwd_id
pd.concat([AA_1, EA]).to_csv("data/bins/idv_bin_1.txt", index=False)
pd.concat([AA_2, EA]).to_csv("data/bins/idv_bin_2.txt", index=False)
pd.concat([AA_3, EA]).to_csv("data/bins/idv_bin_3.txt", index=False)
pd.concat([AA_4, EA]).to_csv("data/bins/idv_bin_4.txt", index=False)
