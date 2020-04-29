sim_res %>% filter(norm=="ind") %>% select(-norm) %>% summarize_all(mean)
sim_res %>% filter(norm=="ind") %>% select(-norm) %>% summarize_all(sd)
sim_res %>% filter(norm=="grp") %>% select(-norm) %>% summarize_all(mean)
sim_res %>% filter(norm=="grp") %>% select(-norm) %>% summarize_all(sd)

n_trial = nrow(sim_res) / 2
type_I_error_ind = sim_res %>% filter(norm=="ind") %>% filter(null_coeff > 1) %>% nrow / n_trial
type_I_error_grp = sim_res %>% filter(norm=="grp") %>% filter(null_coeff > 1) %>% nrow / n_trial
power_ind = sim_res %>% filter(norm=="ind") %>% filter(alt_coeff > 1) %>% nrow / n_trial
power_grp = sim_res %>% filter(norm=="grp") %>% filter(alt_coeff > 1) %>% nrow / n_trial
