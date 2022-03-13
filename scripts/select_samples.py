import pandas as pd
import numpy as np
import argparse

def swap_race(x):
    if x == "Black":
        return("White")
    else:
        return("Black")

def swap_exam(x):
    if x == "1":
        return("5")
    else:
        return("1")

def sort_samples(samp_df, curr_props, counts, prev_dict):
    def update_prop_df(curr_props, counts, id, race, exam, seq_center, sex, prev_dict=None):
        """After choosing an exam for a particular individual, update the 
           dictionary containing covariate proportions accordingly."""
        if prev_dict is not None: # Second or later iteration of algorithm
            # If we decided to swap the exam for the current individual relative
            # to the last iteration of the algorithm, update dictionary 
            # containing covariate proportions; otherwise, do nothing.
            if id not in prev_dict[int(exam)]:
                # Convert covariate proportions to counts
                curr_props[race] = curr_props[race] * counts[race]

                # Update counts based on the exam chosen for the current
                # individual
                curr_props[race, exam, seq_center, sex] += 1
                curr_props[race, swap_exam(exam), seq_center, sex] -= 1

                # Convert covariate counts back to proportions
                curr_props[race] = curr_props[race] / counts[race]
        else: # First iteration of algorithm
            # Convert covariate proportions to counts
            curr_props[race] = curr_props[race] * counts[race]

            # Update counts based on the exam chosen for the current individual
            curr_props[race, exam, seq_center, sex] += 1
            counts[race] += 1

            # Convert covariate counts back to proportions
            curr_props[race] = curr_props[race] / counts[race]

    # Initialize dictionary to store which exam was chosen for each individual
    # with two samples
    exam_choosing_dict = {1: [], 5: []}

    # Create list of all individual IDs we are iterating over (all individuals
    # with two samples)
    samp_df_index = list(pd.concat([samp_df[samp_df.race == "Black"].groupby("nwd_id").size().reset_index(), \
                                    samp_df[samp_df.race == "White"].groupby("nwd_id").size().reset_index()]).sort_index().nwd_id)

    for id in samp_df_index:
        ind = samp_df[samp_df.nwd_id == id]

        # Fetch race and sex of current individual
        race = ind.race.unique()[0]
        sex = ind.sex.unique()[0]

        # Fetch sequencing center for exam 1 and 5 for current individual 
        exam1_center = ind[ind.exam == "1"].seq_center.iloc[0]
        exam5_center = ind[ind.exam == "5"].seq_center.iloc[0]

        # For each exam, fetch covariate proportions for the bin corresponding 
        # to the current individual for the same race (curr) and the opposite
        # race (swap)
        prop1_curr = curr_props[race, "1", exam1_center, sex]
        prop1_swap = curr_props[swap_race(race), "1", exam1_center, sex]
        prop5_curr = curr_props[race, "5", exam1_center, sex]
        prop5_swap = curr_props[swap_race(race), "5", exam5_center, sex]

        # If the difference in covariate proportions is smaller for exam 1, use
        # exam 1 for the current individual, and vice versa. This ultimately
        # means that covariate proportions for each bin should get more similar
        # between the two races.
        diff1 = prop1_curr - prop1_swap
        diff5 = prop5_curr - prop5_swap
        if diff1 < diff5:
            exam_choosing_dict[1].append(id)
            update_prop_df(curr_props, counts, id, race, "1", exam1_center, sex, prev_dict)
        else:
            exam_choosing_dict[5].append(id)
            update_prop_df(curr_props, counts, id, race, "5", exam5_center, sex, prev_dict)

    return(exam_choosing_dict, curr_props, counts)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_data')
    parser.add_argument('--exclusion_list')
    parser.add_argument('--out')
    args = parser.parse_args()
    
    samples = pd.read_csv(args.sample_data, sep='\t')
    
    # Exclude individuals that failed ancestry determination QC
    with open(args.exclusion_list, 'r') as f:
        exclusion_list = [x.strip() for x in f.readlines()]
    samples = samples[-samples.nwd_id.isin(exclusion_list)]
    
    # Recode "Black or African American" in race field as "Black" for simplicity
    samples.loc[samples.race == "Black or African American", "race"] = "Black"
    
    # Select high-quality PBMC RNASeq samples from Black/White genotyped individuals
    samples = samples[(samples.analysis_freeze == True) & 
                      (samples.has_genotype == True) & 
                      (samples.sample_type == "PBMC") & 
                      ((samples.race == "Black") | (samples.race == "White"))]
    
    # Certain individuals in MESA have two gene expression samples (exam 1 and
    # exam 5). We use an iterative greedy approach to choose which exam to use
    # for individuals with two samples, with the ultimate goal of having
    # identical covariate proportions between African-Americans and 
    # European-Americans. By "covariate proportions", we mean the proportion of 
    # individuals of a certain race in each covariate bin (i.e. each 
    # combination of exam/sequencing center/sex metadata).

    # Create table to count number of samples per individual
    ind_sample_counts = samples.groupby("nwd_id").size()

    # Store covariate proportions for individuals with one sample and individuals
    # with two samples. 
    one_samp_index = list(ind_sample_counts[ind_sample_counts == 1].index)
    one_samp = samples.set_index("nwd_id").loc[one_samp_index,:].reset_index()
    one_samp_props = one_samp.groupby(["race", "exam", "seq_center", "sex"]).size().divide(one_samp.groupby(["race"]).size())
    two_samp_index = list(ind_sample_counts[ind_sample_counts == 2].index)
    two_samp = samples.set_index("nwd_id").loc[two_samp_index,:].reset_index()
    two_samp_props = two_samp.groupby(["race", "exam", "seq_center", "sex"]).size().divide(two_samp.groupby(["race"]).size())
    
    # Create dictionary to store counts of individuals by race
    counts = {"White": one_samp[one_samp.race == "White"].shape[0], 
              "Black": one_samp[one_samp.race == "Black"].shape[0]}

    # Initialize table to store final covariate proportions
    samp_props = one_samp_props.copy()

    # Initialize list to store correlation between covariate proportions for
    # European-American and African-American individuals.
    corr_lst = [np.corrcoef(samp_props["Black"], samp_props["White"])[0,1]]
    
    # Perform 10 iterations of sample-choosing algorithm. (Could keep going,
    # but in practice, this is good enough.)
    for i in range(10):
        if i == 0:
            answer_dict, samp_props, counts = \
                sort_samples(two_samp, samp_props, counts, None)
        else:
            answer_dict, samp_props, counts = \
                sort_samples(two_samp, samp_props, counts, answer_dict)
        corr_lst.append(np.corrcoef(samp_props["Black"], samp_props["White"])[0,1])
    
    # Combine correct samples for all individuals based on sample-choosing
    # algorithm.
    two_samp_exam_1 = two_samp[two_samp.exam == "1"].set_index("nwd_id").loc[answer_dict[1],:].reset_index()
    two_samp_exam_5 = two_samp[two_samp.exam == "5"].set_index("nwd_id").loc[answer_dict[5],:].reset_index()
    all_samp = pd.concat([two_samp_exam_1, two_samp_exam_5, one_samp])
    all_samp_props = all_samp.groupby(["race", "seq_center", "exam", "sex"]).size().divide(all_samp.groupby(["race"]).size())
    
    # Write TOR ID and NWDID to file
    all_samp[all_samp.race == "Black"][['tor_id', 'nwd_id']].to_csv(args.out + "_Afr.txt", sep='\t', index=False)
    all_samp[all_samp.race == "White"][['tor_id', 'nwd_id']].to_csv(args.out + "_Eur.txt", sep='\t', index=False)
