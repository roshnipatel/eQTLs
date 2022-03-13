#!/usr/bin/bash

# Run with scripts/slurm_submission/snakemake_submission.sh if re-running
# all eQTL calling; else okay to run in interactive session with sdev.

snakemake --rerun-incomplete --keep-going -j 100 \
        --latency-wait 60 --wait-for-files \
        --cluster-config scripts/sm_slurm_config.json \
        --use-singularity \
        -s SNP_ascertainment_Snakefile \
        --cluster "sbatch -p {cluster.queue} \
                        -t {cluster.time} \
                        --ntasks-per-node={cluster.tasks} \
                        --job-name={cluster.name} \
                        -o {cluster.output} \
                        -e {cluster.error} \
                        --nodes={cluster.nodes} \
                        --mem={cluster.memory}"
