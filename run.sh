
snakemake \
  -j 1 \
  --latency-wait 1000 \
  --use-conda \
  --snakefile 1_coreGenome \
  --cluster "qsub -pe parallel {threads}"
