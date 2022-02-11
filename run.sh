
snakemake \
  -j 1 \
  --latency-wait 300 \
  --use-conda \
  --snakefile 1_coreGenome \
  --cluster "qsub -pe parallel {threads}"
