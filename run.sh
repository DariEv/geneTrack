
snakemake \
  --use-conda \
  --snakefile 1_coreGenome \
  --cluster "qsub -S /bin/bash
            -pe parallel={resources.threads}""
