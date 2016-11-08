# Bioinformatics HW2

### Pairwise Sequence Alignment

### Dependencies
1. R Package [proto](https://cran.r-project.org/web/packages/proto/index.html): for argparse
2. R Package [argparse](https://cran.r-project.org/web/packages/argparse/index.html): for parse argv
3. R Package [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html): for fasta file I/O

### hw2_104753013.R

```
Rscript hw2_104753013.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta
```

usage: hw2_104753013.R [-h] --input INPUT --score SCORE [--aln {global,local}]
                       [--gap_open GAP_OPEN] [--gap_extend GAP_EXTEND]
                       --output OUTPUT [--no_penalty]

Pairwise sequence alignment with gap open penalty and gap extend penalty.

### score.R

```
Rscript score.R --input result.fasta --score PAM250.txt --gap_open -10 --gap_extend -2
```

usage: score.R [-h] --input INPUT --score SCORE [--gap_open GAP_OPEN]
               [--gap_extend GAP_EXTEND]

Calculate pairwise sequence alignment score.
