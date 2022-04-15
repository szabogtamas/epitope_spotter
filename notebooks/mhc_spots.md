---
jupyter:
  jupytext:
    formats: md,ipynb
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.0
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

<!-- #region tags=[] -->
# Parse MHC binding predictions for a set of peptides

## Dependencies
<!-- #endregion -->

```python
### Tools to be used
import matplotlib
from matplotlib import pyplot as plt

import numpy as np
import pandas as pd
```

```python
### Define an antigen sequence

antigen_sequence = "WQERRTASAADFAAAKALAMAMA"
```

```python
### Create 8-mer slices

slices = []

for i in range(8, len(antigen_sequence)):
    slices.append(antigen_sequence[i-8:i])
```

```python
# Create BLAST DB

proteome = "human_proteome.fasta"
!gzip -dk {proteome_dir}/{proteome}.gz
!makeblastdb -in {proteome_dir}/{proteome} -parse_seqids -dbtype prot
```

```python
epi_len = 8
peptide_fasta = "query.fas"
s = [">epi_" + str(i+1) + "\n" + SEQ[i:i+epi_len] +"\n" for i in range(0, len(SEQ)-epi_len+1)]
sm = [">bepi_" + str(i+1) + "\n" + SEQ_MINOR[i:i+epi_len] +"\n" for i in range(0, len(SEQ_MINOR)-epi_len+1)]
with open(peptide_fasta, "w") as file:
    file.write("\n".join(s+sm))
```

```python
proteome = "human_proteome.fasta"
out_file = "human_blast_matches.tsv"
!blastp -query {peptide_fasta} -db {proteome_dir}/{proteome} -task blastp-short -ungapped -comp_based_stats F -evalue 50 -num_threads 4 -outfmt "6 qseqid sseqid sacc staxid sseq qseq length pident evalue bitscore score" -out {out_file}
```