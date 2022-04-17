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
blast_result_fields = ["qseqid", "sacc", "stitle", "staxid", "sseq", "qseq", "length", "pident", "evalue", "bitscore", "score"]
blast_field_string = " ".join(["6"] + blast_result_fields)
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
human_blast_results = "human_blast_matches.tsv"
!blastp -query {peptide_fasta} -db {proteome_dir}/{proteome} -task blastp-short -ungapped -comp_based_stats F -evalue 50 -num_threads 4 -outfmt "6 qseqid sseqid sacc staxid sseq qseq length pident evalue bitscore score" -out {human_blast_results}
```

```python
df_human = pd.read_csv(human_blast_results, sep="\t", names=blast_result_fields).sort_values(by="evalue")
df_human["long_name"] = df_human["stitle"].str.replace(" OS.*", "", regex=True)
df_human["long_name"] = df_human["long_name"].str.replace(" (Fragment)", "", regex=False)
df_human["seq_pos"] = df_human["qseqid"].str.replace("epi_", "", regex=False).astype(int)
df_human["pident"] = df_human["pident"].astype(float)
df_human.head(15)
```

```python
protein_hierarchy = (
    df_human
    .loc[(((df_human["length"] > 6) & (df_human["pident"] == 100)) | (df_human["length"] > 7)),:]
    .loc[df_human["seq_origin"] == "epi_1",:]
    .groupby(["long_name"])
    .agg(dict(evalue=min))
    .sort_values(by="evalue")
    .reset_index()
)
```

```python
sub_df = (
    df_human
    .loc[(((df_human["length"] > 6) & (df_human["pident"] == 100)) | (df_human["length"] > 7)),:]
    .loc[df_human["seq_origin"] == "epi_1",:]
    .groupby(["long_name", "seq_pos"])
    .agg(dict(length=max, pident=max))
    .reset_index()
    .rename(columns=dict(pident="percent_identity"))
)

fig, ax = plt.subplots(figsize=(9.6, 14.4))
ax = sns.scatterplot(
    data=sub_df.set_index("long_name").loc[protein_hierarchy["long_name"].tolist(),:],
    x="seq_pos", y="long_name", hue="percent_identity", size="length",
    sizes=(1, 500), legend="brief", palette="viridis", ax=ax
)
```

```python

```