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
sys.path.append("../src")  # developmental hack, to load the local version of the module
%load_ext autoreload
%autoreload 2

import epiSpotter as epi
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
def aa_mapper(row, N):
    N_pre = row["seq_pos"] + row["qstart"] - 2
    N_post = N - (N_pre + len(row["sseq"]))
    p = "".join(["."]*N_pre) + row["sseq"] + "".join(["."]*N_post)
    return p

top_match_df = (
    df_human
    .loc[(((df_human["length"] > 6) & (df_human["pident"] == 100)) | (df_human["length"] > 7)),:]
    .loc[df_human["seq_origin"] == "epi_1",:]
)

top_match_df["aligned_seq"] = top_match_df.apply(aa_mapper, N=len(epitope_collection)+epi_len, axis=1)

protein_hierarchy = (
    top_match_df
    .groupby(["long_name"])
    .agg(dict(evalue=min))
    .sort_values(by="evalue")
    .reset_index()
)

top_match_df["name_orders"] = pd.Categorical(top_match_df["long_name"], categories=protein_hierarchy["long_name"].tolist(), ordered=True)
```

```python
sub_df = (
    top_match_df
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

ax.set_xticks(range(1, max(sub_df["seq_pos"])))
ax.set_xlabel("Peptide position")
ax.set_ylabel("")
ax.grid()
ax.set_title("Peptide positions with matches in the human proteome")
ax.legend(labelspacing=2, bbox_to_anchor=(1.04,1), borderpad=2)

fig.savefig(os.path.join(fig_dir, "human_blast_matches.pdf"), bbox_inches="tight")
```

```python
sub_df = (
    top_match_df
    .loc[:,["aligned_seq", "name_orders"]]
    .drop_duplicates()
    .sort_values(by="name_orders")
    .reset_index()
)

N = sub_df.shape[0]
fig, ax = plt.subplots(figsize=(9.6, 14.4))

for index, row in sub_df.iterrows():
    yloc = 3*(N-index)
    for i, l in enumerate(row["aligned_seq"]):
        xloc = ((i+1)*0.02)+0.5
        if SEQ[i] == l:
            ax.text(xloc, yloc, l, family="monospace", ha="left", bbox=dict(color="#DCDCDC", pad=0.6))
        else:
            ax.text(xloc, yloc, l, family="monospace", ha="left")
    ax.text(0.5, yloc, row["name_orders"], ha="right")

yloc = 3*(N+1)
for i, l in enumerate(SEQ):
    xloc = ((i+1)*0.02)+0.5
    ax.text(xloc, yloc, l, family="monospace", ha="left", bbox=dict(color="#DCDCDC", pad=0.5))
ax.text(0.5, yloc, "Reference peptide", ha="right")

ax.set_xlim(0, (len(SEQ)+2)*0.04)
ax.set_ylim(-1, 3*(N+2))
ax.set_xlabel("")
ax.set_ylabel("")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.set_frame_on(False)

fig.savefig(os.path.join(fig_dir, "human_blast_alignments.pdf"), bbox_inches="tight")
```

```python
def netmhc_i_parser(fn):
    molten_rows = []
    with open(fn, "r") as f:
        lines = f.read().split("\n")
        alleles = lines[0].split("\t")[3::3]
        N = len(alleles) - 1
        for line in lines[2:-1]:
            r = line.split("\t")
            rowids = r[:3]
            molten_rows += [[alleles[i]] + rowids + r[i*3+3:i*3+6] for i in range(0, N)]
    df = pd.DataFrame.from_records(molten_rows, columns=["Allele", "Position", "Peptide", "ID", "nM", "Rank", "Core"])
    return df
```

```python
def netmhc_ii_parser(fn):
    with open(fn, "r") as f:
        df = (
            pd.DataFrame
            .from_records([[x] for x in f.read().split("\n")], columns=["Raw"])
        )
        df = df.loc[df["Raw"] != "",:]
        df = df.loc[~df["Raw"].str.contains("^#|----------|Allele: "),:]
        df = df.iloc[2:]
        df["Raw"] = df["Raw"].str.replace("\s+", " ", regex=True)
        df = df["Raw"].str.split(" ", expand=True)
        df.columns = df.iloc[0]
        df = df.iloc[1:]
    return df
```

```python
peptides_file = "tmp/peptides.fa"
```

```python
!/usr/cbs/packages/netMHC-4.0/netMHC {peptides_file}
```

```python
!/usr/cbs/packages/netMHCII-2.3 {peptides_file}
```