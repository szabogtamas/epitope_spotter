import pandas as pd

def blast_parser(fn, fields = None):

    if fields is None:
        blast_result_fields = [
            "qseqid", "sacc", "stitle", "staxid", "sseq", "qseq", "length",
            "pident", "evalue", "bitscore", "score", "qstart"
        ]

    df_human = pd.read_csv(human_blast_results, sep="\t", names=blast_result_fields).sort_values(by="evalue")
    df_human["long_name"] = df_human["stitle"].str.replace(" OS=.*", "", regex=True)
    df_human["long_name"] = df_human["long_name"].str.replace(" (Fragment)", "", regex=False)
    df_human["seq_origin"] = df_human["qseqid"].str.replace("(?<=epi_1|epi_2)_\d*", "", regex=True)
    df_human["seq_pos"] = df_human["qseqid"].str.replace("(epi_1|epi_2)_", "", regex=True).astype(int)
    df_human["pident"] = df_human["pident"].astype(float)

    return df