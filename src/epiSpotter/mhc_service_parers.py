import pandas as pd

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
        df = df.loc[~(df["Allele"] == "Allele"),:]
    return df