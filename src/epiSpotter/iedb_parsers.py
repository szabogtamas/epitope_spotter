import pandas as pd

def assay_table_parser(fn):
    df = pd.read_csv(fn, skiprows=1)
    return df

def mhc_table_parser(fn):
    df = pd.read_csv(fn, skiprows=1)
    col_subset = [
        "Description", "Parent Protein Accession", "Parent Protein", "Parent Species",
        "Host ID", "Allele Name", "MHC allele class", "Qualitative Measure", "Units"
    ]
    return df.loc[:, col_subset]
