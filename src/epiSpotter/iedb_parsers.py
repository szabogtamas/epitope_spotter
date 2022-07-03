import pandas as pd

def assay_table_parser(fn):
    df = pd.read_csv(fn, skiprows=1)
    return df
