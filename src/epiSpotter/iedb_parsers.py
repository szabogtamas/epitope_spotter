import pandas as pd

def assay_table_parser(fn):
    df = pd.read_csv(fn)
    return df
