import pandas as pd

def aggregate_expression_to_pathways(expr_df: pd.DataFrame, pathways: dict) -> pd.DataFrame:
    out = {}
    for pwy, genes in pathways.items():
        genes_in = [g for g in genes if g in expr_df.columns]
        if not genes_in: 
            continue
        out[pwy] = expr_df[genes_in].mean(axis=1)
    return pd.DataFrame(out, index=expr_df.index)

def combine_modalities(pwy_expr: pd.DataFrame, mutations: pd.DataFrame) -> pd.DataFrame:
    idx = pwy_expr.index.intersection(mutations.index)
    return pd.concat([pwy_expr.loc[idx], mutations.loc[idx]], axis=1)
