import os, json, pandas as pd

def load_data(data_dir, expr_file, mut_file, labels_file, pathways_file):
    expr = pd.read_csv(os.path.join(data_dir, expr_file), index_col=0)
    mut  = pd.read_csv(os.path.join(data_dir, mut_file), index_col=0)
    labels = pd.read_csv(os.path.join(data_dir, labels_file))
    with open(os.path.join(data_dir, pathways_file),'r') as f:
        pathways = json.load(f)
    # if expression is samples x genes with first col sample_id (from R), handle that:
    if 'sample_id' in expr.columns:
        expr = expr.set_index('sample_id')
    if 'sample_id' in mut.columns:
        mut = mut.set_index('sample_id')
    expr = expr.loc[labels['sample_id']]
    mut = mut.loc[labels['sample_id']]
    y = labels['phenotype'].values
    return expr, mut, y, pathways
