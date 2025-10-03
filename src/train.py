# import os, argparse, joblib, numpy as np
# from sklearn.model_selection import StratifiedKFold
# from sklearn.preprocessing import StandardScaler
# from sklearn.pipeline import Pipeline
# from sklearn.metrics import roc_auc_score, f1_score
# from xgboost import XGBClassifier
# from .data_ingest import load_data
# from .features import aggregate_expression_to_pathways, combine_modalities

# def main(args):
#     expr, mut, y, pathways = load_data(args.data_dir, args.expr_file, args.mut_file, args.labels_file, args.pathways_file)
#     Xp = aggregate_expression_to_pathways(expr, pathways)
#     X  = combine_modalities(Xp, mut)
#     pipe = Pipeline([
#         ("scaler", StandardScaler(with_mean=False)),
#         ("clf", XGBClassifier(n_estimators=200, max_depth=4, learning_rate=0.1,
#                               subsample=0.9, colsample_bytree=0.8, reg_lambda=1.0,
#                               random_state=42, n_jobs=-1, eval_metric='logloss'))
#     ])
#     cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
#     aucs, f1s = [], []
#     for tr, te in cv.split(X, y):
#         pipe.fit(X.iloc[tr], y[tr])
#         proba = pipe.predict_proba(X.iloc[te])[:,1]
#         pred  = (proba>=0.5).astype(int)
#         aucs.append(roc_auc_score(y[te], proba))
#         f1s.append(f1_score(y[te], pred))
#     os.makedirs(args.out_dir, exist_ok=True)
#     with open(os.path.join(args.out_dir,"cv_metrics.txt"),"w") as f:
#         f.write(f"ROC AUC mean±sd: {np.mean(aucs):.3f} ± {np.std(aucs):.3f}\n")
#         f.write(f"F1 mean±sd: {np.mean(f1s):.3f} ± {np.std(f1s):.3f}\n")
#     pipe.fit(X, y)
#     import pandas as pd
#     joblib.dump({"model": pipe, "features": list(X.columns)}, os.path.join(args.out_dir,"best_model.joblib"))

# if __name__ == "__main__":
#     p = argparse.ArgumentParser()
#     p.add_argument("--data_dir", default="data")
#     p.add_argument("--expr_file", default="toy_expression.csv")
#     p.add_argument("--mut_file", default="toy_mutations.csv")
#     p.add_argument("--labels_file", default="toy_labels.csv")
#     p.add_argument("--pathways_file", default="toy_pathways.json")
#     p.add_argument("--out_dir", default="reports")
#     main(p.parse_args())

import os, argparse, joblib, numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
try:
    from sklearn.model_selection import StratifiedGroupKFold  # sklearn >=1.1
    HAVE_SGK = True
except Exception:
    HAVE_SGK = False
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score, f1_score
from xgboost import XGBClassifier

from .data_ingest import load_data
from .features import aggregate_expression_to_pathways, combine_modalities

def patient_id_from_sample(sample_id: str) -> str:
    # TCGA patient = first 12 chars (e.g., TCGA-AB-1234)
    return sample_id[:12]

def main(args):
    expr, mut, y, pathways = load_data(args.data_dir, args.expr_file, args.mut_file,
                                       args.labels_file, args.pathways_file)
    # build pathway features
    Xp = aggregate_expression_to_pathways(expr, pathways)

    # optional modality ablation
    if args.modality == "expr":
        X = Xp.copy()
    elif args.modality == "mut":
        X = mut.copy()
    else:
        X = combine_modalities(Xp, mut)

    # align rows with labels order (already done in load_data, but keep explicit)
    labels = pd.read_csv(os.path.join(args.data_dir, args.labels_file))
    X = X.loc[labels["sample_id"]]
    y = labels["phenotype"].values

    # groups = patient ids for grouped CV
    groups = labels["sample_id"].apply(patient_id_from_sample).values

    # optional permutation test
    if args.permute_labels:
        rng = np.random.RandomState(42)
        y = rng.permutation(y)

    pipe = Pipeline([
        ("scaler", StandardScaler(with_mean=False)),
        ("clf", XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.08,
                              subsample=0.9, colsample_bytree=0.8, reg_lambda=1.0,
                              random_state=42, n_jobs=-1, eval_metric='logloss'))
    ])

    # grouped stratified CV if available; else fallback to standard stratified CV
    if HAVE_SGK:
        cv = StratifiedGroupKFold(n_splits=args.folds, shuffle=True, random_state=42)
        splits = cv.split(X, y, groups=groups)
    else:
        cv = StratifiedKFold(n_splits=args.folds, shuffle=True, random_state=42)
        splits = cv.split(X, y)

    aucs, f1s = [], []
    for tr, te in splits:
        Xtr, Xte = X.iloc[tr], X.iloc[te]
        ytr, yte = y[tr], y[te]
        pipe.fit(Xtr, ytr)
        proba = pipe.predict_proba(Xte)[:, 1]
        pred = (proba >= 0.5).astype(int)
        aucs.append(roc_auc_score(yte, proba))
        f1s.append(f1_score(yte, pred))

    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, "cv_metrics_permute.txt"), "w") as f:
        f.write(f"ROC AUC mean±sd: {np.mean(aucs):.3f} ± {np.std(aucs):.3f}\n")
        f.write(f"F1 mean±sd: {np.mean(f1s):.3f} ± {np.std(f1s):.3f}\n")
        if args.permute_labels:
            f.write("(permutation test enabled)\n")
        f.write(f"Modality: {args.modality}\n")
        f.write(f"CV: {'StratifiedGroupKFold' if HAVE_SGK else 'StratifiedKFold'}\n")

    # final fit on all data (for explanation)
    pipe.fit(X, y)
    joblib.dump({"model": pipe, "features": list(X.columns)},
                os.path.join(args.out_dir, "best_model.joblib"))

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--data_dir", default="data")
    p.add_argument("--expr_file", default="toy_expression.csv")
    p.add_argument("--mut_file", default="toy_mutations.csv")
    p.add_argument("--labels_file", default="toy_labels.csv")
    p.add_argument("--pathways_file", default="toy_pathways.json")
    p.add_argument("--out_dir", default="reports")
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--permute_labels", action="store_true",
                   help="Shuffle labels to detect leakage (AUC should ~0.5).")
    p.add_argument("--modality", choices=["both","expr","mut"], default="both",
                   help="Ablation: use expression only, mutations only, or both.")
    main(p.parse_args())
