import os, argparse, joblib, numpy as np, pandas as pd, shap
from .data_ingest import load_data
from .features import aggregate_expression_to_pathways, combine_modalities
import matplotlib.pyplot as plt

def main(args):
    expr, mut, y, pathways = load_data(args.data_dir, args.expr_file, args.mut_file, args.labels_file, args.pathways_file)
    Xp = aggregate_expression_to_pathways(expr, pathways)
    X  = combine_modalities(Xp, mut)
    # bundle = joblib.load(args.model_path)
    # model = bundle["model"]
    # background = shap.sample(X, min(50, len(X)), random_state=42)
    # explainer = shap.Explainer(model.predict_proba, background)
    # sv = explainer(X)
    # os.makedirs(args.out_dir, exist_ok=True)
    # shap.summary_plot(sv[:,:,1], X, show=False)
    # plt.tight_layout(); plt.savefig(os.path.join(args.out_dir,"shap_summary.png"), dpi=180); plt.close()
    # mean_abs = np.abs(sv.values[:,:,1]).mean(axis=0)
    # top = pd.Series(mean_abs, index=X.columns).sort_values(ascending=False).head(30)
    # top.to_csv(os.path.join(args.out_dir,"top_features.csv"))

    # --- align to training columns ---
    bundle = joblib.load(args.model_path)
    pipe = bundle["model"]
    train_feats = bundle.get("features")
    if train_feats is not None:
        X = X.reindex(columns=train_feats, fill_value=0)

    # --- transform like training (scaler -> clf) ---
    scaler = getattr(pipe, "named_steps", {}).get("scaler", None)
    clf    = getattr(pipe, "named_steps", {}).get("clf", None)
    if clf is None:
        raise RuntimeError("Expected a Pipeline with a 'clf' (XGBClassifier).")
    import numpy as np
    X_trans = scaler.transform(X) if scaler is not None else X.values

    # --- SHAP: provide a background and use interventional + probability ---
    import shap
    bg = shap.sample(X_trans, min(200, len(X_trans)), random_state=42)  # background in transformed space
    explainer = shap.TreeExplainer(
        clf,
        data=bg,
        feature_perturbation="interventional",
        model_output="probability"
    )

    # subsample rows for plotting to keep runtime reasonable
    n_plot = min(500, len(X_trans))
    idx = np.random.RandomState(42).choice(len(X_trans), size=n_plot, replace=False)
    X_plot = X.iloc[idx]          # for feature names
    X_trans_plot = X_trans[idx]

    sv = explainer.shap_values(X_trans_plot)

    # --- save plot & top features ---
    import matplotlib.pyplot as plt, pandas as pd, os
    shap.summary_plot(sv, X_plot, show=False)
    plt.tight_layout(); plt.savefig(os.path.join(args.out_dir, "shap_summary.png"), dpi=180); plt.close()
    mean_abs = np.abs(sv).mean(axis=0)
    pd.Series(mean_abs, index=X_plot.columns).sort_values(ascending=False).head(30)\
    .to_csv(os.path.join(args.out_dir, "top_features.csv"))

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--data_dir", default="data")
    p.add_argument("--expr_file", default="toy_expression.csv")
    p.add_argument("--mut_file", default="toy_mutations.csv")
    p.add_argument("--labels_file", default="toy_labels.csv")
    p.add_argument("--pathways_file", default="toy_pathways.json")
    p.add_argument("--model_path", default="reports/best_model.joblib")
    p.add_argument("--out_dir", default="reports")
    main(p.parse_args())
