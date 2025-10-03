import os, joblib, pandas as pd, streamlit as st
st.title("Multi-Omics Demo")
st.write("Train via CLI, then explore. This is a minimal viewer.")
model_path = "reports/best_model.joblib"
if not os.path.exists(model_path):
    st.warning("Train a model first (see README).")
else:
    bundle = joblib.load(model_path)
    st.success("Model loaded. Features: %d" % len(bundle["features"]))
