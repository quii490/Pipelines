#!/usr/bin/env python3
"""Plot a two-sample RNA-seq count/expression correlation scatter."""

from __future__ import annotations

import argparse
from pathlib import Path

np = None
pd = None
plt = None
sns = None
pearsonr = None
spearmanr = None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Draw a two-sample count/expression scatter with Pearson/Spearman correlation."
    )
    p.add_argument("--matrix", required=True, help="CSV matrix; first column feature ID and sample columns.")
    p.add_argument("--sample-x", required=True, help="Sample column for x-axis.")
    p.add_argument("--sample-y", required=True, help="Sample column for y-axis.")
    p.add_argument("--out-prefix", required=True, help="Output prefix.")
    p.add_argument("--transform", choices=["log2cpm", "log2count", "raw"], default="log2cpm")
    p.add_argument("--feature-col", default=None, help="Feature ID column. Default: first column.")
    p.add_argument("--annotation-column", default=None, help="Optional column used as hue.")
    p.add_argument("--highlight-file", default=None, help="Optional CSV with feature IDs to highlight/label.")
    p.add_argument("--highlight-col", default=None, help="Feature ID column in highlight file. Default: first column.")
    p.add_argument("--label-top-n", type=int, default=20, help="Label top N features by absolute difference.")
    p.add_argument("--alpha", type=float, default=0.45)
    p.add_argument("--point-size", type=float, default=12)
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def read_matrix(path: str, feature_col: str | None) -> tuple[pd.DataFrame, str]:
    df = pd.read_csv(path)
    if df.shape[1] < 3:
        raise ValueError("matrix needs feature ID plus at least two sample columns")
    if feature_col is None:
        feature_col = df.columns[0]
    if feature_col not in df.columns:
        raise ValueError(f"feature column not found: {feature_col}")
    return df, feature_col


def transform_values(df: pd.DataFrame, x: str, y: str, method: str) -> tuple[pd.Series, pd.Series, str]:
    xv = pd.to_numeric(df[x], errors="coerce").fillna(0)
    yv = pd.to_numeric(df[y], errors="coerce").fillna(0)
    if method == "raw":
        return xv, yv, "raw count / signal"
    if method == "log2count":
        return np.log2(xv + 1), np.log2(yv + 1), "log2(count + 1)"
    xlib = xv.sum() or 1
    ylib = yv.sum() or 1
    return np.log2(xv / xlib * 1e6 + 1), np.log2(yv / ylib * 1e6 + 1), "log2(CPM + 1)"


def main() -> None:
    global np, pd, plt, sns, pearsonr, spearmanr
    args = parse_args()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt_local
    import numpy as np_local
    import pandas as pd_local
    import seaborn as sns_local
    from scipy.stats import pearsonr as pearsonr_local, spearmanr as spearmanr_local

    np = np_local
    pd = pd_local
    plt = plt_local
    sns = sns_local
    pearsonr = pearsonr_local
    spearmanr = spearmanr_local

    df, feature_col = read_matrix(args.matrix, args.feature_col)
    for sample in (args.sample_x, args.sample_y):
        if sample not in df.columns:
            raise ValueError(f"sample column not found: {sample}")

    df = df.copy()
    df["x_value"], df["y_value"], axis_label = transform_values(df, args.sample_x, args.sample_y, args.transform)
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["x_value", "y_value"])

    pear = pearsonr(df["x_value"], df["y_value"])[0] if len(df) >= 3 else np.nan
    spear = spearmanr(df["x_value"], df["y_value"])[0] if len(df) >= 3 else np.nan

    hue = args.annotation_column if args.annotation_column in df.columns else None
    if hue is not None:
        top_groups = df[hue].fillna("Unknown").astype(str).value_counts().head(12).index
        df[hue] = np.where(df[hue].fillna("Unknown").astype(str).isin(top_groups), df[hue].fillna("Unknown").astype(str), "Other")

    highlight_ids: set[str] = set()
    if args.highlight_file:
        hdf = pd.read_csv(args.highlight_file)
        hcol = args.highlight_col or hdf.columns[0]
        if hcol not in hdf.columns:
            raise ValueError(f"highlight column not found: {hcol}")
        highlight_ids = set(hdf[hcol].astype(str))
    df["highlight"] = df[feature_col].astype(str).isin(highlight_ids)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    sns.set_theme(style="whitegrid", context="notebook")
    fig, ax = plt.subplots(figsize=(6.8, 6.2))
    if hue:
        sns.scatterplot(
            data=df,
            x="x_value",
            y="y_value",
            hue=hue,
            s=args.point_size,
            alpha=args.alpha,
            linewidth=0,
            ax=ax,
        )
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0, frameon=False, fontsize=8)
    else:
        ax.scatter(df["x_value"], df["y_value"], s=args.point_size, alpha=args.alpha, color="#4C78A8", linewidths=0)

    if df["highlight"].any():
        hd = df[df["highlight"]]
        ax.scatter(hd["x_value"], hd["y_value"], s=args.point_size * 2.2, facecolors="none", edgecolors="#D55E00", linewidths=0.9)

    label_df = df.assign(abs_diff=(df["y_value"] - df["x_value"]).abs()).sort_values("abs_diff", ascending=False).head(args.label_top_n)
    for _, row in label_df.iterrows():
        ax.text(row["x_value"], row["y_value"], str(row[feature_col]), fontsize=7, alpha=0.85)

    lo = min(df["x_value"].min(), df["y_value"].min())
    hi = max(df["x_value"].max(), df["y_value"].max())
    ax.plot([lo, hi], [lo, hi], linestyle="--", color="black", linewidth=0.8)
    ax.set_xlabel(f"{args.sample_x} {axis_label}")
    ax.set_ylabel(f"{args.sample_y} {axis_label}")
    ax.set_title(f"{args.sample_x} vs {args.sample_y}\nPearson r={pear:.3f}; Spearman rho={spear:.3f}; n={len(df):,}")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.scatter.pdf")
    fig.savefig(f"{out_prefix}.scatter.png", dpi=args.dpi)

    df[[feature_col, "x_value", "y_value", "highlight"]].to_csv(f"{out_prefix}.scatter_values.csv", index=False)


if __name__ == "__main__":
    main()
