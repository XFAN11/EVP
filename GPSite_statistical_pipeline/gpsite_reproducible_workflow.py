"""
End-to-end example GPSite statistical workflow.

This script starts from a GPSite output table containing mutant and WT binding
scores, recalculates delta-binding values, performs quantile normalization to
obtain Z-scores, generates a QQ plot, and applies the nominal filtering rule
derived from a two-sided P-value threshold of 0.05.

Recommended environment
-----------------------
This workflow follows the manuscript analysis logic. Environment-dependent
library versions, especially ``scikit-learn``, can slightly change the
quantile-normalized summary values.

Example usage
-------------
python gpsite_reproducible_workflow.py

Optional arguments
------------------
python gpsite_reproducible_workflow.py --input GPSite_data_1.csv --prefix review
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import QuantileTransformer


BINDING_COLS = [
    "DNA_binding",
    "RNA_binding",
    "Peptide_binding",
    "Protein_binding",
    "ATP_binding",
    "ZN_binding",
    "CA_binding",
    "MG_binding",
    "MN_binding",
]

P_THRESHOLD = 0.05
Z_THRESHOLD = stats.norm.ppf(1 - P_THRESHOLD / 2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproducible GPSite delta-binding normalization workflow."
    )
    parser.add_argument(
        "--input",
        default="GPSite_data_1.csv",
        help="GPSite output table with mutant and WT binding columns.",
    )
    parser.add_argument(
        "--prefix",
        default="gpsite_reproducible",
        help="Prefix for exported figure and result tables.",
    )
    return parser.parse_args()


def load_and_validate_input(input_path: Path) -> pd.DataFrame:
    df = pd.read_csv(input_path)

    required_cols = {"Mut_Sample", "Transcript", "mut_info", "AA"}
    required_cols.update(BINDING_COLS)
    required_cols.update(f"{col}_WT" for col in BINDING_COLS)

    missing = sorted(required_cols.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    return df


def build_long_table(df: pd.DataFrame) -> pd.DataFrame:
    long_frames = []

    # Keep the same column-wise concatenation order as cut_off_6.py so the
    # quantile normalization statistics reproduce the manuscript values.
    for binding_col in BINDING_COLS:
        mut_score = df[binding_col]
        wt_score = df[f"{binding_col}_WT"]
        delta = mut_score - wt_score

        both_zero = (mut_score == 0) & (wt_score == 0)
        zero_delta = delta == 0
        keep = (~both_zero) & (~zero_delta)

        long_frames.append(
            pd.DataFrame(
                {
                    "Mut_Sample": df["Mut_Sample"],
                    "Transcript": df["Transcript"],
                    "mut_info": df["mut_info"],
                    "AA": df["AA"],
                    "binding_type": binding_col.replace("_binding", ""),
                    "mut_score": mut_score,
                    "wt_score": wt_score,
                    "delta_binding": delta,
                    "both_zero_filtered": both_zero,
                    "zero_delta_filtered": zero_delta,
                    "kept_for_statistics": keep,
                }
            )
        )

    return pd.concat(long_frames, ignore_index=True)


def quantile_normalize(delta_values: np.ndarray) -> np.ndarray:
    n_quantiles = min(10000, len(delta_values))
    transformer = QuantileTransformer(
        output_distribution="normal",
        n_quantiles=n_quantiles,
        random_state=0,
    )
    return transformer.fit_transform(delta_values.reshape(-1, 1)).ravel()


def derive_delta_threshold(delta_values: np.ndarray, z_scores: np.ndarray) -> float:
    # Reproduce the manuscript logic by mapping the nominal Z cutoff back to the
    # observed delta-binding distribution after quantile normalization.
    threshold_candidates = np.abs(delta_values[np.abs(z_scores) >= Z_THRESHOLD])
    if len(threshold_candidates) == 0:
        raise ValueError("Unable to derive a delta-binding threshold from Z-scores.")
    return float(np.min(threshold_candidates))


def summarize_results(filtered_df: pd.DataFrame, delta_threshold: float) -> pd.Series:
    delta_values = filtered_df["delta_binding"].to_numpy()
    z_scores = filtered_df["z_score"].to_numpy()
    p_values = filtered_df["p_value"].to_numpy()

    summary = pd.Series(
        {
            "n_filtered_variants": int(len(filtered_df)),
            "nominal_p_threshold": P_THRESHOLD,
            "nominal_z_threshold": Z_THRESHOLD,
            "derived_delta_threshold": delta_threshold,
            "n_abs_delta_gt_threshold": int(np.sum(np.abs(delta_values) >= delta_threshold)),
            "n_p_lt_0_05": int(np.sum(p_values < P_THRESHOLD)),
            "n_pass_both_filters": int(np.sum(filtered_df["pass_nominal_filter"])),
            "delta_mean": float(np.mean(delta_values)),
            "delta_median": float(np.median(delta_values)),
            "delta_min": float(np.min(delta_values)),
            "delta_max": float(np.max(delta_values)),
            "z_mean": float(np.mean(z_scores)),
            "z_std": float(np.std(z_scores)),
            "z_skewness": float(stats.skew(z_scores)),
            "z_kurtosis": float(stats.kurtosis(z_scores)),
        }
    )
    return summary


def make_figure(filtered_df: pd.DataFrame, output_path: Path, delta_threshold: float) -> None:
    plt.rcParams["font.family"] = "DejaVu Sans"
    plt.rcParams["axes.unicode_minus"] = False
    sns.set_theme(style="ticks")

    delta_values = filtered_df["delta_binding"].to_numpy()
    z_scores = filtered_df["z_score"].to_numpy()
    skewness = stats.skew(z_scores)
    kurtosis = stats.kurtosis(z_scores)

    fig = plt.figure(figsize=(9, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.28, wspace=0.24)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    sns.histplot(delta_values, bins=80, kde=True, ax=ax_a, color="#74a9cf", edgecolor=None)
    ax_a.axvline(delta_threshold, color="red", linestyle="--", linewidth=1)
    ax_a.axvline(-delta_threshold, color="red", linestyle="--", linewidth=1)
    ax_a.text(
        0.02,
        0.96,
        f"$P < 0.05$ ({delta_threshold:.4f})",
        transform=ax_a.transAxes,
        ha="left",
        va="top",
        fontsize=11,
        style="italic",
    )
    ax_a.set_xlabel("Δbinding")
    ax_a.set_ylabel("Number of mutations")

    sns.histplot(z_scores, bins=80, kde=True, ax=ax_b, color="#74a9cf", edgecolor=None)
    ax_b.axvline(Z_THRESHOLD, color="red", linestyle="--", linewidth=1)
    ax_b.axvline(-Z_THRESHOLD, color="red", linestyle="--", linewidth=1)
    ax_b.text(
        0.05,
        0.96,
        r"$|Z| = 1.96$",
        transform=ax_b.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        style="italic",
    )
    ax_b.set_xlabel("Z-score (quantile normalized)")
    ax_b.set_ylabel("Number of mutations")

    stats.probplot(z_scores, dist="norm", plot=ax_c)
    ax_c.text(
        0.04,
        0.96,
        f"Skewness: {skewness:.4f}\nKurtosis: {kurtosis:.4f}",
        transform=ax_c.transAxes,
        ha="left",
        va="top",
        fontsize=11,
    )
    ax_c.set_xlabel("Theoretical quantiles")
    ax_c.set_ylabel("Observed quantiles")

    for label, ax in zip(["(A)", "(B)", "(C)"], [ax_a, ax_b, ax_c]):
        ax.text(-0.18, 1.08, label, transform=ax.transAxes, fontsize=20, va="top")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    prefix = args.prefix

    df = load_and_validate_input(input_path)
    long_df = build_long_table(df)
    filtered_df = long_df.loc[long_df["kept_for_statistics"]].copy()

    z_scores = quantile_normalize(filtered_df["delta_binding"].to_numpy())
    delta_threshold = derive_delta_threshold(filtered_df["delta_binding"].to_numpy(), z_scores)
    filtered_df["z_score"] = z_scores
    filtered_df["p_value"] = 2 * (1 - stats.norm.cdf(np.abs(z_scores)))
    filtered_df["pass_nominal_filter"] = (
        (filtered_df["delta_binding"].abs() >= delta_threshold)
        & (filtered_df["p_value"] < P_THRESHOLD)
    )

    summary = summarize_results(filtered_df, delta_threshold)

    long_table_path = Path(f"{prefix}_all_long_format.csv")
    filtered_table_path = Path(f"{prefix}_filtered_for_statistics.csv")
    hits_table_path = Path(f"{prefix}_nominal_hits.csv")
    summary_path = Path(f"{prefix}_summary.csv")
    figure_path = Path(f"{prefix}_supplementary_figure8.png")

    long_df.to_csv(long_table_path, index=False)
    filtered_df.to_csv(filtered_table_path, index=False)
    filtered_df.loc[filtered_df["pass_nominal_filter"]].to_csv(hits_table_path, index=False)
    summary.to_frame(name="value").to_csv(summary_path)
    make_figure(filtered_df, figure_path, delta_threshold)

    print("GPSite reproducible workflow completed.")
    print(f"Input file: {input_path}")
    print(f"Total binding records: {len(long_df)}")
    print(f"Records retained for statistics: {len(filtered_df)}")
    print(f"Derived |delta_binding| threshold at P < {P_THRESHOLD}: {delta_threshold:.4f}")
    print(f"Nominal hits (|delta_binding| >= {delta_threshold:.4f}, P < {P_THRESHOLD}): "
          f"{int(filtered_df['pass_nominal_filter'].sum())}")
    print(f"Skewness after quantile normalization: {summary['z_skewness']:.4f}")
    print(f"Kurtosis after quantile normalization: {summary['z_kurtosis']:.4f}")
    print(f"Figure saved to: {figure_path}")
    print(f"Summary saved to: {summary_path}")
    print(f"Filtered hit table saved to: {hits_table_path}")


if __name__ == "__main__":
    main()
