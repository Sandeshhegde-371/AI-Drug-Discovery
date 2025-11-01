import pandas as pd
import numpy as np
import logging
import os

# Configure logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def normalize(series, higher_better=True):
    """Normalize a pandas series to range 0–1."""
    if series.nunique() == 1:
        return np.ones(len(series)) * 0.5
    normalized = (series - series.min()) / (series.max() - series.min())
    return normalized if higher_better else 1 - normalized


def compute_drug_likeness(df):
    """Compute a simple drug-likeness score based on Lipinski’s rules."""
    df["Lipinski_Score"] = 0
    df.loc[df["MW"] <= 500, "Lipinski_Score"] += 1
    df.loc[df["LogP"] <= 5, "Lipinski_Score"] += 1
    df.loc[df["NumHDonors"] <= 5, "Lipinski_Score"] += 1
    df.loc[df["NumHAcceptors"] <= 10, "Lipinski_Score"] += 1
    return df["Lipinski_Score"] / 4.0  # normalize to 0–1


def prioritize_candidates(pred_csv, output_csv, priority_metric="Predicted Activity", top_n=50):
    """Rank and prioritize compounds based on multiple scoring criteria."""

    # Load dataset
    if not os.path.exists(pred_csv):
        raise FileNotFoundError(f"Input file not found: {pred_csv}")

    df = pd.read_csv(pred_csv)
    logger.info(f"Loaded {len(df)} candidates from {pred_csv}")

    # --- Compute additional metrics ---
    # Safety Profile: heuristic (low MW and low LogP often safer)
    df["Safety_Profile"] = normalize((1 / (df["MW"] + 1e-6)) + (1 / (df["LogP"].abs() + 1e-6)))

    # Drug-Likeness Score: rule-based
    df["Drug_Likeness_Score"] = compute_drug_likeness(df)

    # Ease of Synthesis: heuristic (smaller MW, fewer H donors/acceptors)
    df["Ease_of_Synthesis"] = normalize(
        (1 / (df["MW"] + 1e-6)) + (1 / (df["NumHDonors"] + 1)) + (1 / (df["NumHAcceptors"] + 1))
    )

    # Predicted Activity: normalize
    higher_is_better = True
    if df["predicted_activity"].mean() < 0:
        # if negative (like docking scores), lower is better
        higher_is_better = False
    df["Predicted_Activity_Score"] = normalize(df["predicted_activity"], higher_better=higher_is_better)

    # Balanced Multi-Criteria score
    weights = {
        "Predicted_Activity_Score": 0.4,
        "Safety_Profile": 0.2,
        "Drug_Likeness_Score": 0.2,
        "Ease_of_Synthesis": 0.2
    }
    df["Balanced_Score"] = sum(df[k] * w for k, w in weights.items())

    # --- Select sorting metric ---
    metric_map = {
        "Predicted Activity": "Predicted_Activity_Score",
        "Safety Profile": "Safety_Profile",
        "Drug-Likeness Score": "Drug_Likeness_Score",
        "Ease of Synthesis": "Ease_of_Synthesis",
        "Balanced Multi-Criteria": "Balanced_Score"
    }

    metric_col = metric_map.get(priority_metric, "Balanced_Score")

    # --- Rank and save ---
    df = df.sort_values(by=metric_col, ascending=False).reset_index(drop=True)
    top_df = df.head(top_n)

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    top_df.to_csv(output_csv, index=False)

    logger.info(f"Top {top_n} candidates saved to {output_csv}")
    print(f"\n✅ Prioritization complete: top {top_n} saved to {output_csv}")
    print(f"Ranking metric used: {priority_metric}")

    return top_df
