import logging
import sys
from pathlib import Path

# Ensure project root (parent of this pipelines directory) is on sys.path so 'src' can be imported.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import screening, selecting_best_class

logging.basicConfig(
    filename='./results/logs/pipeline_screening_and_selection.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def run_pipeline():
    logging.info("Starting candidate screening pipeline...")
    screening.prioritize_candidates(
        pred_csv='./data/candidates/predicted_candidates.csv',
        output_csv='./data/candidates/top10_hits_for_design.csv',
        priority_metric='pIC50',
        top_n=10
    )
    logging.info("Candidate screening completed.")

    logging.info("Selecting best classification model...")
    selecting_best_class.select_best_model(
        metrics_csv='./results/summaries/classification_metrics.csv',
        output_csv='./results/summaries/best_classification_model.csv',
        metric_name='accuracy'
    )
    logging.info("Best model selection completed.")

if __name__ == "__main__":
    run_pipeline()
