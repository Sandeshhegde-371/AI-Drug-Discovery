import logging
import sys
from pathlib import Path

# Ensure project root (parent of this pipelines directory) is on sys.path so 'src' can be imported.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import model_training, predict

logging.basicConfig(
    filename='./results/logs/pipeline_train_predict.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def run_pipeline():
    logging.info("Starting model training pipeline...")
    model_training.train_models()
    logging.info("Model training completed.")

    logging.info("Starting prediction on generated candidates...")
    predict.predict_activity(
        candidate_csv='./data/candidates/generated_candidates.csv',
        model_path='./models/activity_rf_model.pkl',
        output_csv='./data/candidates/predicted_candidates.csv',
        model_type='rf'
    )
    logging.info("Prediction completed.")

if __name__ == "__main__":
    run_pipeline()
