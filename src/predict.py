import pandas as pd
import joblib
from tensorflow import keras
import logging
import os

logging.basicConfig(
    filename='./results/logs/predict.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def predict_activity(candidate_csv, model_path, output_csv, model_type='rf'):
    try:
        if not os.path.exists(candidate_csv):
            logging.error(f'Candidate file not found: {candidate_csv}')
            print('Candidate file not found.')
            return
        if not os.path.exists(model_path):
            logging.error(f'Model file not found: {model_path}')
            print('Model file not found.')
            return

        cand_df = pd.read_csv(candidate_csv)
        feature_cols = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']

        if not all(col in cand_df.columns for col in feature_cols):
            logging.error(f'Missing feature columns in candidate file: {candidate_csv}')
            print('Missing features in candidate CSV.')
            return

        X = cand_df[feature_cols].fillna(0).values
        logging.info(f'Loaded candidate features from {candidate_csv}.')

        if model_type == 'rf':
            model = joblib.load(model_path)
            logging.info(f'RandomForest model loaded from {model_path}.')
            preds = model.predict(X)
            cand_df['predicted_class'] = preds
        else:
            model = keras.models.load_model(model_path, compile=False)
            logging.info(f'Neural network model loaded from {model_path}.')
            preds = model.predict(X)
            cand_df['predicted_class'] = preds.argmax(axis=1)

        cand_df.to_csv(output_csv, index=False)
        logging.info(f'Predictions saved to {output_csv}.')
        print(f'Predictions written to {output_csv}.')

    except Exception as e:
        logging.error(f'Prediction error: {e}')
        print('An error occurred during prediction. See logs for details.')

if __name__ == "__main__":
    predict_activity(
        './data/candidates/generated_candidates.csv',
        './models/activity_rf_model.pkl',
        './data/candidates/predicted_candidates.csv',
        model_type='rf'
    )
