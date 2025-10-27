import pandas as pd
import logging
import joblib
import os

logging.basicConfig(
    filename='./results/logs/prediction.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def predict_activity(candidate_csv, model_path, output_csv, model_type='rf'):
    """
    Predict pIC50 values for candidate compounds using a trained model and save results.
    """
    try:
        if not os.path.exists(candidate_csv):
            logging.error(f'Candidate file not found: {candidate_csv}')
            print('Candidate file missing.')
            return
        df = pd.read_csv(candidate_csv)
        logging.info(f'Loaded candidates: {candidate_csv}, shape: {df.shape}')

        if not os.path.exists(model_path):
            logging.error(f'Model file not found: {model_path}')
            print('Model file missing.')
            return

        # Load trained model
        model = joblib.load(model_path)
        logging.info(f'Model loaded from: {model_path}')

        # Assume features used for prediction are the following
        features = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
        if not all(f in df.columns for f in features):
            logging.error(f'Missing required features: {features}')
            print('Required features missing in candidate file.')
            return

        # Predict pIC50 values (regression)
        df['pIC50'] = model.predict(df[features])
        logging.info('Predictions made for pIC50.')

        # (Optional) If classification also needed, e.g. activity class
        if hasattr(model, 'predict_proba'):
            df['predicted_class'] = model.predict(df[features])

        df.to_csv(output_csv, index=False)
        logging.info(f"Prediction results saved to {output_csv}")
        print(f"Prediction results saved to {output_csv}")

    except Exception as e:
        logging.error(f'Error during prediction: {e}')
        print('Prediction failed. See log for details.')

if __name__ == "__main__":
    predict_activity(
        candidate_csv='./data/candidates/generated_candidates.csv',
        model_path='./models/activity_rf_model.pkl',
        output_csv='./results/property_predictions.csv',
        model_type='rf'
    )
