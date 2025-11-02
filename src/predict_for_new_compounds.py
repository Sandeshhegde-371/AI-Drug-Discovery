import pandas as pd
import os
import joblib
import logging

from tensorflow.keras.models import load_model

logging.basicConfig(
    filename='./results/logs/prediction.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def predict_activity(candidate_csv, models_dir, output_csv):
    """
    Predicts activity, pIC50, and multitask ADMET properties for candidate compounds.
    Uses only available models.
    """
    try:
        if not os.path.exists(candidate_csv):
            logging.error(f'Candidate file not found: {candidate_csv}')
            print('Candidate file missing.')
            return
        df = pd.read_csv(candidate_csv)
        logging.info(f'Loaded candidates: {candidate_csv}, shape: {df.shape}')
        
        # 1. Activity Prediction (Random Forest)
        act_rf_path = os.path.join(models_dir, 'activity_rf_model.pkl')
        if os.path.exists(act_rf_path):
            act_rf_model = joblib.load(act_rf_path)
            act_features = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
            if all(f in df.columns for f in act_features):
                df['Activity(0-active,1-inactive)'] = act_rf_model.predict(df[act_features])
                logging.info('Random forest activity predictions done.')

        # Activity prediction using neural net
        act_nn_path = os.path.join(models_dir, 'activity_nn_model.h5')
        if os.path.exists(act_nn_path):
            act_nn_model = load_model(act_nn_path)
            nn_features = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
            if all(f in df.columns for f in nn_features):
                y_pred = act_nn_model.predict(df[nn_features])
                print("NN prediction shape:", y_pred.shape)
                if len(y_pred.shape) == 1 or y_pred.shape[1] == 1:
                    df['Activity_NN_Pred'] = y_pred.flatten()
                else:
                    for i in range(y_pred.shape[1]):
                        df[f'Activity_NN_Probability_{i+1}'] = y_pred[:, i]
                logging.info('Neural network activity predictions done.')

        # 3. pIC50 Prediction
        pic50_path = os.path.join(models_dir, 'pIC50_rf_model.pkl')
        if os.path.exists(pic50_path):
            pic50_model = joblib.load(pic50_path)
            pic50_features = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
            if all(f in df.columns for f in pic50_features):
                df['pIC50_RF_Pred'] = pic50_model.predict(df[pic50_features])
                logging.info('pIC50 predictions done.')

        # 4. ADMET (Multitask)
        admet_path = os.path.join(models_dir, 'multitask_admet_model.h5')
        if os.path.exists(admet_path):
            admet_model = load_model(admet_path, compile=False)
            admet_features = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']  
            if all(f in df.columns for f in admet_features):
                admet_preds = admet_model.predict(df[admet_features])
                # If multitask, admet_preds might be a 2D array, one column per endpoint
                for idx, col in enumerate(['Solubility', 'Permeability']):  
                    df[col] = admet_preds[:, idx]
                logging.info('ADMET predictions done.')
        df.to_csv(output_csv, index=False)
        logging.info(f"All available property predictions saved to {output_csv}")
        print(f"All available property predictions saved to {output_csv}")

    except Exception as e:
        logging.error(f'Error during prediction: {e}')
        preds = admet_model.predict(df[admet_features])
        print(preds.shape)
        print('Prediction failed. See log for details.')

if __name__ == "__main__":
    predict_activity(
        candidate_csv='./data/candidates/generated_candidates.csv',
        models_dir='./models/',
        output_csv='./results/property_predictions.csv'
    )
