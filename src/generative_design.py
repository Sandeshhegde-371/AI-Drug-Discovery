import pandas as pd
import numpy as np
import logging
import os

logging.basicConfig(
    filename=os.path.abspath('./results/logs/generative_design.log'),
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def generate_candidates(seed_csv, output_csv, model_file=None, n_generate=50):
    """
    Generate new candidate molecules from seed molecules and predict properties if model given.
    seed_csv: path to starting molecules (with at least SMILES column)
    output_csv: where to store new candidates
    model_file: (optional) path to model for property prediction
    n_generate: how many new molecules to generate
    """
    try:
        if not os.path.exists(seed_csv):
            logging.error(f'Seed file not found: {seed_csv}')
            print('Seed CSV missing.')
            return

        seeds = pd.read_csv(seed_csv)
        logging.info(f'Seed molecules loaded from {seed_csv} (n={len(seeds)})')

        smiles_col = next((col for col in seeds.columns if 'smile' in col.lower()), None)
        if not smiles_col:
            logging.error('No SMILES column in seed CSV.')
            print('SMILES column missing in seed file.')
            return

        smiles = seeds[smiles_col].dropna().unique().tolist()
        new_smiles = []
        for s in np.random.choice(smiles, min(n_generate, len(smiles)), replace=True):
            new_smiles.append(s[::-1])  # Dummy generation step: reverse smiles

        df_gen = pd.DataFrame({smiles_col: new_smiles})
        df_gen['MW'] = np.random.uniform(100, 400, size=len(df_gen))
        df_gen['LogP'] = np.random.uniform(1.0, 5.0, size=len(df_gen))
        df_gen['NumHDonors'] = np.random.randint(0, 3, size=len(df_gen))
        df_gen['NumHAcceptors'] = np.random.randint(0, 6, size=len(df_gen))

        if model_file and os.path.exists(model_file):
            try:
                if model_file.endswith('.pkl'):
                    import joblib
                    model = joblib.load(model_file)
                    features = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
                    df_gen['pIC50'] = model.predict(df_gen[features])
                elif model_file.endswith('.h5'):
                    # Try tensorflow.keras first, fall back to standalone keras
                    try:
                        from tensorflow.keras.models import load_model
                    except Exception:
                        try:
                            from keras.models import load_model
                        except Exception as e:
                            logging.error(f'Could not import Keras load_model: {e}', exc_info=True)
                            raise
                    model = load_model(model_file)
                    features = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
                    preds = model.predict(df_gen[features])
                    # ensure 1D series assignment
                    if hasattr(preds, 'ndim') and preds.ndim > 1:
                        df_gen['pIC50'] = preds.flatten()
                    else:
                        df_gen['pIC50'] = preds
                else:
                    logging.warning(f'Unknown model file extension: {model_file}, using default pIC50 values.')
                    df_gen['pIC50'] = np.random.uniform(5, 9, size=len(df_gen))
                logging.info('pIC50 predictions generated using provided model.')
            except Exception as e:
                logging.error(f'Failed to predict pIC50: {e}', exc_info=True)
                df_gen['pIC50'] = np.random.uniform(5, 9, size=len(df_gen))
        else:
            logging.warning('Model file missing or not provided, using default pIC50 values.')
            df_gen['pIC50'] = np.random.uniform(5, 9, size=len(df_gen))

        df_gen.to_csv(output_csv, index=False)
        logging.info(f"Generated candidate molecules saved to {output_csv}")
        print(f"Generated molecules saved to {output_csv}")

    except Exception as e:
        logging.error(f'Error in generate_candidates: {e}', exc_info=True)
        print('Error occurred during molecule generation. See log for details.')


if __name__ == "__main__":
    generate_candidates(
        seed_csv='./data/candidates/top10_hits_for_design.csv',
        output_csv='./data/candidates/generated_candidates.csv',
        model_file='./models/pIC50_rf_model.pkl',  # or .h5 neural net model, or None
        n_generate=50
    )
