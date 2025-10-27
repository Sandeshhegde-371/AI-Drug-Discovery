import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from tensorflow import keras
import logging
import os

# Set up logging
logging.basicConfig(
    filename='./results/logs/generative.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def randomize_smiles(smiles, n_variants=10):
    variants = set()
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            logging.warning(f'Invalid SMILES for randomization: {smiles}')
            return []
        for _ in range(n_variants):
            randomized = Chem.MolToSmiles(mol, doRandom=True)
            variants.add(randomized)
    except Exception as e:
        logging.error(f'Error generating random SMILES for {smiles}: {e}')
    return list(variants)

def generate_candidates(seed_csv, output_csv, model_file):
    try:
        # Verify input files exist
        if not os.path.exists(seed_csv):
            logging.error(f'Seed file not found: {seed_csv}')
            print('Seed CSV not found. Check provided path.')
            return
        if not os.path.exists(model_file):
            logging.error(f'Model file not found: {model_file}')
            print('Model file not found. Check provided path.')
            return

        seed_df = pd.read_csv(seed_csv)
        logging.info(f'Seed molecules loaded: {seed_csv}')

        unique_smiles = seed_df['canonical_smiles'].dropna().unique()
        new_smiles = []
        for sm in unique_smiles:
            variants = randomize_smiles(sm, n_variants=10)
            new_smiles.extend(variants)
        new_smiles = list(set(new_smiles))
        logging.info(f'Generated {len(new_smiles)} unique molecule SMILES.')

        rows = []
        for sm in new_smiles:
            try:
                mol = Chem.MolFromSmiles(sm)
                if mol:
                    mw = Descriptors.MolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    donors = Descriptors.NumHDonors(mol)
                    acceptors = Descriptors.NumHAcceptors(mol)
                    rows.append({
                        'canonical_smiles': sm, 
                        'MW': mw, 
                        'LogP': logp,
                        'NumHDonors': donors, 
                        'NumHAcceptors': acceptors
                    })
                else:
                    logging.warning(f'Could not featurize SMILES: {sm}')
            except Exception as fe:
                logging.error(f'Feature extraction failed for {sm}: {fe}')
        
        df_features = pd.DataFrame(rows)
        df_features.to_csv(output_csv, index=False)
        logging.info(f'Saved generated candidates to {output_csv}')

        X = df_features[['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']].fillna(0).values
        model = keras.models.load_model(model_file, compile=False)
        logging.info(f'Model loaded from {model_file} for prediction.')

        preds = model.predict(X)
        output_cols = [
            'cyp_3a4_reg_value','hum_mic_cl_reg_value','cyp_2c9_reg_value','logd_reg_value',
            'ames_cls_value','cyp_2d6_reg_value','ppb_reg_value','bbb_cls_value',
            'rat_mic_cl_reg_value','mou_mic_cl_reg_value','water_sol_reg_value'
        ]
        df_preds = pd.DataFrame(preds, columns=output_cols)
        df_final = pd.concat([df_features, df_preds], axis=1)
        pred_path = './data/candidates/generated_candidates_with_preds.csv'
        df_final.to_csv(pred_path, index=False)
        logging.info(f'Predictions saved to {pred_path}')
        print(f'Augmented molecules predicted and saved in {pred_path}')

    except Exception as e:
        logging.error(f'Error in candidate generation: {e}')
        print('An error occurred in candidate generation. See log for details.')

if __name__ == "__main__":
    # Update the following paths for your project structure
    generate_candidates(
        './data/candidates/top10_hits_for_design.csv',
        './data/candidates/generated_candidates.csv',
        './models/multitask_admet_model.h5'
    )
