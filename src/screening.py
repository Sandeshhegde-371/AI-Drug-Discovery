import pandas as pd
import numpy as np
import logging
import os

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, DataStructs

# For ML predictions
import joblib
from tensorflow.keras.models import load_model

logging.basicConfig(
    filename='./results/logs/virtual_screening.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)


def filter_compounds(df, filters):
    """Apply fast property-based filters to compound dataframe."""
    mask = np.ones(len(df), dtype=bool)
    for prop, (low, high) in filters.items():
        mask = mask & (df[prop] >= low) & (df[prop] <= high)
    return df[mask]


def compute_descriptors(smiles):
    """Compute descriptors for a SMILES (can expand as needed)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    desc = {}
    desc['MW'] = Descriptors.MolWt(mol)
    desc['LogP'] = Descriptors.MolLogP(mol)
    desc['NumHDonors'] = Descriptors.NumHDonors(mol)
    desc['NumHAcceptors'] = Descriptors.NumHAcceptors(mol)
    # Add additional descriptors as needed
    return desc


def tanimoto_similarity(smiles, ref_fp):
    """Compute Tanimoto similarity between a compound SMILES and a reference fingerprint."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0.0
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    return DataStructs.TanimotoSimilarity(fp, ref_fp)


def virtual_screen(
    input_csv,
    output_csv,
    model_path,
    reference_smiles=None,
    filter_params=None,
    top_n=10,
    priority_metric="Binding_Pred_1"
):
    """
    Virtual screening pipeline:
    - Applies fast filtering
    - Predicts binding/score with ML
    - Calculates similarity to reference
    - Ranks and saves top compounds
    """
    try:
        if not os.path.exists(input_csv):
            logging.error(f'Input file not found: {input_csv}')
            print('Input file missing.')
            return None

        df = pd.read_csv(input_csv)
        logging.info(f'Loaded {len(df)} compounds for screening.')

        # Compute descriptors (if not present)
        if 'MW' not in df.columns or 'LogP' not in df.columns:
            desc_list = []
            for smi in df['SMILES']:
                desc = compute_descriptors(smi)
                if desc is None:
                    desc = {'MW':np.nan, 'LogP':np.nan, 'NumHDonors':np.nan, 'NumHAcceptors':np.nan}
                desc_list.append(desc)
            desc_df = pd.DataFrame(desc_list)
            df = pd.concat([df, desc_df], axis=1)

        # Fast filtering
        filters = filter_params or {
            'MW': (200, 500),
            'LogP': (-1, 5),
            'NumHDonors': (0, 5),
            'NumHAcceptors': (0, 10)
        }
        df = filter_compounds(df, filters)
        logging.info(f'Filtered down to {len(df)} candidates after property filters.')

        # Target binding prediction
        model = None
        if model_path.endswith('.pkl'):
            model = joblib.load(model_path)
            ml_features = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
        elif model_path.endswith('.h5'):
            model = load_model(model_path, compile=False)
            ml_features = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors'] # Adjust if needed
        else:
            raise ValueError("Unsupported model file format.")

        X = df[ml_features]
        preds = model.predict(X)
        if len(preds.shape) == 1 or preds.shape[1] == 1:
            df['Binding_Pred'] = preds.flatten()
        else:
            for i in range(preds.shape[1]):
                df[f'Binding_Pred_{i+1}'] = preds[:, i]


        # Similarity scoring (if reference given)
        if reference_smiles is not None:
            ref_mol = Chem.MolFromSmiles(reference_smiles)
            ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=2048)
            df['Sim_to_Ref'] = [
                tanimoto_similarity(sm, ref_fp) if Chem.MolFromSmiles(sm) else 0.0 
                for sm in df['SMILES']
            ]
        else:
            df['Sim_to_Ref'] = np.nan

        # Rank: prioritize by binding prediction, then similarity
        if priority_metric not in df.columns:
            logging.error(f"{priority_metric} not found in dataframe columns.")
            print(f"Column {priority_metric} not found.")
            return None

        df = df.sort_values([priority_metric, 'Sim_to_Ref'], ascending=False)
        top_hits = df.head(top_n)
        top_hits.to_csv(output_csv, index=False)

        print(f"Virtual screening complete. Top {top_n} hits saved to {output_csv}.")
        return top_hits

    except Exception as e:
        logging.error(f'Virtual screening error: {e}')
        print('Virtual screening error (check log).')
        return None


if __name__ == "__main__":
    # Example Usage
    virtual_screen(
        input_csv='./data/candidates/candidate_data.csv',
        output_csv='./results/virtual_screen_top10.csv',
        model_path='./models/activity_rf_model.pkl',    # or .h5 model
        reference_smiles='CCOCCO',                      # Replace with a real active SMILES if available
        filter_params=None,                             # Use defaults or pass a dict
        top_n=10
    )
