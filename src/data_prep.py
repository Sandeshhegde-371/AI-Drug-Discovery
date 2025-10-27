import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import logging
import os

# Set up logging
logging.basicConfig(
    filename='./results/logs/data_prep.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def get_descriptors(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            donors = Descriptors.NumHDonors(mol)
            acceptors = Descriptors.NumHAcceptors(mol)
            return mw, logp, donors, acceptors
        else:
            logging.warning(f'Invalid SMILES string: {smiles}')
            return None, None, None, None
    except Exception as e:
        logging.error(f'Failed to process SMILES {smiles}: {e}')
        return None, None, None, None

def prepare_main_data():
    try:
        input_path = './data/processed/normalized_data.csv'
        if not os.path.exists(input_path):
            logging.error(f'Input file not found: {input_path}')
            return None

        norm_df = pd.read_csv(input_path)
        logging.info('Read normalized_data.csv successfully.')

        # Drop duplicates and NAs
        norm_df = norm_df.drop_duplicates(subset=['molecule_chembl_id'])
        norm_df = norm_df.dropna(subset=['canonical_smiles', 'pIC50'])
        if norm_df.empty:
            logging.warning('No data left after dropping NAs and duplicates.')

        norm_df['activityclass'] = norm_df['activity_class'].str.lower()

        # Feature engineering
        norm_df[['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']] = (
            norm_df['canonical_smiles'].apply(get_descriptors).apply(pd.Series)
        )

        output_path = './data/processed/normalized_data_features.csv'
        norm_df.to_csv(output_path, index=False)
        logging.info(f'Successfully saved feature engineered data to {output_path}.')
        print(f"Feature-engineered normalized data saved to {output_path}. Check logs for details.")
        return norm_df

    except Exception as e:
        logging.error(f"Error in prepare_main_data: {e}")
        print("An error occurred during data preparation. See logs for details.")
        return None

if __name__ == "__main__":
    prepare_main_data()
