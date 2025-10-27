import pandas as pd
import logging
import os

# Set up logging
logging.basicConfig(
    filename='./results/logs/merging.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def merge_datasets(input_files, output_path, key='canonical_smiles'):
    """
    Merge multiple input CSV files on a common key, reports merging issues, saves merged result.
    """
    try:
        dfs = []
        for f in input_files:
            if not os.path.exists(f):
                logging.error(f'Dataset file {f} not found.')
                print(f'File missing: {f}')
                continue
            df = pd.read_csv(f)
            # Standardize key column name
            if 'Smiles_unify' in df.columns:
                df.rename(columns={'Smiles_unify': 'canonical_smiles'}, inplace=True)
            dfs.append(df)
            logging.info(f'Loaded dataset: {f} shape={df.shape}, columns={list(df.columns)}')
        if not dfs:
            logging.error('No datasets to merge; exiting.')
            print('No datasets available for merging.')
            return None

        merged = dfs[0]
        for idx, df in enumerate(dfs[1:], 2):
            merged = pd.merge(merged, df, on=key, how='outer')
            logging.info(f'Merged {idx} datasets, current shape: {merged.shape}')
        merged.to_csv(output_path, index=False)
        logging.info(f'Successfully saved merged dataset to {output_path}')
        print(f'Merged dataset saved to {output_path}')
        return merged

    except Exception as e:
        logging.error(f'Error during dataset merging: {e}')
        print('An error occurred while merging datasets. See log for details.')
        return None

if __name__ == "__main__":
    # Example usage; replace with actual dataset paths
    input_csvs = [
        './data/processed/normalized_data_features.csv',
        './data/processed/logd_reg_final_data.csv',
        './data/processed/water_sol_reg_final_data.csv'
    ]
    merge_datasets(input_csvs, './data/processed/merged_multitask_data.csv', key='canonical_smiles')
