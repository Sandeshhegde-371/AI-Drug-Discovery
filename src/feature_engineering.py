import pandas as pd
import logging
import os

logging.basicConfig(
    filename='./results/logs/feature_engineering.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def add_custom_features(input_csv, output_csv):
    """
    Example feature creation: Add a 'Lipinski_Pass' column based on rule-of-five criteria.
    """
    try:
        if not os.path.exists(input_csv):
            logging.error(f'Input data file not found: {input_csv}')
            print('Input data file missing.')
            return None

        df = pd.read_csv(input_csv)
        logging.info(f'Data loaded from {input_csv}, shape={df.shape}')

        # Feature: Lipinski's rule-of-five filter
        df['Lipinski_Pass'] = (
            (df['MW'] <= 500) &
            (df['LogP'] <= 5) &
            (df['NumHDonors'] <= 5) &
            (df['NumHAcceptors'] <= 10)
        )
        logging.info(f'Lipinski_Pass feature added for {df.shape[0]} molecules.')

        # Additional feature engineering can be added here
        df.to_csv(output_csv, index=False)
        logging.info(f'Feature engineered data saved to {output_csv}.')
        print(f'Features added and data written to {output_csv}.')
        return df

    except Exception as e:
        logging.error(f'Feature engineering error: {e}')
        print('Error occurred during feature engineering. See logs for details.')
        return None

if __name__ == "__main__":
    add_custom_features(
        './data/processed/normalized_data_features.csv',
        './data/processed/normalized_data_with_lipinski.csv'
    )
