import argparse
import logging
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src import predict_for_new_compounds

logging.basicConfig(
    filename='./results/logs/pipeline_predict_properties.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Input compound CSV')
    args = parser.parse_args()
    try:
        predict_for_new_compounds.predict_activity(
            input_csv=args.input,
            model_file='./models/multitask_admet_model.h5',
            output_csv='./results/property_predictions.csv'
        )
        logging.info("Property prediction pipeline completed successfully.")
    except Exception as e:
        logging.error(f"Pipeline failed: {e}", exc_info=True)
        print("Prediction error (check log).")

if __name__ == "__main__":
    main()
