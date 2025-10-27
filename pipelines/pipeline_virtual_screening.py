import argparse
import logging
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src import screening

logging.basicConfig(
    filename='./results/logs/pipeline_virtual_screening.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Library CSV')
    parser.add_argument('--num_candidates', type=int, default=100, help='Top N compounds')
    args = parser.parse_args()
    try:
        screening.virtual_screen(
            input_csv=args.input,
            model_file='./models/multitask_admet_model.h5',
            output_csv='./results/virtual_screen_results.csv',
            top_n=args.num_candidates
        )
        logging.info("Virtual screening pipeline completed successfully.")
    except Exception as e:
        logging.error(f"Pipeline failed: {e}", exc_info=True)
        print("Virtual screening error (check log).")

if __name__ == "__main__":
    main()
