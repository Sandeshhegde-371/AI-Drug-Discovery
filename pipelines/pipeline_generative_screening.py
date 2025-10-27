import argparse
import logging
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src import generative_design

logging.basicConfig(
    filename='./results/logs/pipeline_generative_screening.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Input CSV with seed compounds')
    parser.add_argument('--num_candidates', type=int, default=50, help='Compounds to generate')
    args = parser.parse_args()
    try:
        generative_design.generate_candidates(
            seed_csv=args.input,
            output_csv='./results/generated_candidates.csv',
            model_file='./models/multitask_admet_model.h5',
            n_generate=args.num_candidates
        )
        logging.info("Generative screening completed successfully.")
    except Exception as e:
        logging.error(f"Pipeline failed: {e}", exc_info=True)
        print("Generative screening error (check log).")

if __name__ == "__main__":
    main()
