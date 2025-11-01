import argparse
import logging
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src import prioritize_candidates


logging.basicConfig(
    filename='./results/logs/pipeline_prioritize_candidates.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Predicted candidates CSV')
    parser.add_argument('--metric', default='pIC50', help='Property for ranking')
    parser.add_argument('--num_candidates', type=int, default=50, help='Top N compounds')
    args = parser.parse_args()
    try:
        prioritize_candidates.prioritize_candidates(
            pred_csv=args.input,
            output_csv='./data/candidates/prioritized_candidates.csv',
            priority_metric=args.metric,
            top_n=args.num_candidates
        )
        logging.info("Candidate prioritization pipeline completed successfully.")
    except Exception as e:
        logging.error(f"Pipeline failed: {e}", exc_info=True)
        print("Prioritization error (check log).")

if __name__ == "__main__":
    main()
