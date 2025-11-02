import argparse
import logging
import sys
import os
import pandas as pd

sys.stdout.reconfigure(encoding='utf-8')

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src import predict_for_new_compounds

# Configure logging
logging.basicConfig(
    filename='./results/logs/pipeline_predict_properties.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def main():
    parser = argparse.ArgumentParser(description="Predict properties for candidate compounds.")
    parser.add_argument('--input', required=True, help='Input compound CSV')
    parser.add_argument('--num_candidates', type=int, default=None,
                        help='Optionally limit to top N candidates after prediction')
    args = parser.parse_args()

    try:
        output_path = './data/candidates/property_predictions.csv'

        # Run property prediction
        predict_for_new_compounds.predict_activity(
            candidate_csv=args.input,
            models_dir='./models/',
            output_csv=output_path
        )

        logging.info("Property prediction pipeline completed successfully.")

        # If user requested a limit on candidates
        if args.num_candidates:
            df = pd.read_csv(output_path)

            # Try to select a numeric "Score" or "Activity" column if available
            score_cols = [c for c in df.columns if c.lower() in ['score', 'activity', 'prediction']]
            if score_cols:
                sort_col = score_cols[0]
                df = df.sort_values(by=sort_col, ascending=False)

            df.head(args.num_candidates).to_csv(output_path, index=False)
            logging.info(f"Filtered to top {args.num_candidates} candidates.")
            print(f"✅ Filtered to top {args.num_candidates} candidates and saved to {output_path}")

        print("✅ Property prediction pipeline complete.")

    except Exception as e:
        logging.error(f"Pipeline failed: {e}", exc_info=True)
        print("❌ Prediction error (check log).")

if __name__ == "__main__":
    main()
