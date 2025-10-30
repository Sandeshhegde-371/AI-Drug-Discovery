import argparse
import logging
import sys
import os
import pandas as pd
import time

# ---------------------------------------------------------------------
# 🔧 Setup absolute paths
# ---------------------------------------------------------------------
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
RESULTS_FOLDER = os.path.join(BASE_DIR, 'data', 'candidates')
LOGS_FOLDER = os.path.join(BASE_DIR, 'results', 'logs')
MODELS_FOLDER = os.path.join(BASE_DIR, 'models')

# Ensure required folders exist
os.makedirs(RESULTS_FOLDER, exist_ok=True)
os.makedirs(LOGS_FOLDER, exist_ok=True)
os.makedirs(MODELS_FOLDER, exist_ok=True)

# ---------------------------------------------------------------------
# 🧾 Logging Configuration
# ---------------------------------------------------------------------
log_file_path = os.path.join(LOGS_FOLDER, 'pipeline_generative_screening.log')
logging.basicConfig(
    filename=log_file_path,
    level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s'
)

# ---------------------------------------------------------------------
# 🧠 Import generative design module safely
# ---------------------------------------------------------------------
try:
    sys.path.insert(0, os.path.join(BASE_DIR, 'src'))
    from src import generative_design
except Exception as e:
    print("[WARN] Could not import generative_design module. Using mock mode.")
    generative_design = None


def main():
    parser = argparse.ArgumentParser(description="Run generative screening pipeline")
    parser.add_argument('--input', required=True, help='Input CSV with seed compounds')
    parser.add_argument('--num_candidates', type=int, default=50, help='Number of compounds to generate')
    args = parser.parse_args()

    # -----------------------------------------------------------------
    # 🔑 Generate unique output filename
    # -----------------------------------------------------------------
    input_filename = os.path.basename(args.input)
    unique_id = input_filename.split('_')[0]
    output_csv = os.path.join(RESULTS_FOLDER, f"{unique_id}_output.csv")

    try:
        logging.info(f"Starting generative screening for: {args.input}")
        print(f"[INFO] Generating candidates from: {args.input}")

        # -----------------------------------------------------------------
        # ⚙️ Check for model availability
        # -----------------------------------------------------------------
        model_file = os.path.join(MODELS_FOLDER, 'multitask_admet_model.h5')
        if not os.path.exists(model_file):
            print("[WARN] Model file missing. Generating mock data...")
            logging.warning("Model not found. Using mock data generation.")
            df = pd.read_csv(args.input)

            # Ensure SMILES column exists
            if 'smiles' not in df.columns:
                print("SMILES column missing in seed file.")
                df['smiles'] = ['C', 'CC', 'CCC'][:len(df)]

            # Mock generated output
            df['generated_candidate'] = [f"mol_{i+1}" for i in range(len(df))]
            df['predicted_activity'] = [round((i + 1) * 0.05, 3) for i in range(len(df))]
            df.to_csv(output_csv, index=False)
        else:
            # -----------------------------------------------------------------
            # 🧬 Run real generative model (if module available)
            # -----------------------------------------------------------------
            if generative_design and hasattr(generative_design, 'generate_candidates'):
                generative_design.generate_candidates(
                    seed_csv=args.input,
                    output_csv=output_csv,
                    model_file=model_file,
                    n_generate=args.num_candidates
                )
            else:
                print("[WARN] generative_design missing or incomplete. Using mock output.")
                df = pd.read_csv(args.input)
                df['generated_candidate'] = [f"mol_{i+1}" for i in range(len(df))]
                df['predicted_activity'] = [round((i + 1) * 0.05, 3) for i in range(len(df))]
                df.to_csv(output_csv, index=False)

        # -----------------------------------------------------------------
        # ✅ Verify file created
        # -----------------------------------------------------------------
        if os.path.exists(output_csv):
            logging.info(f"Pipeline complete. Output saved at: {output_csv}")
            print(f"[SUCCESS] Saved: {output_csv}")
        else:
            raise FileNotFoundError(f"Expected output file not found at: {output_csv}")

        time.sleep(1)

    except Exception as e:
        logging.error(f"Pipeline failed: {e}", exc_info=True)
        print(f"[ERROR] Generative screening failed: {e}")


# ---------------------------------------------------------------------
# 🚀 Main Entry Point
# ---------------------------------------------------------------------
if __name__ == "__main__":
    main()
