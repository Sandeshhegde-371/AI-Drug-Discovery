import logging
import sys
from pathlib import Path

# Ensure project root (parent of this pipelines directory) is on sys.path so 'src' can be imported.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import data_prep, merging

logging.basicConfig(
    filename='./results/logs/pipeline_data_prep_and_merge.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def run_pipeline():
    logging.info("Pipeline started: Data Preparation and Merging")
    norm_df = data_prep.prepare_main_data()
    if norm_df is None:
        logging.error("Data preparation failed. Exiting pipeline.")
        return
    
    input_files = [
        './data/processed/normalized_data_features.csv',
        './data/processed/logd_reg_final_data.csv',
        './data/processed/water_sol_reg_final_data.csv',
        # Add additional processed dataset paths as required
    ]
    logging.info("Starting dataset merging...")
    merged = merging.merge_datasets(input_files, './data/processed/merged_multitask_data.csv')
    if merged is None:
        logging.error("Dataset merging failed.")
    else:
        logging.info("Pipeline completed successfully.")

if __name__ == "__main__":
    run_pipeline()
