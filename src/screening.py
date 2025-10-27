import pandas as pd
import logging
import os

logging.basicConfig(
    filename='./results/logs/screening.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def prioritize_candidates(pred_csv, output_csv, priority_metric='pIC50', top_n=10):
    """
    Ranks predicted candidates by the given metric and saves top N hits.
    """
    try:
        if not os.path.exists(pred_csv):
            logging.error(f'Predicted candidates file not found: {pred_csv}')
            print('Predicted candidates file missing.')
            return None

        df = pd.read_csv(pred_csv)
        if priority_metric not in df.columns:
            logging.error(f'Priority metric {priority_metric} not in candidates data.')
            print(f'Priority metric {priority_metric} not found.')
            return None

        ranked = df.sort_values(by=priority_metric, ascending=False).head(top_n)
        ranked.to_csv(output_csv, index=False)
        logging.info(f'Top {top_n} candidates ranked by {priority_metric} and saved to {output_csv}.')
        print(f'Top {top_n} prioritized candidates written to {output_csv}.')
        return ranked

    except Exception as e:
        logging.error(f'Error during candidate screening: {e}')
        print('Error in prioritizing candidates. See logs for details.')
        return None

if __name__ == "__main__":
    prioritize_candidates(
        './data/candidates/predicted_candidates.csv',
        './data/candidates/top10_hits_for_design.csv',
        priority_metric='pIC50',
        top_n=10
    )
