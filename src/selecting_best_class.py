import pandas as pd
import logging
import os

logging.basicConfig(
    filename='./results/logs/selecting_best_class.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def select_best_model(metrics_csv, output_csv, metric_name='accuracy'):
    """
    Selects the best classification model based on a specified metric.
    Saves the best model's details to a CSV file.
    """
    try:
        if not os.path.exists(metrics_csv):
            logging.error(f'Metrics file not found: {metrics_csv}')
            print('Metrics file missing.')
            return None

        df = pd.read_csv(metrics_csv)
        if metric_name not in df.columns:
            logging.error(f'Metric {metric_name} not found in metrics file.')
            print(f'Metric {metric_name} missing.')
            return None

        best_model_row = df.loc[df[metric_name].idxmax()]
        best_model_df = pd.DataFrame([best_model_row])
        best_model_df.to_csv(output_csv, index=False)
        logging.info(f'Best model based on {metric_name} saved to {output_csv}.')
        print(f'Best model saved to {output_csv} based on {metric_name}.')
        return best_model_df

    except Exception as e:
        logging.error(f'Error selecting best classification model: {e}')
        print('Error occurred during best model selection. See logs.')
        return None

if __name__ == "__main__":
    select_best_model(
        './results/summaries/classification_metrics.csv',
        './results/summaries/best_classification_model.csv',
        metric_name='accuracy'
    )
