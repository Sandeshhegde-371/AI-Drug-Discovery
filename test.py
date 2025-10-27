import pandas as pd

paths = [
    './data/processed/normalized_data_features.csv',
    './data/processed/logd_reg_final_data.csv',
    './data/processed/water_sol_reg_final_data.csv'
]

for path in paths:
    try:
        df = pd.read_csv(path)
        print(f"{path}: {df.columns.tolist()}")
    except Exception as e:
        print(f"Error reading {path}: {e}")
