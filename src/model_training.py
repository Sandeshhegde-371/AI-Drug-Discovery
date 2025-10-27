import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import classification_report, accuracy_score, precision_score, recall_score, f1_score, r2_score, mean_squared_error
import tensorflow as tf
from tensorflow import keras
import joblib
import logging
import os

# Set up logging
logging.basicConfig(
    filename='./results/logs/model_training.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def train_models():
    try:
        train_path = './data/splits/train_data.csv'
        test_path = './data/splits/test_data.csv'
        if not os.path.exists(train_path) or not os.path.exists(test_path):
            logging.error('Train or test splits not found.')
            print('Train/test data files missing.')
            return

        train_df = pd.read_csv(train_path)
        test_df = pd.read_csv(test_path)
        logging.info('Train and test data loaded.')

        feature_cols = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
        X_train = train_df[feature_cols].values
        X_test = test_df[feature_cols].values

        # Classification
        le = LabelEncoder()
        y_train = le.fit_transform(train_df['activity_class'])
        y_test = le.transform(test_df['activity_class'])
        clf = RandomForestClassifier(n_estimators=100, random_state=42)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        acc = accuracy_score(y_test, y_pred)
        prec = precision_score(y_test, y_pred, average='weighted', zero_division=0)
        rec = recall_score(y_test, y_pred, average='weighted', zero_division=0)
        f1 = f1_score(y_test, y_pred, average='weighted', zero_division=0)
        report = classification_report(y_test, y_pred)
        print("Classification Report:\n", report)
        logging.info('RandomForestClassifier trained and evaluated.')

        # Regression
        y_train_reg = train_df['pIC50'].values
        y_test_reg = test_df['pIC50'].values
        regr = RandomForestRegressor(n_estimators=100, random_state=42)
        regr.fit(X_train, y_train_reg)
        y_pred_reg = regr.predict(X_test)
        r2 = r2_score(y_test_reg, y_pred_reg)
        mse = mean_squared_error(y_test_reg, y_pred_reg)
        print(f"Regression Metrics: R2 Score={r2}, MSE={mse}")
        logging.info(f'Regression metrics: R2 Score={r2}, MSE={mse}.')

        # Neural network
        nn_model = keras.Sequential([
            keras.layers.Dense(64, activation='relu', input_shape=(len(feature_cols),)),
            keras.layers.Dense(32, activation='relu'),
            keras.layers.Dense(len(set(y_train)), activation='softmax')
        ])
        nn_model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])
        nn_model.fit(X_train, y_train, epochs=50, batch_size=16, validation_data=(X_test, y_test))
        logging.info('Neural network model trained.')

        # Saving models
        os.makedirs('./models/', exist_ok=True)
        joblib.dump(clf, './models/activity_rf_model.pkl')
        joblib.dump(regr, './models/pIC50_rf_model.pkl')
        nn_model.save('./models/activity_nn_model.h5')
        logging.info('Models saved: activity_rf_model.pkl, pIC50_rf_model.pkl, activity_nn_model.h5')
        print("Models saved successfully.")

        # Save metrics to CSV for pipeline consumption
        os.makedirs('./results/summaries/', exist_ok=True)
        clf_metrics = pd.DataFrame([{
            'model': 'RandomForest',
            'accuracy': acc,
            'precision': prec,
            'recall': rec,
            'f1': f1
        }])
        clf_metrics.to_csv('./results/summaries/classification_metrics.csv', index=False)

        reg_metrics = pd.DataFrame([{
            'model': 'RandomForestRegressor',
            'r2_score': r2,
            'mse': mse
        }])
        reg_metrics.to_csv('./results/summaries/regression_metrics.csv', index=False)

        print("Metrics files saved.")

    except Exception as e:
        logging.error(f'Error in model training: {e}')
        print('An error occurred during model training. See logs for details.')

if __name__ == "__main__":
    train_models()
