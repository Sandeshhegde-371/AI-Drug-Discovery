import tensorflow as tf
from tensorflow.keras import layers, models
import logging
import os

logging.basicConfig(
    filename=os.path.abspath('./results/logs/multitask_model.log'),
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def build_deep_multitask_model(num_features=4, num_tasks=5):
    """
    Build a deep multitask feedforward neural network model.
    Args:
      num_features (int): Input features (e.g., MW, LogP, etc.)
      num_tasks (int): Number of output properties predicted simultaneously
    Returns:
      keras.Model: compiled deep multitask model
    """
    inputs = layers.Input(shape=(num_features,))
    # Deep architecture with dropout regularization
    x = layers.Dense(256, activation='relu')(inputs)
    x = layers.Dropout(0.3)(x)
    x = layers.Dense(128, activation='relu')(x)
    x = layers.Dropout(0.3)(x)
    x = layers.Dense(64, activation='relu')(x)
    x = layers.Dropout(0.2)(x)
    x = layers.Dense(32, activation='relu')(x)

    # Output layer for multitask regression
    outputs = layers.Dense(num_tasks, activation='linear')(x)

    model = models.Model(inputs=inputs, outputs=outputs)
    model.compile(optimizer='adam', loss='mse', metrics=['mae'])

    logging.info("Deep multitask model built and compiled.")
    return model


def save_model(model, path='./models/multitask_admet_model.h5'):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    model.save(path)
    logging.info(f"Multitask model saved to {path}")


if __name__ == "__main__":
    model = build_deep_multitask_model(num_features=4, num_tasks=5)
    save_model(model)
