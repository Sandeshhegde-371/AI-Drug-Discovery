import os
import logging

logging.basicConfig(
    filename='./results/logs/utils.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)

def check_file_exists(filepath):
    try:
        if not os.path.exists(filepath):
            logging.error(f'File not found: {filepath}')
            return False
        logging.info(f'File exists: {filepath}')
        return True
    except Exception as e:
        logging.error(f'Error checking file {filepath}: {e}')
        return False

def ensure_directory_exists(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)
            logging.info(f'Directory created: {directory}')
        else:
            logging.info(f'Directory exists: {directory}')
    except Exception as e:
        logging.error(f'Error creating directory {directory}: {e}')

def log_and_print(message, level='info'):
    try:
        print(message)
        if level == 'info':
            logging.info(message)
        elif level == 'warning':
            logging.warning(message)
        elif level == 'error':
            logging.error(message)
        else:
            logging.debug(message)
    except Exception as e:
        print(f'Error during logging: {e}')
        logging.error(f'Logging failure: {e}')

# Example usage
if __name__ == "__main__":
    log_and_print('Testing utility functions.')
    ensure_directory_exists('./results/logs')
    file_check = check_file_exists('./data/processed/normalized_data.csv')
    log_and_print(f'File exists: {file_check}')
