from flask import Flask, render_template, request, redirect, url_for, send_from_directory, flash
import os
import uuid
import subprocess

app = Flask(__name__)
app.secret_key = 'your_secret_key_here'  # Replace with your own secret key for session security

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
UPLOAD_FOLDER = os.path.join(BASE_DIR, 'uploads')
RESULTS_FOLDER = os.path.join(BASE_DIR, 'data', 'candidates')

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)

@app.route('/')
def index():
    # Render homepage with form for user inputs
    return render_template('index.html')

@app.route('/run_pipeline', methods=['POST'])
def run_pipeline():
    # Get form inputs
    uploaded_file = request.files.get('input_file')
    pipeline_choice = request.form.get('pipeline_choice')
    num_candidates = request.form.get('num_candidates', type=int)

    if not uploaded_file or uploaded_file.filename == '':
        flash('Please upload an input CSV file.')
        return redirect(url_for('index'))

    if pipeline_choice not in ['generative_screening', 'virtual_screening', 'predict_properties', 'prioritize_candidates']:
        flash('Please select a valid pipeline.')
        return redirect(url_for('index'))

    # Save uploaded file
    unique_id = str(uuid.uuid4())
    filename = f"{unique_id}_{uploaded_file.filename}"
    save_path = os.path.join(UPLOAD_FOLDER, filename)
    uploaded_file.save(save_path)

    # Prepare command to run pipeline
    # This example assumes you have created separate scripts for each pipeline under 'pipelines' folder
    pipeline_script_map = {
        'generative_screening': 'pipeline_generative_screening.py',
        'virtual_screening': 'pipeline_virtual_screening.py',
        'predict_properties': 'pipeline_predict_properties.py',
        'prioritize_candidates': 'pipeline_prioritize_candidates.py'
    }

    script_to_run = pipeline_script_map[pipeline_choice]

    # Parameters to pass to the script can be adjusted
    cmd = [
        'python', os.path.join('pipelines', script_to_run),
        '--input', save_path
        ]
    if pipeline_choice in ['generative_screening', 'virtual_screening', 'prioritize_candidates']:
        if num_candidates:  # Only add if actually set
            cmd += ['--num_candidates', str(num_candidates)]

    # Run pipeline synchronously or asynchronously depending on requirements
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        flash('Pipeline executed successfully.')

        # Assuming pipeline saves output CSV to RESULTS_FOLDER with unique_id prefix
        output_filename = f'{unique_id}_output.csv'
        output_path = os.path.join(RESULTS_FOLDER, output_filename)

        if not os.path.exists(output_path):
            flash('Pipeline completed but output file not found.')
            return redirect(url_for('index'))

        return redirect(url_for('download_file', filename=output_filename))

    except subprocess.CalledProcessError as e:
        flash(f'Pipeline execution failed: {e.stderr}')
        return redirect(url_for('index'))

@app.route('/download/<filename>')
def download_file(filename):
    return send_from_directory(RESULTS_FOLDER, filename, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
