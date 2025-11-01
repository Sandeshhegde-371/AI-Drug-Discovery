from flask import Flask, request, jsonify
try:
    from flask_cors import CORS
except Exception:
    # Flask-CORS not available; provide a no-op CORS function so the app can run without the package.
    def CORS(app, **kwargs):
        return None

import os
import uuid
import subprocess
import glob
import pandas as pd
import csv

app = Flask(__name__)
CORS(app)

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
UPLOAD_FOLDER = os.path.join(BASE_DIR, 'uploads')
RESULTS_FOLDER = os.path.join(BASE_DIR, 'results')
PIPELINE_FOLDER = os.path.join(BASE_DIR, 'pipelines')

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)


@app.route('/')
def home():
    return "✅ Flask backend is running successfully!"

# ----------------- FRONTEND API -----------------
@app.route('/api/run_pipeline', methods=['POST'])
def api_run_pipeline():
    """
    Endpoint called from React frontend to run the pipeline.
    Accepts JSON data with pipeline type, params, and optionally uploaded data.
    Returns real results from generated CSV.
    """
    try:
        data = request.get_json()
        pipeline = data.get('pipeline', 'generative')
        params = data.get('params', {})
        num_candidates = int(params.get('numCandidates', 50))

        # Save input data if provided
        uploaded_data = data.get('data')
        input_file = os.path.join(UPLOAD_FOLDER, f"{uuid.uuid4()}_input.csv")

        if uploaded_data:
            keys = uploaded_data[0].keys()
            with open(input_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=keys)
                writer.writeheader()
                writer.writerows(uploaded_data)
        else:
            return jsonify({"error": "No input data provided"}), 400

        # ----------------- Pipeline Map -----------------
        pipeline_script_map = {
            'generative': 'pipeline_generative_screening.py',
            'virtual': 'pipeline_virtual_screening.py',
            'properties': 'pipeline_predict_properties.py',
            'prioritize': 'pipeline_prioritize_candidates.py'
        }

        script = pipeline_script_map.get(pipeline)
        if not script:
            return jsonify({"error": f"Unknown pipeline: {pipeline}"}), 400

        script_path = os.path.join(PIPELINE_FOLDER, script)
        if not os.path.exists(script_path):
            return jsonify({"error": f"Pipeline script not found: {script_path}"}), 500

        # ----------------- Run Pipeline -----------------
        cmd = [
            'python', script_path,
            '--input', input_file,
            '--num_candidates', str(num_candidates)
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        print("PIPELINE STDOUT:", result.stdout)
        print("PIPELINE STDERR:", result.stderr)

        # ----------------- Locate Output File -----------------
        output_files = glob.glob(os.path.join(BASE_DIR, 'data', 'candidates', '*_output.csv'))
        if not output_files:
            return jsonify({"error": "Pipeline completed but no output file found"}), 500

        output_file = max(output_files, key=os.path.getctime)

        # ----------------- Read & Return Real Results -----------------
        df = pd.read_csv(output_file)
        results_json = df.to_dict(orient='records')

        return jsonify({
            "status": "success",
            "results": results_json,
            "message": f"Pipeline executed successfully with {len(results_json)} results."
        })

    except Exception as e:
        print("Error running pipeline:", e)
        return jsonify({"error": str(e)}), 500


# ----------------- MAIN ENTRY -----------------
if __name__ == '__main__':
    app.run(debug=True)
