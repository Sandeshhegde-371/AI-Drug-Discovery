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
import sys

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
    Accepts multipart/form-data:
    - 'file': uploaded CSV file
    - 'pipeline': which pipeline to run
    - 'numCandidates': optional integer
    """
    try:
        # Check for file upload
        if 'file' not in request.files:
            return jsonify({"error": "No file uploaded"}), 400

        file = request.files['file']
        if file.filename == '':
            return jsonify({"error": "Empty filename"}), 400

        # Save uploaded file
        input_file = os.path.join(UPLOAD_FOLDER, f"{uuid.uuid4()}_{file.filename}")
        file.save(input_file)
        print(f"✅ Saved input file: {input_file}")

        # Get form data
        pipeline = request.form.get('pipeline', 'generative')
        num_candidates = int(request.form.get('numCandidates', 50))

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
            sys.executable,  # uses same Python environment
            script_path,
            '--input', input_file,
            '--num_candidates', str(num_candidates)
        ]


        # ✅ Run inside project root so relative imports (like from src import generative_design) work
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=BASE_DIR)

        print("PIPELINE STDOUT:", result.stdout)
        print("PIPELINE STDERR:", result.stderr)


        # ----------------- Locate Output File -----------------
        output_files = glob.glob(os.path.join(BASE_DIR, 'data', 'candidates', '*_output.csv'))
        if not output_files:
            return jsonify({"error": "Pipeline completed but no output file found"}), 500

        output_file = max(output_files, key=os.path.getctime)
        print(f"✅ Output file: {output_file}")

        # ----------------- Read & Return Results -----------------
        df = pd.read_csv(output_file)
        results_json = df.to_dict(orient='records')

        return jsonify({
            "status": "success",
            "results": results_json,
            "message": f"Pipeline executed successfully with {len(results_json)} results."
        })

    except Exception as e:
        print("❌ Error running pipeline:", e)
        return jsonify({"error": str(e)}), 500


# ----------------- MAIN ENTRY -----------------
if __name__ == '__main__':
    app.run(debug=True)
