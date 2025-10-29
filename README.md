🧬 AI Drug Discovery Platform

An end-to-end AI-powered Drug Discovery Platform that integrates machine learning, Flask-based APIs, and a React frontend to streamline compound screening, candidate generation, and predictive analysis for drug design.

🚀 Features

✅ Generative Candidate Screening
Automatically generates new molecular candidates using seed compound data and predictive models.

✅ ADMET Prediction Pipeline
Evaluates generated molecules for Absorption, Distribution, Metabolism, Excretion, and Toxicity (ADMET).

✅ Interactive Dashboard (React)
Frontend built using React for uploading files, triggering pipelines, and visualizing results.

✅ Flask Backend API
Backend handles file uploads, invokes the ML pipeline, and returns results dynamically.

✅ Seamless Integration with ML Models
Supports TensorFlow/Keras .h5 model loading for multitask ADMET property prediction.

🏗️ Project Structure
AI-Drug-Discovery/
│
├── backend_pipeline/
│   ├── pipeline_generative_screening.py    # Main ML pipeline
│   ├── src/
│   │   └── generative_design.py            # Candidate generation logic
│
├── frontend/
│   ├── src/
│   │   ├── App.js / App.tsx                # React app entry
│   │   ├── components/                     # React UI components
│   ├── public/
│   └── package.json
│
├── data/
│   ├── candidates/                         # Generated candidates output
│
├── models/
│   └── multitask_admet_model.h5            # Pretrained model (ignored in .gitignore)
│
├── results/
│   ├── logs/                               # Log files
│
├── uploads/                                # Uploaded CSV files
│
├── app.py                                  # Flask backend server
├── .gitignore
└── README.md

⚙️ Setup Instructions
🧩 1. Clone the Repository
git clone https://github.com/<your-username>/AI-Drug-Discovery.git
cd AI-Drug-Discovery

🐍 2. Set Up the Backend (Flask)
cd backend_pipeline
python -m venv venv
source venv/bin/activate     # On Windows: venv\Scripts\activate
pip install -r requirements.txt

▶️ Run Flask Server
python app.py


The Flask API will start at:
👉 http://127.0.0.1:5000

⚛️ 3. Set Up the Frontend (React)
cd ../frontend
npm install
npm start


The React app will start at:
👉 http://localhost:3000

📤 Usage Flow

Upload a CSV file containing seed compounds (e.g., compound_name, SMILES).

Click Run Pipeline to start the generative screening process.

The backend runs your ML pipeline (pipeline_generative_screening.py).

Once completed, results are stored in:

data/candidates/<unique_id>_output.csv


View or download results directly from the dashboard.

🧠 Technologies Used
Layer	Tech Stack
Frontend	React, TailwindCSS
Backend	Flask, Python
Machine Learning	TensorFlow, Pandas, NumPy, scikit-learn
Data Storage	CSV (local), easily extendable to PostgreSQL
Version Control	Git + GitHub
🧪 Example Input File
compound_name,SMILES
Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
Paracetamol,CC(=O)NC1=CC=C(O)C=C1

🧾 Example Output File
compound_name,SMILES,generated_candidate,predicted_activity
Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O,mol_0,0.123
Paracetamol,CC(=O)NC1=CC=C(O)C=C1,mol_1,0.256

🛡️ .gitignore Highlights

The project ignores:

node_modules/, __pycache__/

uploads/, data/candidates/

models/, results/logs/

.env, virtual environments, build artifacts

This ensures your repo stays clean and lightweight.

👨‍💻 Author

Sathwik N H ,Sandesh Hegde,Sinchan M S
AI & ML Enthusiasts 


⭐ Contribute

Contributions are welcome!
To contribute:

Fork this repo

Create your feature branch (git checkout -b feature-xyz)

Commit changes (git commit -m "Added xyz")

Push (git push origin feature-xyz)

Open a Pull Request

📜 License

This project is licensed under the MIT License — feel free to use, modify, and distribute with attribution.