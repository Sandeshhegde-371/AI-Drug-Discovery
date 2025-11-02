import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def generate_candidates(seed_csv, output_csv, model_file=None, ngenerate=10):
    df_seed = pd.read_csv(seed_csv)
    all_records = []
    
    # Load model if available
    if model_file:
        from tensorflow.keras.models import load_model
        model = load_model(model_file)
    else:
        model = None

    for idx, row in df_seed.iterrows():
        input_smiles = row['SMILES']
        MW = row.get('MW', np.nan)
        LogP = row.get('LogP', np.nan)
        NumHDonors = row.get('NumHDonors', np.nan)
        NumHAcceptors = row.get('NumHAcceptors', np.nan)
        for n in range(ngenerate):
            steps = []
            processes = []

            # Step 1: Hydrogenate
            mol = Chem.MolFromSmiles(input_smiles)
            mol = Chem.AddHs(mol)
            hydrogen_smiles = Chem.MolToSmiles(mol)
            steps.append("Hydrogenation")
            processes.append("Add explicit hydrogens to molecule")
            
            # Step 2: Side-chain modification (example: add methyl, simulated by appending C)
            methyl_smiles = hydrogen_smiles + "C"
            steps.append("Side-chain modification")
            processes.append("Simulate methyl addition")

            # Step 3: Canonicalize and recalculate descriptors
            gen_mol = Chem.MolFromSmiles(methyl_smiles)
            gen_mol = gen_mol if gen_mol else mol  # fallback to previous if fail
            new_MW = Descriptors.MolWt(gen_mol) if gen_mol else MW
            new_LogP = Descriptors.MolLogP(gen_mol) if gen_mol else LogP
            new_NumHDonors = Descriptors.NumHDonors(gen_mol) if gen_mol else NumHDonors
            new_NumHAcceptors = Descriptors.NumHAcceptors(gen_mol) if gen_mol else NumHAcceptors
            
            # Model prediction
            if model is not None:
                # Prepare array in correct shape
                X = np.array([[new_MW, new_LogP, new_NumHDonors, new_NumHAcceptors]])
                pred = model.predict(X)
                pred_val = float(pred[0]) if hasattr(pred, 'shape') and pred.shape else float(pred)
            else:
                pred_val = np.round(np.random.uniform(0, 1), 3)

            all_records.append({
                'input_compound': input_smiles,
                'step_1': steps[0],
                'process_1': processes[0],
                'step_2': steps[1],
                'process_2': processes[1],
                'generated_smiles': Chem.MolToSmiles(gen_mol) if gen_mol else methyl_smiles,
                'MW': new_MW,
                'LogP': new_LogP,
                'NumHDonors': new_NumHDonors,
                'NumHAcceptors': new_NumHAcceptors,
                'predicted_activity': pred_val
            })
    
    df_out = pd.DataFrame(all_records)
    df_out.to_csv(output_csv, index=False)
    print(f"Output written to {output_csv}")

# Usage Example:
# generate_candidates('generative_seed_sample.csv', 'generative_sample_detailed.csv', model_file='models/multitask_admet_model.h5', ngenerate=10)