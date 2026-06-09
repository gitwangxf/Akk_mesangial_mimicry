import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pymol import cmd


def sliding_window_smooth(data: pd.DataFrame, window: int, calc_col: str) -> pd.DataFrame:
    """Smoothes the specified column using a Hamming window convolution."""
    data = data.copy()
    if len(data) >= window:
        kernel = np.hamming(window)
        data[calc_col] = np.convolve(data[calc_col], kernel / kernel.sum(), mode='same')
    return data


def extract_consecutive_numbers(numbers: list) -> list:
    """
    Extracts consecutive number sequences from a sorted list.
    Only retains sequences with a length greater than 5.
    """
    if not numbers:
        return []

    sorted_numbers = sorted(numbers)
    result = []
    current_sequence = [sorted_numbers[0]]

    for i in range(1, len(sorted_numbers)):
        if sorted_numbers[i] == current_sequence[-1] + 1:
            current_sequence.append(sorted_numbers[i])
        else:
            if len(current_sequence) > 5:
                result.append(current_sequence)
            current_sequence = [sorted_numbers[i]]
            
    if len(current_sequence) > 5:
        result.append(current_sequence)

    return result


def process_discotope3_epitopes(csv_dir: str, output_dir: str, window_size: int = 20, threshold: float = 0.7):
    """
    Processes Discotope3 prediction results, smoothes scores, extracts 
    consecutive epitope fragments, and saves them into individual PDB files.
    
    Parameters:
        csv_dir: Directory containing Discotope3 CSV results and matching PDB files.
        output_dir: Directory where the extracted epitope PDB files will be saved.
        window_size: Size of the sliding window for smoothing.
        threshold: Calibrated score threshold to identify epitope residues.
    """
    os.makedirs(output_dir, exist_ok=True)
    files = glob.glob(os.path.join(csv_dir, '*/*.csv'))
    
    print(f"[Discotope3] Found {len(files)} files to process.")

    for file in files:
        filename = os.path.splitext(os.path.basename(file))[0]
        df_discotope = pd.read_csv(file)
        
        try:
            df_smooth = sliding_window_smooth(df_discotope, window_size, 'calibrated_score')
        except Exception as e:
            print(f"Warning: Smoothing failed for {filename} ({e}). Using raw data instead.")
            df_smooth = df_discotope
            
        # Filter residue IDs above the threshold
        passed_resids = list(df_smooth.query(f'calibrated_score > {threshold}')['res_id'])
        epi_range_list = extract_consecutive_numbers(passed_resids)

        if not epi_range_list:
            continue

        # Load the corresponding original PDB file
        pred_pdb_path = file.replace('.csv', '.pdb')
        if not os.path.exists(pred_pdb_path):
            print(f"Error: Matching PDB file not found at {pred_pdb_path}")
            continue

        for epi_range in epi_range_list:
            epi_tag = f"{epi_range[0]}-{epi_range[-1]}"
            pred_epi_path = os.path.join(output_dir, f"{filename}.{epi_tag}.pdb")
            
            # PyMOL structure extraction logic
            cmd.delete("all")
            cmd.load(pred_pdb_path, "pred_PDB")
            selection_cmd = f"resi {epi_tag} and mode pred_PDB"
            cmd.create("pred_epi", selection_cmd)
            cmd.save(pred_epi_path, "pred_epi")


def process_netmhciipan_epitopes(csv_dir: str, pdb_dir: str, output_pdb_dir: str, output_fasta_path: str, rank_threshold: float = 10.0):
    """
    Processes netMHCiiPan prediction results, extracts low-rank epitope fragments 
    as PDB files, and aggregates all epitope sequences into a single FASTA file.
    
    Parameters:
        csv_dir: Directory containing tab-separated netMHCiiPan CSV files.
        pdb_dir: Source directory containing the full-length PDB database.
        output_pdb_dir: Directory where the extracted epitope PDB files will be saved.
        output_fasta_path: File path where the combined FASTA file will be exported.
        rank_threshold: Percentile rank threshold for strong/weak binders.
    """
    os.makedirs(output_pdb_dir, exist_ok=True)
    os.makedirs(os.path.dirname(output_fasta_path), exist_ok=True)
    
    files = glob.glob(os.path.join(csv_dir, '*.csv'))
    print(f"[netMHCiiPan] Found {len(files)} files to process.")
    
    fasta_records = []

    for file in files:
        filename = os.path.splitext(os.path.basename(file))[0]
        df_out = pd.read_csv(file, sep='\t')
        
        # Filter and drop duplicate peptides safely
        df_rank = df_out.query(f'percentile_rank < {rank_threshold}').copy()
        df_rank.drop_duplicates(subset='peptide', inplace=True)
        df_rank.reset_index(drop=True, inplace=True)

        # 1. Extract structural fragments using PyMOL
        pred_pdb_path = os.path.join(pdb_dir, f"{filename}.pdb")
        if os.path.exists(pred_pdb_path):
            for i, row in df_rank.iterrows():
                epi_range = f"{row['start']}-{row['end']}"
                pred_epi_path = os.path.join(output_pdb_dir, f"{filename}_{epi_range}.pdb")
                
                cmd.delete("all")
                cmd.load(pred_pdb_path, "pred_PDB")
                selection_cmd = f"resi {epi_range} and mode pred_PDB"
                cmd.create("pred_epi", selection_cmd)
                cmd.save(pred_epi_path, "pred_epi")
        else:
            print(f"Warning: Database PDB file not found at {pred_pdb_path}. Skipping structural extraction.")

        # 2. Build FASTA records
        for i, row in df_rank.iterrows():
            seq_id = f"{row['allele']}_{filename}"
            fasta_records.append((seq_id, row['peptide']))

    # 3. Write records to a FASTA file
    if fasta_records:
        with open(output_fasta_path, 'w') as f:
            for seq_id, peptide in fasta_records:
                f.write(f">{seq_id}\n{peptide}\n")
        print(f"Successfully exported FASTA file to: {output_fasta_path}")


# ==========================================
# Execution Pipeline Example
# ==========================================
if __name__ == "__main__":
    
    # --- Task 1: Execute Discotope3 Pipeline ---
    DISCOTOPE_CSV_DIR = "./data/discotope_results" 
    DISCOTOPE_OUT_DIR = "./data/query_epitopes"
    
    process_discotope3_epitopes(
        csv_dir=DISCOTOPE_CSV_DIR, 
        output_dir=DISCOTOPE_OUT_DIR
    )
    
    # --- Task 2: Execute netMHCiiPan Pipeline ---
    NETMHCIIPAN_CSV_DIR = "./data/netmhciipan_results"
    GLOBAL_PDB_DATABASE = "./data/foldseek_db"
    NETMHC_OUT_PDB_DIR  = "./data/netmhc_epi_pdb"
    NETMHC_OUT_FASTA    = "./data/output_fasta/results.fasta"
    
    process_netmhciipan_epitopes(
        csv_dir=NETMHCIIPAN_CSV_DIR,
        pdb_dir=GLOBAL_PDB_DATABASE,
        output_pdb_dir=NETMHC_OUT_PDB_DIR,
        output_fasta_path=NETMHC_OUT_FASTA
    )
