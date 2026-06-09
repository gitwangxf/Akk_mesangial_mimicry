import os
import argparse
import pandas as pd
from pymol import cmd


def calculate_average_b_factor(selection: str) -> float:
    """
    Calculates the average B-factor (or pLDDT score) for a given PyMOL selection.
    Utilizes PyMOL atom object attributes natively for better efficiency.
    """
    try:
        atoms = cmd.get_model(selection).atom
        if not atoms:
            return 0.0
        b_factors = [atom.b for atom in atoms]
        return sum(b_factors) / len(b_factors)
    except Exception as e:
        print(f"Warning: Failed to calculate B-factors for selection '{selection}': {e}")
        return 0.0


def calculate_structural_alignment(
    akk_dir: str, 
    target_dir: str,
    akk_id: str, 
    target_id: str, 
    qstart: int, 
    qend: int, 
    tstart: int, 
    tend: int
) -> tuple:
    """
    Loads structural files, extracts interacting local residue segments, 
    calculates localized/global pLDDT scores, and registers structural RMSD.
    
    Parameters:
        akk_dir: Base directory path containing the query PDB files.
        target_dir: Base directory path containing the target PDB files.
        akk_id: File identifier name for the microbial query structure.
        target_id: File identifier name for the mesangial target structure.
        qstart, qend: Sequence alignment range markers for the query.
        tstart, tend: Sequence alignment range markers for the target.
    """
    # Build robust structural paths across systems
    akk_path = os.path.join(akk_dir, f"{akk_id}.pdb")
    target_path = os.path.join(target_dir, f"{target_id}.pdb")
    
    # Initialize PyMOL workspace
    cmd.delete("all")
    
    # Validate structure existence before loading
    if not os.path.exists(target_path) or not os.path.exists(akk_path):
        print(f"Warning: Missing files for pair {akk_id} <-> {target_id}. Skipping.")
        return 0.0, 0.0, 0.0, -1.0

    # 1. Load target and slice aligned epitope/domain segment
    cmd.load(target_path, "target_obj")
    target_range = f"{tstart}-{tend}"
    target_selection = f"resi {target_range} and model target_obj"
    cmd.create("target_align", target_selection)

    # Compute pLDDT metrics for the target
    plddt_align = calculate_average_b_factor("target_align")
    plddt_all = calculate_average_b_factor("target_obj")

    # 2. Load query and slice corresponding localized region
    cmd.load(akk_path, "akk_obj")
    plddt_all_akk = calculate_average_b_factor("akk_obj")
    akk_range = f"{qstart}-{qend}"
    akk_selection = f"resi {akk_range} and model akk_obj"
    cmd.create("akk_align", akk_selection)

    # 3. Structural alignment and calculation of alpha-carbon RMSD
    rmsd_align_ca = -1.0
    try:
        # Preferred structural alignment tool for sequence-independent matching
        result = cmd.cealign("akk_align", "target_align", object="aln")
        rmsd_align_ca = result['RMSD']
    except Exception:
        try:
            # Fallback conventional sequence-dependent alignment method
            result = cmd.align("akk_align", "target_align", object="aln")
            rmsd_align_ca = result[0]
        except Exception:
            print(f"Error: Alignment failed between {akk_id} and {target_id}")
            rmsd_align_ca = -1.0

    return plddt_align, plddt_all, plddt_all_akk, rmsd_align_ca


def main():
    # Command-line execution setup
    parser = argparse.ArgumentParser(description="Calculate regional structural metrics and alignment RMSD values.")
    parser.add_argument('--akk_dir', type=str, required=True, help="Directory path to Akk query structures.")
    parser.add_argument('--target_dir', type=str, required=True, help="Directory path to target mesangial structures.")
    parser.add_argument('--out_path', type=str, required=True, help="Path to the input Foldseek blast-like result file.")
    args = parser.parse_args()

    # Define structured tabular names matching standard Foldseek out-format specifications
    foldseek_cols = [
        "query", "target", "fident", "alnlen", "mismatch", "gapopen",
        "qstart", "qend", "tstart", "tend", "evalue", "bits"
    ]

    print(f"Reading structural match outputs from: {args.out_path}")
    df_out = pd.read_csv(args.out_path, sep='\t', header=None, names=foldseek_cols)

    # Perform structural calculation mapping loop row-by-row
    print("Processing structures and extracting metric features...")
    df_out['cal'] = df_out.apply(
        lambda row: calculate_structural_alignment(
            args.akk_dir, args.target_dir, row['query'], row['target'], 
            int(row['qstart']), int(row['qend']), int(row['tstart']), int(row['tend'])
        ), axis=1
    )

    # Expand calculations into distinctive feature matrix parameters
    df_out['plddt_align'] = df_out['cal'].apply(lambda x: x[0])
    df_out['plddt_all'] = df_out['cal'].apply(lambda x: x[1])
    df_out['plddt_all_akk'] = df_out['cal'].apply(lambda x: x[2])
    df_out['rmsd_align_ca'] = df_out['cal'].apply(lambda x: x[3])

    # Remove temporary coordinate mappings
    df_out.drop(columns=['cal'], inplace=True)

    # Format output filename target logically and save values safely
    result_path = args.out_path.replace('out', 'rmsd_all.csv')
    df_out.to_csv(result_path, index=False)
    print(f"Data mapping process complete! Final records output safely to: {result_path}")


if __name__ == "__main__":
    main()
