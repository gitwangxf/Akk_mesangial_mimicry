import os
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices


def parse_mouse_uniprot_id(target_name: str) -> str:
    """Extracts the UniProt ID from a standard target string format (e.g., 'prefix-ID-suffix')."""
    if pd.isna(target_name) or '-' not in str(target_name):
        return ""
    return str(target_name).split('-')[1]


def parse_akk_query_id(query_name: str) -> str:
    """Converts query IDs into a standardized uppercase format (e.g., 'Akk_1234')."""
    if pd.isna(query_name) or '_' not in str(query_name):
        return ""
    parts = str(query_name).split('_')
    return f"{parts[0].upper()}_{parts[1]}"


def run_pairwise_global_alignment(seq1: str, seq2: str) -> tuple:
    """
    Performs global sequence alignment using the modern Biopython PairwiseAligner API.
    Uses the BLOSUM62 matrix with standard gap penalties (-10 open, -0.5 extend).
    """
    if pd.isna(seq1) or pd.isna(seq2) or str(seq1) == 'nan' or str(seq2) == 'nan':
        return "", ""
    
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    aligner.mode = 'global'
    
    alignments = aligner.align(str(seq1), str(seq2))
    if not alignments:
        return "", ""
        
    # Extract aligned sequences with gaps inserted (-)
    best_alignment = alignments[0]
    return best_alignment.seqA, best_alignment.seqB


def find_letter_indices(aligned_seq: str, start: int, end: int) -> tuple:
    """
    Maps 1-based residue numbers of the original sequence to their 0-based index 
    positions in the gap-containing aligned sequence.
    """
    if not aligned_seq:
        return None, None
        
    # Gather alignment string indices that are actual amino acid characters (not gaps)
    letter_indices = [i for i, char in enumerate(aligned_seq) if char.isalpha()]
    
    if not letter_indices:
        return None, None
        
    # Validate bound constraints
    if start < 1 or end > len(letter_indices) or start > end:
        return None, None

    # start-1 matches 0-based conversion. end maps to the slice position safely.
    epi_start = letter_indices[start - 1]
    epi_end = letter_indices[end - 1] + 1  # Inclusive range boundary for slicing
    return epi_start, epi_end


def extract_epitope_slice(aligned_seq_tuple: tuple, index_tuple: tuple) -> tuple:
    """Extracts the structural structural alignment segment belonging to the epitope slice."""
    if not aligned_seq_tuple or not index_tuple or index_tuple[0] is None:
        return "", ""
        
    seq1_aligned, seq2_aligned = aligned_seq_tuple
    start_idx, end_idx = index_tuple
    
    epi_seq1 = seq1_aligned[start_idx:end_idx]
    epi_seq2 = seq2_aligned[start_idx:end_idx]
    return epi_seq1, epi_seq2


def calculate_blosum_similarity(aligned_seqA: str, aligned_seqB: str, matrix_name: str = 'BLOSUM62') -> float:
    """
    Calculates the BLOSUM62-weighted biological similarity percentage 
    between two pre-aligned sequence segments, skipping gap comparisons.
    """
    if not aligned_seqA or not aligned_seqB:
        return 0.0

    matrix = substitution_matrices.load(matrix_name)
    similar_positions = 0
    total_positions = 0

    for a, b in zip(aligned_seqA, aligned_seqB):
        # Skip position if either contains an alignment gap
        if a == '-' or b == '-':
            continue

        total_positions += 1
        if a == b:
            similar_positions += 1
        # Check if the amino acid substitution score is positive (physicochemically conservative)
        elif matrix.get((a, b), matrix.get((b, a), -1)) > 0:
            similar_positions += 1

    if total_positions == 0:
        return 0.0
        
    return (similar_positions / total_positions) * 100


def run_homolog_screening_pipeline(
    manual_selected_csv: str,
    mouse_id_mapping_tsv: str,
    human_id_mapping_tsv: str,
    akk_metadata_csv: str,
    output_csv: str,
    similarity_threshold: float = 0.0
):
    """
    Main orchestration function to merge mouse and human data, align homologs,
    map epitope regions, calculate alignment similarities, and filter candidate proteins.
    """
    print("Step 1: Loading input datasets...")
    df_info = pd.read_csv(manual_selected_csv)
    df_mouse_id = pd.read_csv(mouse_id_mapping_tsv, sep='\t')
    df_human_id = pd.read_csv(human_id_mapping_tsv, sep='\t')
    df_akk_info = pd.read_csv(akk_metadata_csv)

    # Step 2: Extract UniProt IDs and merge Mouse properties
    df_info['mouse_uniprotid'] = df_info['target'].apply(parse_mouse_uniprot_id)
    df_merge = pd.merge(
        df_info, 
        df_mouse_id[['Entry', 'Entry Name', 'Protein names', 'Sequence']],
        left_on='mouse_uniprotid', right_on='Entry', how='left'
    )
    df_merge.rename(columns={'Sequence': 'Sequence_mouse'}, inplace=True)

    # Step 3: Match Human Ortholog sequence targets 
    df_merge['Entry_name_human'] = df_merge['Entry Name'].apply(
        lambda x: str(x).replace('MOUSE', 'HUMAN') if pd.notna(x) else ""
    )
    df_merge = pd.merge(
        df_merge, 
        df_human_id[['Entry Name', 'Sequence']], 
        left_on='Entry_name_human', right_on='Entry Name', how='left'
    )
    df_merge.rename(columns={'Sequence': 'Sequence_human'}, inplace=True)

    # Step 4: Perform Global Multi-sequence Alignment
    print("Step 2: Aligning Full-length Homologs...")
    df_merge['global_align'] = df_merge.apply(
        lambda row: run_pairwise_global_alignment(row['Sequence_mouse'], row['Sequence_human']), axis=1
    )

    # Step 5: Slice Specific Epitope Coordinates
    print("Step 3: Correlating Epitope regional sequences...")
    df_merge['epi_index'] = df_merge.apply(
        lambda row: find_letter_indices(row['global_align'][0], int(row['tstart']), int(row['tend'])) 
        if pd.notna(row['tstart']) and pd.notna(row['tend']) else (None, None), axis=1
    )
    
    df_merge['epi_align'] = df_merge.apply(
        lambda row: extract_epitope_slice(row['global_align'], row['epi_index']), axis=1
    )

    # Step 6: Score Similarities using Scoring Matrices
    print("Step 4: Evaluating scoring matrix similarity ratios...")
    df_merge['global_similarity'] = df_merge['global_align'].apply(
        lambda x: calculate_blosum_similarity(x[0], x[1])
    )
    df_merge['epi_similarity'] = df_merge['epi_align'].apply(
        lambda x: calculate_blosum_similarity(x[0], x[1])
    )

    # Step 7: Map Back to AKK Microbe Structural Source Annotations
    df_akk_info['qseqid_clean'] = df_akk_info['qseqid_x'].apply(
        lambda x: str(x).split('.')[0] if pd.notna(x) else ""
    )
    df_merge['Akk_id'] = df_merge['query'].apply(parse_akk_query_id)
    
    df_final = pd.merge(
        df_merge, 
        df_akk_info[['qseqid_clean', 'sseqid_x', 'tag']], 
        left_on='Akk_id', right_on='qseqid_clean', how='left'
    )
    df_final.rename(columns={'sseqid_x': 'Akk_uniprotid'}, inplace=True)

    # Step 8: Apply Screening Filter Thresholds
    if similarity_threshold > 0:
        df_final = df_final[df_final['epi_similarity'] >= similarity_threshold]
        print(f"Filtered results using threshold: >= {similarity_threshold}% epitope similarity.")

    # Drop redundant merge artifacts 
    if 'Entry Name_y' in df_final.columns:
        df_final.drop(columns=['Entry Name_y'], inplace=True)

    # Step 9: Save output results
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df_final.to_csv(output_csv, index=False)
    print(f"Pipeline finished! Saved file output securely to: {output_csv}")


# ==========================================
# Script Execution Configuration
# ==========================================
if __name__ == "__main__":
    
    # Define custom environment run directories
    DATA_DIR = "./data/foldseek_result"
    AKK_DIR  = "./data/AKK_metadata"
    
    # Launch pipeline using modular arguments
    run_homolog_screening_pipeline(
        manual_selected_csv = os.path.join(DATA_DIR, "Akk_ATCC_mouse_mesangial_selected_mannual.csv"),
        mouse_id_mapping_tsv = os.path.join(DATA_DIR, "idmapping_ATCC_mouse_mesangial_seq.tsv"),
        human_id_mapping_tsv = os.path.join(DATA_DIR, "idmapping_ATCC_human_mesangial_seq.tsv"),
        akk_metadata_csv     = os.path.join(AKK_DIR,  "AKK.ATCC.protein.membrane.secretion.phi.MS_overlap.csv"),
        output_csv           = os.path.join(DATA_DIR, "Akk_mouse_human_screened_similarity_results.csv"),
        similarity_threshold = 80.0  # Optional screening criteria percentage limit (e.g., set to 80% matching)
    )
