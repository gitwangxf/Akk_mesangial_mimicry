# Molecular Mimicry Identification Pipeline
A bioinformatics pipeline designed to identify molecular mimicry pairs between _Akkermansia muciniphila_ (_Akk_) microbial proteins and human mesangial cell host proteins. By integrating immunology predictors with 3D structural alignment tools (Foldseek), this framework screens for bacterial epitopes capable of driving cross-reactive autoimmune responses.

## Environment & Dependencies

### 1. System Requirements & Command-Line Tools
Ensure the following binaries are installed and accessible in your system's `PATH`:
* **Foldseek**: For ultra-fast structural database searching. [Installation Guide](https://github.com/steineggerlab/foldseek)
* **PyMOL (Open-Source or Corporate)**: Used for structural manipulation, coordinate slicing, and RMSD verification.

### 2. Python Packages
The pipeline requires Python 3.8+. Install the required dependencies using `pip`:
```bash
pip install pandas numpy biopython matplotlib seaborn
```

## Step-by-Step Execution Guide

### Step 1: B-Cell and T-Cell Epitope Extraction

#### 1. Run Epitope Predictions:
Maps immunogenic regions within _Akkermansia muciniphila_ surface/secreted structures. It utilizes structural B-cell predictors (e.g., Discotope3) and sequence-based MHC-II binders (e.g., netMHCIIpan) to find vulnerable regions, then clips corresponding coordinates out of full-length PDB models.
```bash
bash scripts/predict_epitopes.sh
```

#### 2. Extract Epitope 3D Substructures:
Parse raw prediction outputs, locate consecutive high-scoring residue spans, and crop structural sub-segments.
```bash
python scripts/parse_extract_epitope_structures.py
```

### Step 2: Host Ortholog Sequence Identity Filtering
Matches mouse mesangial targets to their human orthologs, performs global alignments, and filters for pairs showing > 80% sequence identity.
```bash
python scripts/mouse_human_mesangial_map.py
```

### Step 3: Foldseek Local Structural Alignment
Aligns the structural _Akk_ epitope query models against the curated host mesangial database to identify similar spatial topologies.
```bash
bash scripts/foldseek_db_search.sh
```

### Step 4: Multi-Attribute Quality Filtering
Raw structural matches are subjected to strict geometric and biological filters to minimize false positives.
```bash
python scripts/get_Akk_mesangial_filterings.py
```
