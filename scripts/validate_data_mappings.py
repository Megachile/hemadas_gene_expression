#!/usr/bin/env python3
"""
validate_data_mappings.py
Created: 2025-03-05
Author: Claude

Purpose:
This script validates the consistency of mappings between different data files 
in the Hemadas gene expression project. It checks for missing mappings,
inconsistent IDs, and ensures traceability.

Inputs:
- data/consolidated_mads_proteins.csv
- data/consolidated_domain_boundaries.csv
- data/pair_to_protein_mapping.csv
- data/mads_proteins/mads_genes_final.csv
- data/mads_proteins/mads_genes_with_expression.csv
- data/mads_proteins/i_domain_boundaries.csv
- models/dimers/dimer_metadata_mapping.csv
- models/dimers/consolidated_dimer_results.csv
- models/dimers/standardized_files_manifest.csv
- models/tetramers/tetramer_analysis.csv
- models/tetramers/tetramer_status.csv
- models/tetramers/potential_tetramers.csv

Outputs:
- Console output with validation results
- Missing mappings report
"""

import os
import sys
import pandas as pd
import numpy as np

# Base paths
DATA_PATH = "/carc/scratch/users/akranz8174/hemadas_gene_expression/data"
MODELS_PATH = "/carc/scratch/users/akranz8174/hemadas_gene_expression/models"

# Define file paths
files = {
    "consolidated_mads_proteins": os.path.join(DATA_PATH, "consolidated_mads_proteins.csv"),
    "consolidated_domain_boundaries": os.path.join(DATA_PATH, "consolidated_domain_boundaries.csv"),
    "pair_to_protein_mapping": os.path.join(DATA_PATH, "pair_to_protein_mapping.csv"),
    "mads_genes_final": os.path.join(DATA_PATH, "mads_proteins", "mads_genes_final.csv"),
    "mads_genes_with_expression": os.path.join(DATA_PATH, "mads_proteins", "mads_genes_with_expression.csv"),
    "i_domain_boundaries": os.path.join(DATA_PATH, "mads_proteins", "i_domain_boundaries.csv"),
    "dimer_metadata_mapping": os.path.join(MODELS_PATH, "dimers", "dimer_metadata_mapping.csv"),
    "consolidated_dimer_results": os.path.join(MODELS_PATH, "dimers", "consolidated_dimer_results.csv"),
    "standardized_files_manifest": os.path.join(MODELS_PATH, "dimers", "standardized_files_manifest.csv"),
    "tetramer_analysis": os.path.join(MODELS_PATH, "tetramers", "tetramer_analysis.csv"),
    "tetramer_status": os.path.join(MODELS_PATH, "tetramers", "tetramer_status.csv"),
    "potential_tetramers": os.path.join(MODELS_PATH, "tetramers", "potential_tetramers.csv"),
    "gene_expression": os.path.join(DATA_PATH, "blueberry_cleaned_gene_expression.csv"),
    "integrated_domain_analysis": os.path.join(DATA_PATH, "integrated_domain_analysis.csv")
}

def check_file_exists(file_path, file_name):
    """Check if a file exists and print status"""
    exists = os.path.exists(file_path)
    status = "✓" if exists else "✗"
    print(f"{status} {file_name}: {file_path}")
    return exists

def load_dataframe(file_path):
    """Load data from CSV into a pandas DataFrame"""
    try:
        return pd.read_csv(file_path)
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None

def check_protein_ids_consistency():
    """Check consistency of protein IDs across different files"""
    print("\n== Checking protein ID consistency ==")
    
    # Load necessary files
    consolidated_proteins = load_dataframe(files["consolidated_mads_proteins"])
    domain_boundaries = load_dataframe(files["consolidated_domain_boundaries"])
    
    if consolidated_proteins is None or domain_boundaries is None:
        print("Cannot check protein IDs consistency - files missing")
        return
    
    # Get sets of protein IDs
    consolidated_protein_ids = set(consolidated_proteins["protein_id"])
    domain_boundary_protein_ids = set(domain_boundaries["protein_id"])
    
    # Check missing IDs
    missing_in_domains = consolidated_protein_ids - domain_boundary_protein_ids
    missing_in_consolidated = domain_boundary_protein_ids - consolidated_protein_ids
    
    print(f"Proteins in consolidated file: {len(consolidated_protein_ids)}")
    print(f"Proteins in domain boundaries: {len(domain_boundary_protein_ids)}")
    print(f"Proteins missing from domain boundaries: {len(missing_in_domains)}")
    print(f"Proteins missing from consolidated file: {len(missing_in_consolidated)}")
    
    if len(missing_in_domains) > 0:
        print("\nExample proteins missing from domain boundaries:")
        for pid in list(missing_in_domains)[:5]:
            print(f"  - {pid}")
    
    if len(missing_in_consolidated) > 0:
        print("\nExample proteins missing from consolidated file:")
        for pid in list(missing_in_consolidated)[:5]:
            print(f"  - {pid}")

def check_gene_ids_consistency():
    """Check consistency of gene IDs across different files"""
    print("\n== Checking gene ID consistency ==")
    
    # Load necessary files
    consolidated_proteins = load_dataframe(files["consolidated_mads_proteins"])
    mads_genes_final = load_dataframe(files["mads_genes_final"])
    mads_genes_with_expression = load_dataframe(files["mads_genes_with_expression"])
    
    if (consolidated_proteins is None or mads_genes_final is None 
            or mads_genes_with_expression is None):
        print("Cannot check gene IDs consistency - files missing")
        return
    
    # Get sets of gene IDs (filter out non-gene IDs like "Not in HMMER mapping")
    consolidated_gene_ids = set([g for g in consolidated_proteins["gene_id"] 
                                if isinstance(g, str) and not g.startswith("Not")])
    mads_genes_final_ids = set(mads_genes_final["Gene_ID"])
    mads_genes_with_expression_ids = set(mads_genes_with_expression["Gene_ID"])
    
    # Check missing IDs
    missing_in_final = consolidated_gene_ids - mads_genes_final_ids
    missing_in_expression = consolidated_gene_ids - mads_genes_with_expression_ids
    missing_in_consolidated_from_final = mads_genes_final_ids - consolidated_gene_ids
    
    print(f"Genes in consolidated file: {len(consolidated_gene_ids)}")
    print(f"Genes in final gene list: {len(mads_genes_final_ids)}")
    print(f"Genes in expression data: {len(mads_genes_with_expression_ids)}")
    print(f"Genes missing from final list: {len(missing_in_final)}")
    print(f"Genes missing from expression data: {len(missing_in_expression)}")
    print(f"Genes missing from consolidated but in final: {len(missing_in_consolidated_from_final)}")
    
    if len(missing_in_final) > 0:
        print("\nExample genes missing from final list:")
        for gid in list(missing_in_final)[:5]:
            print(f"  - {gid}")

def check_pair_consistency():
    """Check consistency of pair IDs across different files"""
    print("\n== Checking pair ID consistency ==")
    
    # Load necessary files
    pair_to_protein = load_dataframe(files["pair_to_protein_mapping"])
    dimer_metadata = load_dataframe(files["dimer_metadata_mapping"])
    dimer_results = load_dataframe(files["consolidated_dimer_results"])
    manifest = load_dataframe(files["standardized_files_manifest"])
    
    if (pair_to_protein is None or dimer_metadata is None 
            or dimer_results is None or manifest is None):
        print("Cannot check pair IDs consistency - files missing")
        return
    
    # Extract pair IDs
    pair_protein_ids = set(pair_to_protein["monomer_id"].str.split("_monomer").str[0].unique())
    dimer_metadata_ids = set(dimer_metadata["pair_id"].astype(str))
    dimer_result_ids = set(dimer_results["pair_id"].astype(str))
    manifest_ids = set(manifest["pair_id"].astype(str))
    
    # Check missing IDs
    missing_in_metadata = pair_protein_ids - dimer_metadata_ids
    missing_in_results = pair_protein_ids - dimer_result_ids
    missing_in_manifest = dimer_result_ids - manifest_ids
    
    print(f"Pairs in protein mapping: {len(pair_protein_ids)}")
    print(f"Pairs in dimer metadata: {len(dimer_metadata_ids)}")
    print(f"Pairs in dimer results: {len(dimer_result_ids)}")
    print(f"Pairs in manifest: {len(manifest_ids)}")
    print(f"Pairs missing from metadata: {len(missing_in_metadata)}")
    print(f"Pairs missing from results: {len(missing_in_results)}")
    print(f"Pairs missing from manifest: {len(missing_in_manifest)}")

def check_tetramer_consistency():
    """Check consistency of tetramer IDs across different files"""
    print("\n== Checking tetramer ID consistency ==")
    
    # Load necessary files
    tetramer_analysis = load_dataframe(files["tetramer_analysis"])
    tetramer_status = load_dataframe(files["tetramer_status"])
    potential_tetramers = load_dataframe(files["potential_tetramers"])
    
    if tetramer_analysis is None or tetramer_status is None or potential_tetramers is None:
        print("Cannot check tetramer IDs consistency - files missing")
        return
    
    # Extract tetramer IDs
    analysis_ids = set(tetramer_analysis["tetramer_name"])
    status_ids = set(tetramer_status["tetramer_name"])
    potential_ids = set(potential_tetramers["tetramer_id"])
    
    # Check missing IDs
    missing_in_status = analysis_ids - status_ids
    missing_in_analysis = status_ids - analysis_ids
    missing_in_potential = analysis_ids - potential_ids
    
    print(f"Tetramers in analysis: {len(analysis_ids)}")
    print(f"Tetramers in status: {len(status_ids)}")
    print(f"Tetramers in potential list: {len(potential_ids)}")
    print(f"Tetramers missing from status: {len(missing_in_status)}")
    print(f"Tetramers missing from analysis: {len(missing_in_analysis)}")
    print(f"Tetramers missing from potential list: {len(missing_in_potential)}")

def check_expression_gene_mapping():
    """Check if genes in expression data can be mapped to MADS proteins"""
    print("\n== Checking expression to MADS protein mapping ==")
    
    # Load necessary files
    gene_expression = load_dataframe(files["gene_expression"])
    mads_genes_with_expression = load_dataframe(files["mads_genes_with_expression"])
    consolidated_proteins = load_dataframe(files["consolidated_mads_proteins"])
    
    if (gene_expression is None or mads_genes_with_expression is None 
            or consolidated_proteins is None):
        print("Cannot check expression mapping - files missing")
        return
    
    # Get sets of gene IDs
    expression_gene_ids = set(gene_expression["Gene_ID"])
    mads_expression_gene_ids = set(mads_genes_with_expression["Gene_ID"])
    
    # MADS genes with no expression data
    mads_missing_expression = mads_expression_gene_ids - expression_gene_ids
    
    # Check which genes have non-zero expression in any condition
    mads_with_expression = mads_genes_with_expression[mads_genes_with_expression["Has_NonZero_Expression"] == True]
    mads_expressed_count = len(mads_with_expression)
    
    print(f"Total genes in expression data: {len(expression_gene_ids)}")
    print(f"MADS genes with expression data: {len(mads_expression_gene_ids)}")
    print(f"MADS genes missing from expression data: {len(mads_missing_expression)}")
    print(f"MADS genes with non-zero expression: {mads_expressed_count}")

def check_model_files_existence():
    """Check if the standardized PDB model files exist"""
    print("\n== Checking model file existence ==")
    
    manifest = load_dataframe(files["standardized_files_manifest"])
    
    if manifest is None:
        print("Cannot check model files - manifest missing")
        return
    
    # Sample a few files to check
    sample_size = min(10, len(manifest))
    samples = manifest.sample(sample_size)
    
    missing_count = 0
    for _, row in samples.iterrows():
        file_path = row["file_path"]
        exists = os.path.exists(file_path)
        if not exists:
            missing_count += 1
            print(f"Missing file: {file_path}")
    
    # Calculate estimated missing percentage
    missing_pct = (missing_count / sample_size) * 100
    print(f"Checked {sample_size} model files")
    print(f"Missing files: {missing_count} ({missing_pct:.1f}%)")

def print_schema_recommendations():
    """Print recommended schema for future data organization"""
    print("\n== Schema Recommendations ==")
    
    schema = """
    Key ID types and naming conventions:
    
    1. Gene IDs: Format 'maker-VaccDscaff{X}-{type}-gene-{Y}.{Z}-mRNA-{N}'
       - Primary identifier for genes in the Vaccinium corymbosum genome
       - Used to link genomic data with expression data
       
    2. Protein IDs: Format 'XP_{digits}.{version}'
       - Primary identifier for MADS-box proteins
       - Used to link structure models with functional annotations
       
    3. Pair IDs: Format 'pair{N}' (numeric IDs from 1-999+)
       - Primary identifier for protein dimers
       - Used to link dimer models with their constituent proteins
       
    4. Tetramer IDs: Format 'tetramer_{pair1}_{pair2}'
       - Compound identifier for tetramers formed from two dimers
       - Used to link tetramer models with their constituent dimers
    
    File linking relationships:
    
    Gene expression data → Gene IDs → Protein IDs → Pair IDs → Tetramer IDs
    
    Core linking files:
    
    1. consolidated_mads_proteins.csv: Maps Protein IDs to Gene IDs
    2. pair_to_protein_mapping.csv: Maps Pair IDs to Protein IDs
    3. dimer_metadata_mapping.csv: Maps Pair IDs to Gene IDs and Protein IDs
    4. mads_genes_with_expression.csv: Maps Gene IDs to expression data
    5. tetramer_status.csv: Maps Tetramer IDs to Pair IDs
    """
    
    print(schema)
    
def main():
    """Main function to run all validation checks"""
    print("Hemadas Gene Expression Data Validation")
    print("=======================================")
    
    # Check if required files exist
    print("\n== Checking file existence ==")
    all_files_exist = True
    for name, path in files.items():
        exists = check_file_exists(path, name)
        all_files_exist = all_files_exist and exists
    
    if not all_files_exist:
        print("\nSome required files are missing. Please check paths.")
        return
    
    # Run validation checks
    check_protein_ids_consistency()
    check_gene_ids_consistency()
    check_pair_consistency()
    check_tetramer_consistency()
    check_expression_gene_mapping()
    check_model_files_existence()
    
    # Print schema recommendations
    print_schema_recommendations()
    
if __name__ == "__main__":
    main()