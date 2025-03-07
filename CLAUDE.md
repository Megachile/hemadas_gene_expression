# Hemadas nubilipennis Gene Expression Project

## Project Overview
This repository contains consolidated data for studying gene expression changes induced by the gall wasp Hemadas nubilipennis. The project focuses on analyzing protein interactions of MADS-box transcription factors (dimers and tetramers) and their relationship to gene expression changes in blueberry (Vaccinium corymbosum).

## Directory Structure

```
/carc/scratch/users/akranz8174/hemadas_gene_expression/
├── data/                 # Core data files and mapping information
│   ├── blueberry_cleaned_gene_expression.csv  # Gene expression data
│   ├── consolidated_domain_boundaries.csv     # Domain boundary information for proteins
│   ├── consolidated_mads_proteins.csv         # Master protein list with gene IDs
│   ├── genome/                                # Reference genome files
│   │   ├── V_corymbosum_genome_v1.0.fasta     # Blueberry genome sequence
│   │   ├── V_corymbosum_genome_v1.0.fasta.fai # Genome index
│   │   ├── V_corymbosum_mRNA.gff              # mRNA annotations
│   │   ├── V_corymbosum_v1.0_geneModels.gff   # Gene model annotations
│   │   └── V_corymbosum_v1.0_proteins_.fasta  # Protein sequences
│   ├── integrated_domain_analysis.csv         # Domain-specific structure analysis
│   ├── mads_proteins/                         # MADS-box protein information
│   │   ├── i_domain_boundaries.csv            # I-domain boundaries
│   │   ├── mads_genes_final.csv               # Final MADS gene list
│   │   └── mads_genes_with_expression.csv     # MADS genes with expression data
│   └── pair_to_protein_mapping.csv            # Maps pair IDs to protein IDs
├── models/               # Protein structure models
│   ├── dimers/                                # MADS-box dimer models
│   │   ├── consolidated_dimer_results.csv     # Comprehensive dimer results
│   │   ├── dimer_metadata_mapping.csv         # Maps pairs to proteins/genes
│   │   ├── dimer_structure_metrics.csv        # Structural metrics for each model
│   │   ├── standardized_files/                # Standardized PDB model files
│   │   └── standardized_files_manifest.csv    # Maps pair IDs to model files
│   └── tetramers/                             # MADS-box tetramer models
│       ├── potential_tetramers.csv            # Potential tetramer combinations
│       ├── tetramer_analysis.csv              # Analysis of tetramer structures
│       └── tetramer_status.csv                # Status of tetramer modeling
├── docs/                 # Documentation and notes
├── outputs/              # Analysis results and visualizations
└── scripts/              # Analysis and processing scripts
```

## Data Organization and ID Mapping Structure

### Key ID Types and Naming Conventions

1. **Gene IDs**: Format `maker-VaccDscaff{X}-{type}-gene-{Y}.{Z}-mRNA-{N}`
   - Primary identifier for genes in the Vaccinium corymbosum genome
   - Example: `maker-VaccDscaff3-augustus-gene-198.30-mRNA-1`
   - Used to link genomic data with expression data

2. **Protein IDs**: Format `XP_{digits}.{version}`
   - Primary identifier for MADS-box proteins
   - Example: `XP_021970064.1`
   - Used to link structure models with functional annotations

3. **Pair IDs**: Format `pair{N}` (numeric IDs)
   - Primary identifier for protein dimers
   - Example: `pair714`
   - Used to link dimer models with their constituent proteins

4. **Tetramer IDs**: Format `tetramer_{pair1}_{pair2}`
   - Compound identifier for tetramers formed from two dimers
   - Example: `tetramer_796_1313`
   - Used to link tetramer models with their constituent dimers

### Mapping Relationships

The data is organized in a hierarchical structure with the following relationships:

```
Gene expression data → Gene IDs → Protein IDs → Pair IDs → Tetramer IDs
```

### Core Mapping Files

1. **consolidated_mads_proteins.csv**: Maps Protein IDs to Gene IDs
   - Key columns: `protein_id`, `gene_id`

2. **pair_to_protein_mapping.csv**: Maps Pair IDs to Protein IDs
   - Key columns: `monomer_id` (contains pair ID), `protein_id`

3. **dimer_metadata_mapping.csv**: Maps Pair IDs to Gene IDs and Protein IDs
   - Key columns: `pair_id`, `gene1_id`, `gene2_id`, `protein1_id`, `protein2_id`

4. **mads_genes_with_expression.csv**: Maps Gene IDs to expression data
   - Key columns: `Gene_ID`, expression columns (tissue-specific FPKM values)

5. **tetramer_status.csv**: Maps Tetramer IDs to Pair IDs
   - Key columns: `tetramer_name`, `pair1_id`, `pair2_id`

### Current Data Status

Based on validation checks, there are some inconsistencies in the datasets:

1. **Protein Coverage**: 
   - Some proteins appear in domain boundaries but are missing from the consolidated proteins list
   - Recommended fix: Update consolidated_mads_proteins.csv to include all proteins

2. **Gene Coverage**:
   - The consolidated proteins file contains fewer genes than the final gene list
   - Good coverage between gene lists and expression data
   - Recommended fix: Ensure all genes in mads_genes_final.csv are represented in consolidated_mads_proteins.csv

3. **Tetramer Models**:
   - Few tetramers have complete analysis results (2 of 9 in status list)
   - Most potential tetramers have not been modeled
   - Recommended fix: Add metadata about which tetramers are highest priority for modeling

## Key Data Files

### Genome and Expression Data
- `data/genome/V_corymbosum_genome_v1.0.fasta`: Blueberry (Vaccinium corymbosum) reference genome
- `data/genome/V_corymbosum_v1.0_geneModels.gff`: Gene model annotations
- `data/blueberry_cleaned_gene_expression.csv`: Gene expression data for analysis

### MADS Protein Information
- `data/mads_proteins/mads_genes_final.csv`: Curated list of MADS-box genes
- `data/mads_proteins/mads_genes_with_expression.csv`: MADS genes with expression data
- `data/mads_proteins/i_domain_boundaries.csv`: I-domain boundaries in MADS proteins
- `data/consolidated_mads_proteins.csv`: Master list of MADS-box proteins with gene IDs and annotations
- `data/consolidated_domain_boundaries.csv`: Defines the domain boundaries for each protein
- `data/pair_to_protein_mapping.csv`: Maps pair IDs to their constituent protein IDs
- `data/integrated_domain_analysis.csv`: Domain-specific analysis of protein structures

### Protein Structure Models
- `models/dimers/consolidated_dimer_results.csv`: Comprehensive analysis of all dimer models
- `models/dimers/dimer_metadata_mapping.csv`: Maps pair IDs to protein/gene information
- `models/dimers/dimer_structure_metrics.csv`: Quality metrics for each dimer model
- `models/dimers/standardized_files_manifest.csv`: Catalog of all standardized model files

### Tetramer Models
- `models/tetramers/tetramer_analysis.csv`: Analysis of completed tetramer structures
- `models/tetramers/tetramer_status.csv`: Status information on all tetramers
- `models/tetramers/potential_tetramers.csv`: List of potential tetramer combinations

## Best Practices for Future Work

### File Organization
1. **Directories**: 
   - `/data`: Raw and processed data files
   - `/scripts`: All analysis scripts
   - `/outputs`: Results and visualizations
   - `/models`: Structure models and related files
   - `/docs`: Documentation and notes

2. **File Naming Conventions**:
   - For data files: `[entity_type]_[detail]_[version].[extension]`
   - For output files: `[descriptive_name]_[script_name]_[YYYYMMDD].[extension]`
   - Always use underscores, not spaces or hyphens
   - Include version numbers or dates in filenames

### Script Organization and Documentation
1. All scripts should be stored in the `/scripts` directory
2. Every script should include a detailed header with:
   ```python
   """
   script_name.py
   Created: YYYY-MM-DD
   Author: [Author Name]
   
   Purpose:
   [Brief description of what the script does]
   
   Input Files:
   - /path/to/input1.csv: [Description of input1]
   - /path/to/input2.csv: [Description of input2]
   
   Output Files:
   - /path/to/output1.csv: [Description of output1]
   - /path/to/output2.png: [Description of output2]
   
   Based on previous version:
   - previous_script.py (YYYY-MM-DD)
   
   Changes from previous version:
   - [List major changes]
   """
   ```

3. **Code Documentation**:
   - Standardize column names across all files (use snake_case)
   - Add validation checks to ensure data integrity
   - Include comments explaining complex operations and data transformations

### Data Integrity Guidelines
1. **ID Consistency**: Always use consistent ID formats across all files:
   - Gene IDs should always be full `maker-VaccDscaff` format
   - Protein IDs should always be full `XP_` format with version
   - Pair IDs should always include the `pair` prefix
   - Tetramer IDs should always follow `tetramer_X_Y` format

2. **Required Metadata**: Each analysis file should include:
   - Creation date
   - Script that generated it
   - Version/parameters information
   - Brief description of the data

3. **Validation**: Run the `validate_data_mappings.py` script after any major data update

### Version Control and GitHub Workflow
This project now uses GitHub for version control. The repository is available at:
https://github.com/Megachile/hemadas_gene_expression

#### Local and Remote Repository
1. The local repository is at: `/carc/scratch/users/akranz8174/hemadas_gene_expression/`
2. The remote repository is: `https://github.com/Megachile/hemadas_gene_expression`
3. Authentication uses a Personal Access Token (PAT)

#### Git and GitHub Workflow
1. Commit message format: `[script/data]: [brief description of changes]`
2. Use branches for experimental features
3. Include updated CLAUDE.md in commits when project organization changes
4. **IMPORTANT**: All analysis outputs should be committed and pushed to GitHub
5. Large data files are excluded via .gitignore and remain only in the local repository

#### Commands for Common Tasks
```bash
# Check status
git status

# Add files
git add filename.ext

# Commit changes
git commit -m "[category]: Brief description of changes"

# Push to GitHub
git push origin main

# Pull latest changes
git pull origin main
```

## Accessing Legacy Files
If needed, original files can be found in the following locations:
- Original MADS-box dimer files: `/carc/scratch/users/akranz8174/mads_box_dimers/`
- Original chromatin repression analysis: `/carc/scratch/users/akranz8174/chromatin_repression/`

## Output File Management

### GitHub Storage Guidelines
The following files should be committed and pushed to GitHub:
1. **Analysis Scripts**: All scripts in the `scripts/` directory
2. **Documentation**: All files in the `docs/` directory
3. **Analysis Results**: All files in the `outputs/` directory including:
   - CSV data tables
   - JSON analysis results
   - Summary text files
   - Figures and plots (PNG, JPG, PDF)
   - Statistical analysis outputs
4. **Configuration Files**: Any parameter or configuration files
5. **Metadata**: All data mapping and annotation files (not raw data)

### Local-Only Storage
The following files should NOT be committed to GitHub and should remain local only:
1. **Genome Files**: All genome FASTA and annotation files
2. **Structure Models**: All PDB files for protein structures
3. **Raw Sequencing Data**: Any raw sequence data
4. **Large Datasets**: Any file larger than 50MB

## Next Steps
- Complete missing protein and gene mappings
- Integrate gene expression analysis with protein structure information
- Analyze how MADS-box dimer formation correlates with gene expression changes
- Identify key regulatory interactions in gall formation
- Upload initial analysis results to GitHub repository