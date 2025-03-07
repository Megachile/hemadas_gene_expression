# Hemadas nubilipennis Gene Expression Project

## Project Overview
This repository contains data and analysis tools for studying gene expression changes induced by the gall wasp *Hemadas nubilipennis*. The project focuses on analyzing protein interactions of MADS-box transcription factors (dimers and tetramers) and their relationship to gene expression changes in blueberry (*Vaccinium corymbosum*).

## Repository Structure

```
├── data/                 # Core data files and mapping information
│   ├── gene_expression/  # Gene expression data and annotations
│   ├── genome/           # Reference genome files (not in repo - documented only)
│   ├── mads_proteins/    # MADS-box protein information
│   └── mapping/          # ID mapping files between datasets
├── models/               # Protein structure model metadata
│   ├── dimers/           # MADS-box dimer model information
│   └── tetramers/        # MADS-box tetramer model information
├── docs/                 # Documentation and notes
├── outputs/              # Analysis results and visualizations
└── scripts/              # Analysis and processing scripts
```

## Data Overview

This repository contains metadata and analysis code only. Large data files (genome files, protein structure models) are stored separately and referenced in the documentation.

### Key Data Types

1. **Gene Expression Data**: Tissue-specific expression profiles from gall-forming and non-gall tissues
2. **MADS-box Protein Data**: Annotations and domain information for MADS transcription factors
3. **Protein Structure Models**: Metadata for dimer and tetramer protein interaction models

## Getting Started

### Dependencies
- Python 3.8+
- pandas
- numpy
- matplotlib
- biopython

### Installation
```bash
git clone https://github.com/Megachile/hemadas_gene_expression.git
cd hemadas_gene_expression
```

## Data Validation

Run the validation script to check data consistency:
```bash
python scripts/validate_data_mappings.py
```

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

Adam Kranz - adamjameskranz@gmail.com

Project Link: [https://github.com/Megachile/hemadas_gene_expression](https://github.com/Megachile/hemadas_gene_expression)