## Quick start

Run the program like this:
```
mkdir segdup_lookup_folder && cd segdup_lookup_folder && git clone https://github.com/ScottSauers/segdup_lookup.git . && python3 main.py
```

# Segmental Duplication Lookup Tool

A Python tool to identify overlaps between gene coordinates and genomic segmental duplications.

## Overview

This tool performs two main functions:
1. Fetches gene coordinates from Ensembl for a predefined list of genes
2. Checks for overlaps between these genes and known segmental duplications

## Requirements

- Python 3
- Required packages: `pandas`, `requests`, `tqdm`

```pip3 install pandas tqdm requests```

## Installation

```bash
mkdir segdup_lookup_folder
cd segdup_lookup_folder
git clone https://github.com/ScottSauers/segdup_lookup.git .
```

## Usage

Run the main script:

```bash
python3 main.py
```

This will:
1. Download gene coordinates from Ensembl
2. Fetch segmental duplication data
3. Give information about overlaps between genes and duplications

## Output

- `overlap_results.csv`: Overlap information
- `gene_coordinates_with_group_name.tsv`: Gene coordinates
- `gene_errors.log`: Any errors encountered during execution

## Data sources
- Segmental duplications: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.txt.gz
- Gene list: https://www.genenames.org/cgi-bin/genegroup/download?id=2054&type=branch
- Gene coordinates: https://rest.ensembl.org/lookup/id/{ensembl_id}?content-type=application/json

## Note

This script is untested and may contain inaccuracies or bugs. If an overlap type is primary but also other, it will be classified as primary. The output file may have duplicates.
