# Segmental Duplication Lookup Tool

A Python tool to identify overlaps between gene coordinates and genomic segmental duplications from the UCSC genomicSuperDups database ("Schema for Segmental Dups - Duplications of >1000 Bases of Non-RepeatMasked Sequence," data provided by Ginger Cheng, Xinwei She, Archana Raja, Tin Louie and Evan Eichler).

## Quick Start

Run the program with the following command:
```bash
mkdir segdup_lookup_folder && cd segdup_lookup_folder && git clone https://github.com/ScottSauers/segdup_lookup.git . && python3 main.py
```

## Overview

This tool performs two main functions:
1. Fetches gene coordinates from Ensembl for a predefined list of genes.
2. Checks for overlaps between these genes and known **segmental duplications** based on the UCSC `genomicSuperDups` database.

## Requirements

- **Python 3** is required.
- The following Python packages are needed: `pandas`, `requests`, and `tqdm`.
  
Install the required packages using:
```bash
pip3 install pandas tqdm requests
```

## Installation

To install the tool, use the following steps:
```bash
mkdir segdup_lookup_folder
cd segdup_lookup_folder
git clone https://github.com/ScottSauers/segdup_lookup.git .
```

## Usage

Once installed, run the main script with the following command:
```bash
python3 main.py
```

The script:
1. **Downloads gene coordinates** from Ensembl for predefined genes.
2. **Fetches segmental duplication data** from UCSC.
3. **Identifies overlaps** between the gene coordinates and the duplication regions.

## Output Files

- `overlap_results.csv`: Contains information about overlaps between genes and duplication regions.
- `GPCR_output.csv`: Overlap information specifically for GPCR-related genes.
- `gene_coordinates_with_group_name.tsv`: The list of gene coordinates used in the analysis.
- `gene_errors.log`: Logs any errors encountered during the process.

## Explanation of Overlap Matches

### Direct Match:
A **direct match** occurs when a gene physically overlaps a duplication region **on the same chromosome**.

### Indirect Match:
A **non-direct (indirect) match** occurs if:
- The gene overlaps a duplication on **a different chromosome**
- The gene overlaps a duplication on the **same chromosome** but in a **different, non-overlapping region**.

For **indirect matches**, the reported overlap percentage may refer to how much of the gene would overlap with the duplicated region, but it does not indicate a physical overlap on the gene's region.

## Data Sources

- **Segmental duplications**: [UCSC GenomicSuperDups](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.txt.gz)
- **Gene list**: [HGNC gene groups](https://www.genenames.org/cgi-bin/genegroup/download?id=2054&type=branch)
- **Gene coordinates**: [Ensembl REST API](https://rest.ensembl.org/lookup/id/{ensembl_id}?content-type=application/json)

## Notes
- This tool is untested, still under development, and may contain bugs.
