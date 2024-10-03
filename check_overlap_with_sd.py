import pandas as pd
from tqdm import tqdm

def check_overlap_with_sd(sd_file, save_output=True):
    """
    Identifies overlaps between gene coordinates and segmental duplications (SDs).

    Parameters:
    - gene_file (str): Path to the gene coordinates file (TSV with columns including 'Chromosome', 'Start', 'End', 'Ensembl gene ID', 'Group name').
    - sd_file (str): Path to the SD file (BED-like with columns: chrom, start, end, sd_identity).
    - save_output (bool): If True, saves the overlap results to CSV files.

    Returns:
    - overlap_df (pd.DataFrame): DataFrame containing all overlap information.
    - gpcr_df (pd.DataFrame): DataFrame containing overlaps for GPCR-related genes.
    """

    gene_file = "gene_coordinates_with_group_name.tsv"

    # Load gene coordinates
    print("Loading gene coordinates...")
    try:
        genes_df = pd.read_csv(gene_file, sep='\t')
    except Exception as e:
        print(f"Error loading gene file: {e}")
        return

    # Validate necessary columns in gene file
    required_gene_cols = {'Chromosome', 'Start', 'End', 'Ensembl gene ID', 'Group name'}
    if not required_gene_cols.issubset(genes_df.columns):
        print(f"Gene file is missing one or more required columns: {required_gene_cols}")
        return

    # Load SD regions
    print("Loading segmental duplications data...")
    try:
        sd_df = pd.read_csv(sd_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'sd_identity'])
    except Exception as e:
        print(f"Error loading SD file: {e}")
        return

    # Initialize list to store overlap results
    overlap_results = []
    total_matches = 0

    print("Checking overlaps between genes and segmental duplications...")

    # Iterate over genes with progress bar
    for _, gene in tqdm(genes_df.iterrows(), total=genes_df.shape[0], desc="Processing Genes"):
        gene_chrom = f"chr{gene['Chromosome']}"
        gene_start = gene['Start']
        gene_end = gene['End']
        gene_name = gene['Ensembl gene ID']
        gene_group = gene['Group name']

        # Filter SDs on the same chromosome
        sd_on_chrom = sd_df[sd_df['chrom'] == gene_chrom]

        # Iterate over SDs on the chromosome and check for overlap
        for _, sd in sd_on_chrom.iterrows():
            sd_start = sd['start']
            sd_end = sd['end']
            sd_identity = sd['sd_identity']

            # Check if gene overlaps with SD
            if max(gene_start, sd_start) <= min(gene_end, sd_end):
                overlap_length = min(gene_end, sd_end) - max(gene_start, sd_start)
                overlap_region = f"{max(gene_start, sd_start)}-{min(gene_end, sd_end)}"
                gene_length = gene_end - gene_start
                overlap_percentage = (overlap_length / gene_length) * 100 if gene_length > 0 else 0

                overlap_results.append({
                    "Gene": gene_name,
                    "Group Name": gene_group,
                    "Gene Chromosome": gene['Chromosome'],
                    "Gene Start": gene_start,
                    "Gene End": gene_end,
                    "Duplication Chromosome": sd['chrom'],
                    "Duplication Start": sd_start,
                    "Duplication End": sd_end,
                    "Overlap Region": overlap_region,
                    "Overlap Length": overlap_length,
                    "Gene Length": gene_length,
                    "Overlap Percentage": overlap_percentage,
                    "SD Identity": sd_identity
                })

                total_matches += 1

    print(f"Total overlaps found: {total_matches}")

    # Create DataFrame from overlap results
    overlap_df = pd.DataFrame(overlap_results)

    # Remove exact duplicate rows
    overlap_df = overlap_df.drop_duplicates()

    # Create GPCR output: filter Group Name containing 'G protein' and sort by Overlap Percentage
    gpcr_df = overlap_df[overlap_df['Group Name'].str.contains("G protein", case=False, na=False)]

    # Sort by Overlap Percentage in descending order
    gpcr_df = gpcr_df.sort_values(by="Overlap Percentage", ascending=False)

    # Save output files if required
    if save_output:
        try:
            overlap_df.to_csv("overlap_results.csv", index=False)
            print("Saved overlap results to 'overlap_results.csv'")
        except Exception as e:
            print(f"Error saving overlap_results.csv: {e}")

        try:
            gpcr_df.to_csv("GPCR_output.csv", index=False)
            print("Saved GPCR overlaps to 'GPCR_output.csv'")
        except Exception as e:
            print(f"Error saving GPCR_output.csv: {e}")
    else:
        print("All Overlaps:")
        print(overlap_df)
        print("\nGPCR Overlaps:")
        print(gpcr_df)

    return overlap_df, gpcr_df
