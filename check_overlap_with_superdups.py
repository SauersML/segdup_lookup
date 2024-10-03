import pandas as pd
import gzip
import urllib.request
from tqdm import tqdm

# Function to download and unzip the genomicSuperDups file
def download_and_unzip(url, output_path):
    print(f"Downloading from {url}...")
    urllib.request.urlretrieve(url, output_path)
    print(f"Downloaded to {output_path}. Unzipping...")
    with gzip.open(output_path, 'rb') as f_in:
        with open(output_path.replace('.gz', ''), 'wb') as f_out:
            f_out.write(f_in.read())
    print(f"Unzipped to {output_path.replace('.gz', '')}.")

# Function to check for overlap
def is_overlap(start1, end1, start2, end2):
    return max(start1, start2) <= min(end1, end2)

# Function to calculate overlap size
def overlap_size(start1, end1, start2, end2):
    if is_overlap(start1, end1, start2, end2):
        return min(end1, end2) - max(start1, start2)
    return 0

# Main function to process overlap detection
def check_overlap_with_superdups(gene_file="gene_coordinates_with_group_name.tsv", 
                                 superdups_url="https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.txt.gz", 
                                 save_output=True):
    # Download and unzip genomicSuperDups file
    download_and_unzip(superdups_url, "genomicSuperDups.txt.gz")
    
    # Load the gene coordinates file
    print("Loading gene coordinates...")
    genes_df = pd.read_csv(gene_file, sep='\t')
    
    # Load the genomicSuperDups file
    print("Loading genomicSuperDups data...")
    dups_df = pd.read_csv("genomicSuperDups.txt", sep='\t', header=None,
                          names=["bin", "chrom", "chromStart", "chromEnd", "name", "score", "strand", 
                                 "otherChrom", "otherStart", "otherEnd", "otherSize", "uid", "posBasesHit", 
                                 "testResult", "verdict", "chits", "ccov", "alignfile", "alignL", "indelN", 
                                 "indelS", "alignB", "matchB", "mismatchB", "transitionsB", "transversionsB", 
                                 "fracMatch", "fracMatchIndel", "jcK", "k2K"])
    
    # Combine chrom and otherChrom into a single duplication regions DataFrame
    print("Combining primary and secondary duplication regions...")
    primary_dups = dups_df[['chrom', 'chromStart', 'chromEnd']].copy()
    primary_dups = primary_dups.rename(columns={'chrom': 'dup_chrom', 'chromStart': 'dup_start', 'chromEnd': 'dup_end'})
    
    secondary_dups = dups_df[['otherChrom', 'otherStart', 'otherEnd']].copy()
    secondary_dups = secondary_dups.rename(columns={'otherChrom': 'dup_chrom', 'otherStart': 'dup_start', 'otherEnd': 'dup_end'})
    
    # Concatenate primary and secondary duplications
    all_dups = pd.concat([primary_dups, secondary_dups], ignore_index=True)
    print(f"Total duplication regions to check: {len(all_dups)}")
    
    # Process and find overlaps
    overlap_results = []
    total_matches = 0
    print("Checking overlaps...")
    
    # tqdm for real-time progress tracking
    for _, gene_row in tqdm(genes_df.iterrows(), total=genes_df.shape[0], desc="Processing Genes"):
        gene_chrom = gene_row['Chromosome']
        gene_start = gene_row['Start']
        gene_end = gene_row['End']
        gene_name = gene_row['Ensembl gene ID']
        gene_group = gene_row['Group name']

        # Ensure chromosome labels match (e.g., 'chr1')
        gene_chrom_label = f"chr{gene_chrom}"
        
        # Filter duplications on the same chromosome as the gene
        dups_on_same_chrom = all_dups[all_dups['dup_chrom'] == gene_chrom_label]
        
        for _, dup_row in dups_on_same_chrom.iterrows():
            dup_start = dup_row['dup_start']
            dup_end = dup_row['dup_end']
            
            # Check overlap
            if is_overlap(gene_start, gene_end, dup_start, dup_end):
                overlap_len = overlap_size(gene_start, gene_end, dup_start, dup_end)
                overlap_region = f"{max(gene_start, dup_start)}-{min(gene_end, dup_end)}"
                gene_length = gene_end - gene_start
                overlap_percentage = (overlap_len / gene_length) * 100 if gene_length > 0 else 0
                overlap_type = "Duplication"
                match_type = "Duplication"
    
                total_matches += 1
                overlap_results.append({
                    "Gene": gene_name,
                    "Group Name": gene_group,
                    "Gene Chromosome": gene_chrom,
                    "Gene Start": gene_start,
                    "Gene End": gene_end,
                    "Duplication Chromosome": dup_row['dup_chrom'],
                    "Duplication Start": dup_start,
                    "Duplication End": dup_end,
                    "Overlap Region": overlap_region,
                    "Overlap Length": overlap_len,
                    "Gene Length": gene_length,
                    "Overlap Percentage": overlap_percentage,
                    "Overlap Type": overlap_type,
                    "Match Type": match_type
                })
    
    print(f"Total matches found: {total_matches}")
    
    overlap_df = pd.DataFrame(overlap_results)
    
    # Remove exact duplicate rows
    overlap_df = overlap_df.drop_duplicates()
    
    # Create GPCR output: filter Group Name containing 'G protein' and sort by Overlap Percentage
    gpcr_df = overlap_df[overlap_df['Group Name'].str.contains("G protein", case=False, na=False)]
    
    # Sort by Overlap Percentage in descending order
    gpcr_df = gpcr_df.sort_values(by="Overlap Percentage", ascending=False)
                                        
    if save_output:
        # Save the main output file
        overlap_df.to_csv("overlap_results.csv", index=False)
        print("Results saved to overlap_results.csv")
                 
        # Save the GPCR output file
        gpcr_df.to_csv("GPCR_output.csv", index=False)
        print("GPCR results saved to GPCR_output.csv")
    else:
        print(overlap_df)
        print(gpcr_df)
    
    return overlap_results

if __name__ == "__main__":
    check_overlap_with_superdups()
