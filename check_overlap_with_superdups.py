import pandas as pd
import gzip
import urllib.request
from tqdm import tqdm

# Function to download and unzip the genomicSuperDups file
def download_and_unzip(url, output_path):
    urllib.request.urlretrieve(url, output_path)
    with gzip.open(output_path, 'rb') as f_in:
        with open(output_path.replace('.gz', ''), 'wb') as f_out:
            f_out.write(f_in.read())

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
    print("Downloading and unzipping genomicSuperDups file...")
    download_and_unzip(superdups_url, "genomicSuperDups.txt.gz")
    print("Download complete.")

    # Load the gene coordinates file
    print("Loading gene coordinates file...")
    genes_df = pd.read_csv(gene_file, sep='\t')
    
    # Load the genomicSuperDups file
    print("Loading genomicSuperDups file...")
    dups_df = pd.read_csv("genomicSuperDups.txt", sep='\t', header=None,
                          names=["bin", "chrom", "chromStart", "chromEnd", "name", "score", "strand", 
                                 "otherChrom", "otherStart", "otherEnd", "otherSize", "uid", "posBasesHit", 
                                 "testResult", "verdict", "chits", "ccov", "alignfile", "alignL", "indelN", 
                                 "indelS", "alignB", "matchB", "mismatchB", "transitionsB", "transversionsB", 
                                 "fracMatch", "fracMatchIndel", "jcK", "k2K"])

    # Process and find overlaps
    overlap_results = []
    total_matches = 0
    print("Checking overlaps between gene coordinates and genomicSuperDups...")

    # Use tqdm for real-time progress tracking
    for _, gene_row in tqdm(genes_df.iterrows(), total=genes_df.shape[0], desc="Gene Progress"):
        gene_chrom = gene_row['Chromosome']
        gene_start = gene_row['Start']
        gene_end = gene_row['End']
        gene_name = gene_row['Ensembl gene ID']
        gene_group = gene_row['Group name']

        # Filter superdups for matching chromosome
        dups_on_same_chrom = dups_df[dups_df['chrom'] == f"chr{gene_chrom}"]
        
        for _, dup_row in dups_on_same_chrom.iterrows():
            dup_start = dup_row['chromStart']
            dup_end = dup_row['chromEnd']
            other_chrom = dup_row['otherChrom']
            other_start = dup_row['otherStart']
            other_end = dup_row['otherEnd']
            
            # Check overlap with primary region
            matched = False
            overlap_type = ""
            if is_overlap(gene_start, gene_end, dup_start, dup_end):
                overlap_len = overlap_size(gene_start, gene_end, dup_start, dup_end)
                overlap_region = f"{max(gene_start, dup_start)}-{min(gene_end, dup_end)}"
                gene_length = gene_end - gene_start
                overlap_percentage = (overlap_len / gene_length) * 100
                overlap_type = "Primary"
                matched = True
  
            # Check overlap with otherChrom region
            elif is_overlap(gene_start, gene_end, other_start, other_end):
                overlap_len = overlap_size(gene_start, gene_end, other_start, other_end)
                overlap_region = f"{max(gene_start, other_start)}-{min(gene_end, other_end)}"
                gene_length = gene_end - gene_start
                overlap_percentage = (overlap_len / gene_length) * 100
                overlap_type = "Other"
                matched = True

            if matched:
                total_matches += 1
                overlap_results.append({
                    "Gene": gene_name,
                    "Group Name": gene_group,
                    "Chromosome of Gene": gene_chrom,
                    "Overlap Region": overlap_region,
                    "Overlap Length": overlap_len,
                    "Gene Length": gene_length,
                    "Overlap Percentage": overlap_percentage,
                    "Overlap Type": overlap_type
                })

                # Print progress info for every 100 matches found
                if total_matches % 100 == 0:
                    print(f"Matches found: {total_matches}")
    
    print(f"Total matches found: {total_matches}")

    if save_output:
        overlap_df = pd.DataFrame(overlap_results)
        overlap_df.to_csv("overlap_results.csv", index=False)
        print("Results saved to overlap_results.csv")
    else:
        print(pd.DataFrame(overlap_results))
    
    return overlap_results

if __name__ == "__main__":
    check_overlap_with_superdups()
