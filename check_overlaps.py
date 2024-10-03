import csv

def trim_chr(chrom):
    return chrom.replace('chr', '')

def load_genes(tsv_file):
    genes = []
    g_protein_genes = []
    with open(tsv_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            gene = {
                'id': parts[0],
                'chrom': trim_chr(parts[1]),
                'start': int(parts[2]),
                'end': int(parts[3]),
                'group': parts[4] if len(parts) > 4 else "",
                'is_g_protein': "G protein" in parts[4] if len(parts) > 4 else False
            }
            genes.append(gene)
            if gene['is_g_protein']:
                g_protein_genes.append(gene)
    return genes, g_protein_genes

def load_duplication_regions(bed_file):
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            regions.append({
                'chrom': trim_chr(parts[0]),
                'start': int(parts[1]),
                'end': int(parts[2])
            })
    return regions

def check_overlap(gene, region):
    if gene['chrom'] == region['chrom']:
        overlap_start = max(gene['start'], region['start'])
        overlap_end = min(gene['end'], region['end'])
        if overlap_start <= overlap_end:
            return {
                'chrom': gene['chrom'],
                'start': overlap_start,
                'end': overlap_end
            }
    return None

def save_results_to_file(all_overlaps, output_file):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        # Write header
        header = [
            "Gene_ID", "Is_G_Protein", "Gene_Group",
            "Gene_Chr", "Gene_Start", "Gene_End",
            "Dup_Chr", "Dup_Start", "Dup_End",
            "Overlap_Chr", "Overlap_Start", "Overlap_End"
        ]
        writer.writerow(header)
        
        # Sort results to put G protein genes at the top
        sorted_results = sorted(all_overlaps, 
                               key=lambda x: (not x['gene']['is_g_protein'], x['gene']['id']))
        
        # Write data
        for match in sorted_results:
            gene = match['gene']
            dup = match['dup_region']
            overlap = match['overlap']
            row = [
                gene['id'],
                'Yes' if gene['is_g_protein'] else 'No',
                gene['group'],
                gene['chrom'],
                gene['start'],
                gene['end'],
                dup['chrom'],
                dup['start'],
                dup['end'],
                overlap['chrom'],
                overlap['start'],
                overlap['end']
            ]
            writer.writerow(row)

# Main execution
gene_file = 'gene_coordinates_with_group_name.tsv'
duplication_file = 'wgac_filtered.no_alt.bed'
output_file = 'gene_duplication_overlaps.csv'

print("Loading genes...")
genes, g_protein_genes = load_genes(gene_file)
print(f"Loaded {len(genes)} total genes")
print(f"Found {len(g_protein_genes)} G protein genes")

print("Loading duplication regions...")
duplication_regions = load_duplication_regions(duplication_file)
print(f"Loaded {len(duplication_regions)} duplication regions")

# Process all genes
all_overlaps = []
genes_with_overlaps = set()
g_protein_genes_with_overlaps = set()
progress_interval = max(1, len(genes) // 20)

print("Checking overlaps...")
for i, gene in enumerate(genes):
    if (i + 1) % progress_interval == 0:
        print(f"Processed {i+1}/{len(genes)} genes ({((i+1)/len(genes))*100:.1f}%)")
    
    found_overlap = False
    for region in duplication_regions:
        overlap = check_overlap(gene, region)
        if overlap:
            all_overlaps.append({
                'gene': gene,
                'dup_region': region,
                'overlap': overlap
            })
            genes_with_overlaps.add(gene['id'])
            if gene['is_g_protein']:
                g_protein_genes_with_overlaps.add(gene['id'])
            found_overlap = True

print(f"\nResults:")
print(f"Total genes: {len(genes)}")
print(f"Total unique genes overlapping with duplications: {len(genes_with_overlaps)}")
print(f"G protein genes overlapping with duplications: {len(g_protein_genes_with_overlaps)}")
print(f"Total number of overlaps found: {len(all_overlaps)}")
print(f"Average overlaps per overlapping gene: {len(all_overlaps) / len(genes_with_overlaps):.2f}")

# Save results to file
print(f"Saving results to {output_file}...")
save_results_to_file(all_overlaps, output_file)
print("Done!")
