import requests
import csv
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

def fetch_gene_coordinates():
    # Download the gene data file
    url = "https://www.genenames.org/cgi-bin/genegroup/download?id=2054&type=branch"
    response = requests.get(url)
    file_content = response.content.decode('utf-8')

    # Save the file content to a local file
    file_path = "gene_data.tsv"
    with open(file_path, 'w') as f:
        f.write(file_content)

    print("Gene data file downloaded and saved.")

    # Load already existing coordinates from gene_coordinates_with_group_name.tsv
    existing_coordinates = {}
    output_file = "gene_coordinates_with_group_name.tsv"
    if os.path.exists(output_file):
        with open(output_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                existing_coordinates[row['Ensembl gene ID']] = row['Group name']  # Store the group name along with the gene ID
        print(f"Loaded {len(existing_coordinates)} existing records from {output_file}")

    # Extract Ensembl gene IDs and Group name from the file
    gene_data = []
    with open(file_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ensembl_gene_id = row['Ensembl gene ID']
            group_name = row['Group name']
            if ensembl_gene_id:  # Gene ID is not empty
                gene_data.append({
                    'ensembl_gene_id': ensembl_gene_id,
                    'group_name': group_name
                })

    print(f"Extracted {len(gene_data)} gene records.")

    # Exponential backoff
    def get_coordinates(ensembl_id, retries=4, backoff_factor=2):
        if not ensembl_id:  # Gene ID is valid
            print("Skipping empty gene ID.")
            return None, None, None
        
        url = f"https://rest.ensembl.org/lookup/id/{ensembl_id}?content-type=application/json"
        backoff = 4  # Initial backoff time in seconds
        for attempt in range(retries):
            try:
                response = requests.get(url, timeout=10)
                if response.status_code == 200:
                    data = response.json()
                    chrom = data.get('seq_region_name')
                    start = data.get('start')
                    end = data.get('end')
                    if chrom and start and end:
                        return chrom, start, end
                elif response.status_code == 429:  # Handle rate limiting
                    print(f"Rate limit hit for {ensembl_id}, sleeping for {backoff} seconds")
                    time.sleep(backoff)
                    backoff *= backoff_factor  # Exponential backoff
            except requests.exceptions.RequestException as e:
                print(f"Error fetching data for {ensembl_id}: {str(e)}")
            print(f"Retry {attempt + 1}/{retries} for {ensembl_id}")
        return None, None, None  # Return None if all retries fail

    # Fetch coordinates concurrently
    def fetch_coordinates_concurrently(gene_data):
        total = len(gene_data)
        coordinates = {}
        errors = {}

        print("Fetching coordinates...")

        with ThreadPoolExecutor(max_workers=50) as executor:  # Using 50 threads for faster fetching
            futures = {
                executor.submit(get_coordinates, gene['ensembl_gene_id']): gene
                for gene in gene_data if gene['ensembl_gene_id'] not in existing_coordinates
            }
            for i, future in enumerate(as_completed(futures)):
                gene = futures[future]
                try:
                    ensembl_id = gene['ensembl_gene_id']
                    chrom, start, end = future.result()
                    if chrom:
                        coordinates[ensembl_id] = {
                            'chromosome': chrom,
                            'start': start,
                            'end': end,
                            'group_name': gene['group_name']
                        }
                        print(f"Found: {ensembl_id} -> {chrom}:{start}-{end}")
                    else:
                        print(f"No data for {ensembl_id}")
                        errors[ensembl_id] = "No data returned"
                except Exception as e:
                    print(f"Error processing {ensembl_id}: {str(e)}")
                    errors[ensembl_id] = str(e)

                # Print progress every 10 genes
                if (i + 1) % 10 == 0 or i == total - 1:
                    percent_done = (i + 1) / total * 100
                    print(f"Progress: {i + 1}/{total} ({percent_done:.2f}%)")

        return coordinates, errors

    # Call the function and get coordinates
    coordinates, errors = fetch_coordinates_concurrently(gene_data)

    # Calculate success rate for all genes and "G protein" group names (including cached ones)
    total_genes = len(gene_data)
    success_genes = len(coordinates) + len(existing_coordinates)  # Count existing as success
    
    # Case-insensitive match for "G protein" in the group name (includes newly fetched and cached genes)
    g_protein_genes = [gene for gene in gene_data if "g protein" in gene['group_name'].lower()]
    
    # Include cached successes for "G protein" genes by checking the group name
    g_protein_success = len([gene for gene in coordinates if "g protein" in coordinates[gene]['group_name'].lower()]) + \
                        len([gene for gene, group_name in existing_coordinates.items() if "g protein" in group_name.lower()])

    
    g_protein_success_rate = (g_protein_success / len(g_protein_genes)) * 100 if g_protein_genes else 0
    
    # Print the final success rates
    print(f"Successfully retrieved data for {success_genes} genes ({(success_genes / total_genes) * 100:.2f}%)")
    print(f"Successfully retrieved data for {g_protein_success} G protein genes ({g_protein_success_rate:.2f}%)")
    print(f"Errors encountered for {len(errors)} genes.")

    # Append new coordinates to the existing file without deleting it
    with open(output_file, 'a') as f:
        writer = csv.writer(f, delimiter='\t')
        if os.stat(output_file).st_size == 0:  # If the file is empty, write header
            writer.writerow(['Ensembl gene ID', 'Chromosome', 'Start', 'End', 'Group name'])
        for ensembl_id, coord in coordinates.items():
            if ensembl_id not in existing_coordinates:
                writer.writerow([ensembl_id, coord['chromosome'], coord['start'], coord['end'], coord['group_name']])

    print(f"New coordinates appended to {output_file}")

    # Optional: Save errors to a separate file
    error_file = "gene_errors.log"
    with open(error_file, 'w') as f:
        for ensembl_id, error in errors.items():
            f.write(f"{ensembl_id}: {error}\n")

    print(f"Error log saved to {error_file}")

    # Function to verify all genes are in the output file
    def verify_genes():
        with open(output_file, 'r') as f:
            existing_gene_ids = {row['Ensembl gene ID'] for row in csv.DictReader(f, delimiter='\t')}
        
        missing_genes = [gene for gene in gene_data if gene['ensembl_gene_id'] not in existing_gene_ids]

        if missing_genes:
            print(f"Missing genes: {len(missing_genes)}")
            for gene in missing_genes:
                print(f"Missing gene: {gene['ensembl_gene_id']}, Group name: {gene['group_name']}")
        else:
            print("All genes are present in the output file.")

    # Verify genes
    verify_genes()

if __name__ == "__main__":
    fetch_gene_coordinates()
