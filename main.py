import time
import argparse
from fetch_gene_coordinates import fetch_gene_coordinates
from check_overlap_with_superdups import check_overlap_with_superdups
from check_overlap_with_sd import check_overlap_with_sd

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run the overlap check")
    parser.add_argument("--segdup_file", type=str, help="Path to the segdup file")
    args = parser.parse_args()

    # Run fetch_gene_coordinates
    print("Running fetch_gene_coordinates()...")
    fetch_gene_coordinates()

    time.sleep(3)

    # Run the function again, in case anything was missed
    print("Running fetch_gene_coordinates() again...")
    fetch_gene_coordinates()

    time.sleep(1)

    # Decide which function to run based on the input flag
    if args.segdup_file:
        print(f"Running check_overlap_with_sd() with {args.segdup_file}...")
        check_overlap_with_sd(args.segdup_file)
    else:
        print("Running check_overlap_with_superdups()...")
        check_overlap_with_superdups()

    print("Execution complete.")

if __name__ == "__main__":
    main()
