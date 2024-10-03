import time
from fetch_gene_coordinates import fetch_gene_coordinates
from check_overlap_with_superdups import check_overlap_with_superdups

if __name__ == "__main__":
    # Run fetch_gene_coordinates
    print("Running fetch_gene_coordinates()...")
    fetch_gene_coordinates()

    time.sleep(3)

    # Run the function again, in case anything was missed
    print("Running fetch_gene_coordinates() again...")
    fetch_gene_coordinates()

    time.sleep(1)

    # Run tcheck_overlap_with_superdups
    print("Running check_overlap_with_superdups()...")
    check_overlap_with_superdups()

    print("Execution complete.")
