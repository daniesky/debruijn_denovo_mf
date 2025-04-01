import argparse
from motif_discovery import find_motifs

def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Run motif discovery on a given FASTA file.")
    parser.add_argument('fasta_file', type=str, default="data/output.fa", help="Path to the input FASTA file.")
    parser.add_argument('--limit', type=int, default=10_000, help="Limit on the number of sequences to read (default: 300,000).")
    parser.add_argument('--k', type=int, default=5, help="k-mer length (default: 5).")
    parser.add_argument('--gaps', type=bool, default=False, help=" Allow gaps in the k-mers (default: False).")
    return parser.parse_args()

def main():
    # Parse command-line arguments
    args = parse_arguments()
    # Set k-mer length and FASTA file
    motif_counts = find_motifs(file = args.fasta_file, allow_gaps = args.gaps, k= args.k, max_read = args.limit, threshold=0.3, overlap_factor=0.4, occurance_threshold = 0.0, apply_hamming_distance=True)
    for key, value, _ in motif_counts:
        print(key, value, value/len(key))


if __name__ == "__main__":
    main()