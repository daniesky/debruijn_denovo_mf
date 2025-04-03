import argparse
from motif_discovery import find_motifs
from result_analysis import score_motifs
from fasta_parser import FastaParser

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
    parser = FastaParser(args.fasta_file, args.limit)
    sequences = parser.sequences

    candidate_motifs = find_motifs(sequences, allow_gaps = args.gaps, k= args.k, threshold=0.3, overlap_factor=0.3, apply_hamming_distance=True, limit=10)
    score_motifs(candidate_motifs, sequences, k=args.k, entropy_weight=1, occurance_weight=0, length_weight=0, save_logos=3, print_motifs=5)

if __name__ == "__main__":
    main()