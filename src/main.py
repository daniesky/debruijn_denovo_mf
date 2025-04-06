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
    parser.add_argument('--inst_limit', type=float, default=0.0, help="Instance limit for motif scoring (default: 0.0).")
    parser.add_argument('--save_logos', type=int, default=10, help="Number of logos to save (default: 10).")
    parser.add_argument('--print_motifs', type=int, default=10, help="Number of motifs to print (default: 10).")
    parser.add_argument('--overlap_factor', type=float, default=0.25, help="Overlap factor for motif discovery (default: 0.25).")
    parser.add_argument('--threshold', type=float, default=0.25, help="Threshold for motif discovery (default: 0.25).")
    parser.add_argument('--scoring_limit', type=int, default=10, help="Limit for number of motifs to score (default: 10). High impact on runtime.")
    return parser.parse_args()

def main():
    # Parse command-line arguments
    args = parse_arguments()
    # Set k-mer length and FASTA file
    parser = FastaParser(args.fasta_file, args.limit)
    sequences = parser.sequences

    candidate_motifs = find_motifs(sequences, allow_gaps = args.gaps, k= args.k, threshold=args.threshold, overlap_factor=args.overlap_factor, limit = args.scoring_limit)
    score_motifs(candidate_motifs, sequences, save_logos=args.save_logos, print_motifs=args.print_motifs, inst_limit = args.inst_limit)

if __name__ == "__main__":
    main()