def score_motifs(occurance_w, length_w, motif_counts):
    # Sorting the motifs by count (descending order)
    # Sort motifs by score: (Occurrences * Length of motif)
    sorted_motifs = sorted(motif_counts.items(), key=lambda x: occurance_w*x[1] * length_w*len(x[0]), reverse=True)
    print(len(sorted_motifs))

    # Print the sorted motifs with the custom score
    for motif, count in sorted_motifs:
        print(f"Motif: {motif[0]}, Count: {count}, Length: {len(motif)}, Score: {count * len(motif)}")
