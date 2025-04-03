import logomaker
from collections import defaultdict
import numpy as np
import pandas as pd
def shannon_entropy(pfm):
    """
    Calculate the Shannon entropy of a position frequency matrix (PFM).
    :param pfm: Position frequency matrix (DataFrame)
    :return: Shannon entropy value
    """
    # Normalize the PFM to get probabilities
    pfm_normalized = pfm.div(pfm.sum(axis=1), axis=0)
    
    # Calculate Shannon entropy
    entropy = (pfm_normalized * pfm_normalized.apply(lambda x: x[x > 0].apply(lambda y: y * np.log2(y)))).sum(axis=1).sum()    
    return entropy


def score_motif(pfm, occurance_factor, k, motif_len, entropy_weight, occurance_weight, length_weight):
    return entropy_weight * (1-shannon_entropy(pfm)) + occurance_weight * occurance_factor + length_weight * (1 - (motif_len / k))
    
def score_motifs(motifs, sequences, k, entropy_weight, occurance_weight, length_weight, save_logos=3, print_motifs=5):
    scored_motifs = []
    for motif, _, _ in motifs:
        pfm, occurrences = compute_occurance_pfm(motif, sequences)
        score = score_motif(pfm, occurrences / len(sequences), k, len(motif), entropy_weight, occurance_weight, length_weight)
        scored_motifs.append((motif, score, pfm))
    
    # Sort motifs by score
    scored_motifs.sort(key=lambda x: x[1], reverse=True)
    # Print top motifs
    print("Top motifs:")
    for i, (motif, score, pfm) in enumerate(scored_motifs[:print_motifs]):
        print(f"Motif {i+1}: {motif} (Score: {score})")
    
    if save_logos > 0:
        print("Saving logos for top motifs...")
        for i, (motif, score, pfm) in enumerate(scored_motifs[:save_logos]):
            motif_logo_alignment_version(motif, sequences, pfm, occurrences)   
    return scored_motifs


def align_best_fit(sequence, pattern):
    if len(pattern) > len(sequence):
        raise ValueError("Pattern length cannot be greater than sequence length.")
    
    best_score = -1
    best_position = 0
    
    seq_len = len(sequence)
    pat_len = len(pattern)

    for i in range(seq_len - pat_len + 1):
        segment = sequence[i:i + pat_len]
        score = sum(1 for a, b in zip(segment, pattern) if a == b)

        if score > best_score:
            best_score = score
            best_position = i
    
    # Create alignment visualization
    alignment = ["-"] * len(sequence)
    alignment[best_position:best_position + pat_len] = pattern
    return best_position, best_score


def compute_occurance_pfm(motif, chip_seq_data):
    dataframe = pd.DataFrame(0, index=range(len(motif)), columns=['A', 'C', 'G', 'T'])
    occurances  = 0
    for seq in chip_seq_data:
        if not seq:
            continue
        # Align the motif with the sequence
        alignment, _ = align_best_fit(seq, motif)
        occurances += 1
        for i in range(len(motif)):
            if seq[alignment + i] in ('A', 'C', 'G', 'T'):
                dataframe.at[i, seq[alignment + i]] = dataframe.at[i, seq[alignment + i]] + 1

    # Normalize the counts by the sum of each row (not max)
    dataframe = dataframe.div(dataframe.sum(axis=1), axis=0)

    # Check if any elements are non finite
    if not np.isfinite(dataframe.values).all():
        print("Non-finite values found in the DataFrame.")
        # Handle non-finite values (e.g., replace with 0 or NaN)
    return dataframe, occurances


def motif_logo_alignment_version(motif, chip_seq_data, dataframe, occurances):
    # create and style logo
    logo = logomaker.Logo(df=dataframe,
                fade_below=0.5,
                shade_below=0.5,
                figsize=(10,3))

    # set axes labels
    logo.ax.set_xlabel('Position',fontsize=14)
    logo.ax.set_ylabel("Frequency", labelpad=-1,fontsize=14)

    fig = logo.ax.get_figure()
    fig.savefig(f"logos/{motif}.png", dpi=300, bbox_inches='tight')