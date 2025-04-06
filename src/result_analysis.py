import logomaker
import numpy as np
import pandas as pd
def shannon_entropy(pfm):
    """
    Calculate the Shannon entropy of a position frequency matrix (PFM).
    :param pfm: Position frequency matrix (DataFrame)
    :return: Shannon entropy value
    """
    # Add a small epsilon to avoid log2(0)
    eps = 1e-10
    entropy_per_position = -np.sum(pfm * np.log2(pfm + eps), axis=1)
    max_entropy = np.log2(pfm.shape[1])  # log2(4) = 2 for DNA
    normalized_entropy = (entropy_per_position / max_entropy).sum() / pfm.shape[0]
    return normalized_entropy

def score_motif(pfm, consensus_score):
    return (1-shannon_entropy(pfm))*(consensus_score)

def score_motifs(motifs, sequences, save_logos, print_motifs, inst_limit):
    scored_motifs = []
    for motif, _, _ in motifs:
        pfm, occurrences, consensus_score = compute_occurance_pfm(motif, sequences, inst_limit)
        score = score_motif(pfm, consensus_score)
        scored_motifs.append((motif, score, pfm, consensus_score))
    
    # Sort motifs by score
    scored_motifs.sort(key=lambda x: x[1], reverse=True)
    # Print top motifs
    print("Top motifs:")
    for i, (motif, score, pfm, consensus_score) in enumerate(scored_motifs[:print_motifs]):
        print(f"Motif {i+1}: {motif} (Score: {score}), Consensus Score: {consensus_score}")
    
    if save_logos > 0:
        print("Saving logos for top motifs...")
        for i, (motif, score, pfm, _) in enumerate(scored_motifs[:save_logos]):
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


def compute_occurance_pfm(motif, chip_seq_data, inst_limit):
    dataframe = pd.DataFrame(0, index=range(len(motif)), columns=['A', 'C', 'G', 'T'])
    threshold = int(len(motif)*inst_limit)
    occurances  = 0
    consensus_score = 0
    for seq in chip_seq_data:
        if not seq:
            continue
        # Align the motif with the sequence
        alignment, score = align_best_fit(seq, motif)
        if(score < threshold):
            continue
        consensus_score += score
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
    return dataframe, occurances, consensus_score


def motif_logo_alignment_version(motif, chip_seq_data, dataframe, occurances):
    # create and style logo
    logo = logomaker.Logo(df=dataframe,
                fade_below=0.5,
                shade_below=0.5,
                figsize=(10,3))

    # set axes labels
    logo.ax.set_xlabel('Position',fontsize=14)
    logo.ax.set_ylabel("Frequency", labelpad=-1,fontsize=14)
    logo.ax.set_title(f"Motif: {motif} occurrances/sequences: {occurances/len(chip_seq_data)}", fontsize=14)
    fig = logo.ax.get_figure()
    fig.savefig(f"logos/{motif}.png", dpi=300, bbox_inches='tight')