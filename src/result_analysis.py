import pandas as pd
import logomaker
from collections import defaultdict


#def score_motifs(motifs, sequences, occurrence_threshold, save_logos=3):
#   for motif in motifs:
#        pfm, occurrences = compute_occurance_pfm(motif, sequences, occurrence_threshold)

        
   
def align_best_fit(sequence, pattern):
    if len(pattern) > len(sequence):
        print(pattern)
        print(sequence)
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


def compute_occurance_pfm(motif, chip_seq_data, occurance_threshold):
    threshold = len(motif) * occurance_threshold
    dataframe = pd.DataFrame(0, index=range(len(motif)), columns=['A', 'C', 'G', 'T'])
    occurances  = 0
    for seq in chip_seq_data:
        if not seq:
            continue
        # Align the motif with the sequence
        alignment, score = align_best_fit(seq, motif)
        if score > threshold:
            occurances += 1
            for i in range(len(motif)):
                if seq[alignment + i] in ('A', 'C', 'G', 'T'):
                    dataframe.at[i, seq[alignment + i]] = dataframe.at[i, seq[alignment + i]] + 1

    # Normalize the counts by the sum of each row (not max)
    dataframe = dataframe.div(dataframe.sum(axis=1), axis=0)

    return dataframe, occurances


def motif_logo_alignment_version(motif, chip_seq_data, occurance_threshold):
    # create and style logo
    dataframe, occurances = compute_occurance_pfm(motif, chip_seq_data, occurance_threshold)
    logo = logomaker.Logo(df=dataframe,
                fade_below=0.5,
                shade_below=0.5,
                figsize=(10,3))

    # set axes labels
    logo.ax.set_xlabel('Position',fontsize=14)
    logo.ax.set_ylabel("Frequency", labelpad=-1,fontsize=14)

    logo.ax.set_title(f"Motif Logo for {motif}; Occurance percentage: {occurances / (len(chip_seq_data)-1)}", fontsize=16, fontweight='bold')

    fig = logo.ax.get_figure()
    fig.savefig(f"logos/{motif}.png", dpi=300, bbox_inches='tight')