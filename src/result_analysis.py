import pandas as pd
import logomaker
def motif_logo(motif, chip_seq_data, occurances):
    """
    The first node in the traversal path will be the defining points for the sequences 
    which contain the motifs.
    This means that the occurance list for the starting node, will define these 
    sequences, regardless if they contain the entire motif. 
    """
    dataframe = pd.DataFrame(0, index=range(len(motif)), columns=['A', 'C', 'G', 'T'])
    for seq_index, pos in occurances:
        sample = chip_seq_data[seq_index]
        for i in range(len(motif)):
            if len(motif) + pos <= len(sample):
                if sample[pos + i] in ('A', 'C', 'G', 'T'):
                    dataframe.at[i, sample[pos + i]] = dataframe.at[i, sample[pos + i]] + 1

    
    # Normalize the counts to get frequencies
    dataframe = dataframe.div(dataframe.sum(axis=1), axis=0)
    dataframe = dataframe.astype('float64')

    # create and style logo
    logo = logomaker.Logo(df=dataframe,
                fade_below=0.5,
                shade_below=0.5,
                figsize=(10,3))

    # set axes labels
    logo.ax.set_xlabel('Position',fontsize=14)
    logo.ax.set_ylabel("Frequency", labelpad=-1,fontsize=14)
    fig = logo.ax.get_figure()
    fig.savefig(f"{motif}.png", dpi=300, bbox_inches='tight')
        