from collections import Counter
import re
from de_bruijn_graph import DeBruijnGraph
from fasta_parser import FastaParser


def find_motifs(file, allow_gaps, k, max_read, apply_hamming_distance = False, kmer_mismatch_length = None, threshold = 0.5):
    parser = FastaParser(file, max_read)
    sequences = parser.sequences
    print(len(sequences))
    # Construct the De Bruijn graph
    graph = DeBruijnGraph(sequences, k=k, allow_gaps=allow_gaps, kmer_mismatch_length=kmer_mismatch_length).graph
    if apply_hamming_distance:
        apply_hamming_reward_after_creation(graph)
    # Discover motifs in the graph
    reconstructed_seq = motif_discovery(graph, threshold)
    for i, (seq, weight) in enumerate(reconstructed_seq):
        print(f"{i+1}. Sequence: {seq}, Total Weight: {weight/len(seq)}")


    # Count the occurrences of the motifs in the sequences
    motif_counts = count_occurances(sequences, reconstructed_seq)
    # Sorting the motifs by count (descending order)
    # Sort motifs by score: (Occurrences * Length of motif)
    sorted_motifs = sorted(motif_counts.items(), key=lambda x: x[1] * len(x[0]), reverse=True)

    # Print the sorted motifs with the custom score
    for motif, count in sorted_motifs:
        print(f"Motif: {motif}, Count: {count}, Length: {len(motif)}, Score: {count * len(motif)}")

    return motif_counts

def motif_discovery(graph, threshold):

    visited_edges = set()
    reconstructed_seq = []
    max_weight = max([d['weight'] for u, v, d in graph.edges(data=True)])
    threshold = max_weight * threshold

    while True:
        # Find all unvisited edges
        unvisited_edges = [(u, v, d) for u, v, d in graph.edges(data=True) if (u, v) not in visited_edges]
        
        if not unvisited_edges:  # Stop if all edges have been visited
            break

        # Get the max-weight edge from unvisited edges
        max_edge = max(unvisited_edges, key=lambda x: x[2]['weight'])
        if max_edge[2]['weight'] < threshold:
            break
        root_node = max_edge[0]
        seq = root_node  # Start the sequence from the source of the max edge
        accumulated_weight = 0
        
        node = root_node
        while True:
            # Get the highest-weight outgoing edge
            neighbors = [(neighbor, graph[node][neighbor]['weight']) for neighbor in graph.successors(node) if (node, neighbor) not in visited_edges]
            
            if not neighbors:
                break  # Stop when no unvisited edges remain
            
            # Select the highest-weight edge
            next_node, edge_weight = max(neighbors, key=lambda x: x[1])
            visited_edges.add((node, next_node))  # Mark edge as visited
            if edge_weight < threshold:
                break
            accumulated_weight += edge_weight
            seq += next_node[-1]  # Add last character of the next node
            node = next_node  # Move to the next node

        # Store the sequence and its weight
        reconstructed_seq.append((seq, accumulated_weight))

    # Sort sequences based on accumulated weight per nucleotide in descending order
    reconstructed_seq.sort(key=lambda x: x[1] / len(x[0]), reverse=True)

    return reconstructed_seq



def count_occurances(sequences, motifs):
    motif_counts = Counter()
    # Create regex patterns by replacing '*' with '.' (though there may be no '*' in this case)
    motif_patterns = {motif: re.compile(motif.replace('*', '.')) for motif, _ in motifs}
    
    for seq in sequences:
        for motif, pattern in motif_patterns.items():
            if pattern.search(seq):  # Use regex search to match the motif
                motif_counts[motif] += 1
                
    return motif_counts

def hamming_distance(kmer1, kmer2):
    """Calculate the Hamming distance between two k-mers."""
    return sum(c1 != c2 for c1, c2 in zip(kmer1, kmer2))

def reward_based_on_hamming(kmer1, kmer2, max_distance=3, reward_factor=2.0):
    """
    Reward the edge weight based on the Hamming distance between two k-mers.
    - max_distance: The maximum Hamming distance for a reward to be applied.
    - reward_factor: The factor by which to increase the edge weight.
    """
    hamming_dist = hamming_distance(kmer1, kmer2)
    
    if hamming_dist == 0:
        return reward_factor  # Identical k-mers, maximum reward
    elif hamming_dist == 1:
        return reward_factor * 1.5  # Small reward for one-character difference
    elif hamming_dist <= max_distance:
        # Gradually decrease the reward as the distance grows
        return reward_factor * (1.0 / (1 + hamming_dist))
    else:
        return 0  # No reward for larger Hamming distances
    
def apply_hamming_reward_after_creation(graph, max_distance=3, reward_factor=2.0):
    """
    Apply the Hamming distance-based reward to the edge weights of the graph.
    This function should be called after the graph has been created.
    
    - max_distance: The maximum Hamming distance for a reward to be applied.
    - reward_factor: The factor by which to increase the edge weight.
    """
    for u, v, _ in graph.edges(data=True):
        k_mer_start = u  # start k-mer of the edge
        k_mer_end = v    # end k-mer of the edge
        
        # Apply Hamming distance-based reward to the edge weight
        additional_weight = reward_based_on_hamming(k_mer_start, k_mer_end, max_distance, reward_factor)
        
        # Update the edge weight (add the reward to the existing weight)
        graph[u][v]['weight'] += additional_weight
