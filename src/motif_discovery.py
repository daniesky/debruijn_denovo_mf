from collections import Counter, defaultdict
import re
from de_bruijn_graph import DeBruijnGraph
from fasta_parser import FastaParser
from itertools import product


def find_motifs(file, allow_gaps, k, max_read, apply_hamming_distance = False, threshold = 0.5):
    parser = FastaParser(file, max_read)
    sequences = parser.sequences

    # Construct the De Bruijn graph
    graph_obj = DeBruijnGraph(sequences, k=k, allow_gaps=False)
    graph = graph_obj.graph
    if apply_hamming_distance:
        apply_hamming_reward_after_creation(graph)
    if not allow_gaps:
        reconstructed_seq = motif_discovery(graph, threshold)
    else:
        reconstructed_seq = iuapac_motif_discovery(graph, threshold)

    return reconstructed_seq

def motif_discovery(graph, threshold, weight_reduction_factor=0.5):
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

        previous_positions = None    
        while True:
            # Get the highest-weight outgoing edge
            neighbors = [(neighbor, graph[node][neighbor]['weight']) 
                         for neighbor in graph.successors(node) if (node, neighbor) not in visited_edges]
            
            if not neighbors:
                break  # Stop when no unvisited edges remain


            # Select the highest-weight edge
            next_node, edge_weight = max(neighbors, key=lambda x: x[1])

            #Now we want to check if this edge has an actual overlap with the previous edge. If it does not, we do not want to keep traversing.
            overlap_count = 0
            if previous_positions:
                overlap_count = position_set_overlap(previous_positions, graph[node][next_node]['occurances'])

            visited_edges.add((node, next_node))  # Mark edge as visited
            graph[node][next_node]['weight'] *= weight_reduction_factor
            if edge_weight < threshold or (overlap_count != 0 and overlap_count < edge_weight*0.3):
                break
            accumulated_weight += edge_weight
            seq += next_node[-1]  # Add last character of the next node
            previous_positions = graph[node][next_node]['occurances']
            node = next_node  # Move to the next node

        # Store the sequence and its weight
        reconstructed_seq.append((seq, accumulated_weight))

    # Sort sequences based on accumulated weight per nucleotide in descending order keeping only the top 20
    reconstructed_seq.sort(key=lambda x: x[1] / len(x[0]), reverse=True)
    reconstructed_seq = reconstructed_seq[:20]

    return reconstructed_seq


IUPAC_CODES = {
    frozenset(["A"]): "A",
    frozenset(["C"]): "C",
    frozenset(["G"]): "G",
    frozenset(["T"]): "T",
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W",
    frozenset(["G", "T"]): "K",
    frozenset(["A", "C"]): "M",
    frozenset(["C", "G", "T"]): "B",
    frozenset(["A", "G", "T"]): "D",
    frozenset(["A", "C", "T"]): "H",
    frozenset(["A", "C", "G"]): "V",
    frozenset(["A", "C", "G", "T"]): "N",
}

def iuapac_motif_discovery(graph, threshold, similarity_threshold=0.95, weight_reduction_factor=0.5):

    reconstructed_seq = []
    max_weight = max([d['weight'] for u, v, d in graph.edges(data=True)])
    threshold = max_weight * threshold

    while True:
        # Get the max-weight edge from unvisited edges
        max_edge = max(graph.edges(data=True), key=lambda x: x[2]['weight'])
        if max_edge[2]['weight'] < threshold:
            break
        node = max_edge[0]
        seq = node  # Start the sequence from the source of the max edge
        accumulated_weight = 0
        previous_positions = None
        
        while True:
            # Get unvisited outgoing edges
            neighbors = [(neighbor, graph[node][neighbor]['weight']) 
                         for neighbor in graph.successors(node) if (node, neighbor)]
            
            if not neighbors:
                break  # Stop when no unvisited edges remain


            base_weight_map = defaultdict(float)
            for neighbor, weight in neighbors:
                base_weight_map[neighbor] = weight  # Sum weights for each nucleotide
            if previous_positions:
                valid_neighbors = {base for base, weight in base_weight_map.items() if weight >= threshold and position_set_overlap(previous_positions, graph[node][base]['occurances']) > weight*0.3}
            else:
                valid_neighbors = {base for base, weight in base_weight_map.items() if weight >= threshold}
            # Check if any IUPAC code passes the threshold
            if not valid_neighbors:
                sorted_bases = sorted(base_weight_map.items(), key=lambda x: x[1], reverse=True)
                total_weight = 0
                candidate_bases = set()
                positions = list()
                for base, weight in sorted_bases:
                    total_weight += weight
                    positions.extend(graph[node][base]['occurances'])
                    candidate_bases.add(base)
                    if total_weight >= threshold and position_set_overlap(previous_positions, positions) > total_weight*0.3:
                        # Check similarity condition
                        min_weight = min(base_weight_map[b] for b in candidate_bases)
                        max_weight = max(base_weight_map[b] for b in candidate_bases)
                        similarity_ratio = min_weight / max_weight if max_weight > 0 else 1

                        if similarity_ratio >= similarity_threshold:
                            valid_neighbors = candidate_bases

                        break
            
            
            # Use IUPAC encoding if multiple bases are valid
            if not valid_neighbors:
                break

            # Reduce the weight of the traversed edge after visiting it
            for neighbor in valid_neighbors:
                    graph[node][neighbor]['weight'] *= weight_reduction_factor
            # Select the next node based on max edge weight among valid bases
            next_node = max(valid_neighbors,
                            key=lambda x: graph[node][x]['weight'])
            previous_positions = graph[node][next_node]['occurances']
            valid_neighbor_bases = {seq[-1] for seq in valid_neighbors}
            iupac_code = IUPAC_CODES[frozenset(valid_neighbor_bases)] if len(valid_neighbor_bases) > 1 else next(iter(valid_neighbor_bases))
            seq += iupac_code  # Add to sequence

            accumulated_weight += graph[node][next_node]['weight']

            node = next_node

        # Store the sequence and its weight
        reconstructed_seq.append((seq, accumulated_weight))

    # Sort sequences based on accumulated weight per nucleotide in descending order
    reconstructed_seq.sort(key=lambda x: x[1] / len(x[0]), reverse=True)

    return reconstructed_seq



# Reverse lookup from nucleotide set to IUPAC code
IUPAC_MAP = {value: key for key, value in IUPAC_CODES.items()}

def expand_iupac(motifs):
    """Expands a sequence with IUPAC codes into all possible combinations of ACGT."""
    
    # Go through each motif
    for motif, _ in motifs:
        # Create a list of nucleotide sets for each character in the motif
        nucleotide_sets = [IUPAC_MAP.get(base, [base]) for base in motif]
        
        # Generate and yield each expanded motif on the fly (avoiding memory overhead)
        for comb in product(*nucleotide_sets):
            yield ''.join(comb)

def count_occurances(sequences, motifs):
    motif_counts = Counter()  # To hold the counts of original motifs
    
    # Step 1: Generate the expanded motifs
    expanded_motifs = {}
    for motif,_ in motifs:
        expanded_motifs[motif] = list(expand_iupac([(motif, None)]))  # Create expanded motifs for each original motif
    
    # Step 2: Iterate through the sequences and check for matches
    for seq in sequences:
        for motif, expanded in expanded_motifs.items():
            # For each expanded motif, check if it matches the sequence
            for expanded_motif in expanded:
                if expanded_motif in seq:  # If the expanded motif is found in the sequence
                    motif_counts[motif] += 1  # Count it under the original motif
                    break  # Break after the first match to avoid double counting lines
    
    return motif_counts


def position_set_overlap(set1, set2):
    shifted_prev = {(seq, pos+1) for seq, pos in set1}
    return sum(1 for seq, pos in set2 if (seq, pos) in shifted_prev)

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
