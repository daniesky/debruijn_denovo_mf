from collections import Counter, defaultdict
from de_bruijn_graph import DeBruijnGraph
from fasta_parser import FastaParser
from itertools import product

def find_motifs(sequences, allow_gaps, k, threshold, overlap_factor, limit):

    # Construct the De Bruijn graph
    graph_obj = DeBruijnGraph(sequences, k=k)
    graph = graph_obj.graph
    if not allow_gaps:
        reconstructed_seq = strict_motif_discovery(graph, threshold, overlap_factor=overlap_factor)
    else:
        reconstructed_seq = ambig_motif_discovery(graph, threshold, overlap_factor=overlap_factor)

    return reconstructed_seq[:limit]

def strict_motif_discovery(graph, threshold, weight_reduction_factor=0.5, overlap_factor=0.3):
    reconstructed_seq = []
    max_weight = max([d['weight'] for u, v, d in graph.edges(data=True)])
    threshold = max_weight * threshold

    while True:
        # Get the max-weight edge from unvisited edges
        edge = max(graph.edges(data=True), key=lambda x: x[2]['weight'])
        node = backtrack_path(graph, edge[0], threshold)
        edge = max(
                ((neighbor, graph[node][neighbor]['weight']) for neighbor in graph.successors(node)),
                key=lambda x: x[1]
            )
        edge_weight = edge[1]
        if edge_weight < threshold:
            break
        seq = node
        accumulated_weight = 0
        traversed_nodes = [node]

        previous_positions = None    
        while True:
            # Select the highest-weight edge
            next_node, edge_weight = max(
                ((neighbor, graph[node][neighbor]['weight']) for neighbor in graph.successors(node)),
                key=lambda x: x[1]
            )

            # Compute the overlap between occurances in the previous and current edge. There needs to be a given overlap for traversal to continue
            overlap_count = -1
            if previous_positions:
                overlap_count = position_set_overlap(previous_positions, graph[node][next_node]['occurances'])

            # Reduce the weight of the traversed edge after visiting it
            graph[node][next_node]['weight'] *= weight_reduction_factor

            # Check if the edge weight is below the threshold or the overlap is too low to continue
            if edge_weight < threshold or (overlap_count != -1 and overlap_count < edge_weight*overlap_factor):
                break

            accumulated_weight += edge_weight
            seq += next_node[-1]  # Add last character of the next node
            previous_positions = graph[node][next_node]['occurances']
            traversed_nodes.append(next_node)  # Keep track of traversed nodes
            node = next_node  # Move to the next node

        # Store the sequence and its weight
        reconstructed_seq.append((seq, accumulated_weight, traversed_nodes))

    # Sort sequences based on accumulated weight per nucleotide in descending order keeping only the top 20
    reconstructed_seq.sort(key=lambda x: x[1] / len(x[0]), reverse=True)
    reconstructed_seq = reconstructed_seq[:20]

    return reconstructed_seq

def backtrack_path(graph, node, threshold, overlap_factor=0.3):
    """
    Backtrack the path from the given node in the graph until the edge weight falls below the threshold.
    """
    previous_positions = None
    visited_nodes = set()
    visited_nodes.add(node)

    while True:
        if not list(graph.predecessors(node)):
            return node
        # Select the highest-weight edge
        next_node, edge_weight = max(
            ((neighbor, graph[neighbor][node]['weight']) for neighbor in graph.predecessors(node) if neighbor not in visited_nodes),
            key=lambda x: x[1]
        )
        if not next_node:
            return node
        
        # Compute the overlap between occurances in the previous and current edge. There needs to be a given overlap for traversal to continue
        overlap_count = -1
        if previous_positions:
            overlap_count = position_set_overlap(previous_positions, graph[next_node][node]['occurances'] )

        # Check if the edge weight is below the threshold or the overlap is too low to continue
        if edge_weight < threshold or (overlap_count != -1 and overlap_count < edge_weight*overlap_factor):
            return node
        visited_nodes.add(next_node)
        previous_positions = graph[next_node][node]['occurances']
        node = next_node  # Move to the next node
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

def ambig_motif_discovery(graph, threshold, similarity_threshold=0.75, weight_reduction_factor=0.5, overlap_factor=0.3):

    reconstructed_seq = []
    max_weight = max([d['weight'] for u, v, d in graph.edges(data=True)])
    threshold = max_weight * threshold

    while True:
        edge = max(graph.edges(data=True), key=lambda x: x[2]['weight'])
        node = backtrack_path(graph, edge[0], threshold)
        edge = max(
                ((neighbor, graph[node][neighbor]['weight']) for neighbor in graph.successors(node)),
                key=lambda x: x[1]
            )
        edge_weight = edge[1]
        if edge_weight < threshold:
            break
        seq = node
        accumulated_weight = 0
        edge_origin = None
        traversed_nodes = [node]
        while True:
            # Get outgoing edges
            neighbors = [(neighbor, graph[node][neighbor]['weight']) 
                         for neighbor in graph.successors(node) if (node, neighbor)]
            
            if not neighbors:
                break  # Stop when no unvisited edges remain

            base_weight_map = defaultdict(float)
            for neighbor, weight in neighbors:
                base_weight_map[neighbor] = weight  # Sum weights for each nucleotide
            if edge_origin:
                valid_neighbors = {base for base, weight in base_weight_map.items() if weight >= threshold and position_set_overlap(edge_origin, graph[node][base]['occurances']) > weight*overlap_factor}
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
                    if total_weight >= threshold and position_set_overlap(edge_origin, positions) > total_weight*overlap_factor:
                        # Check similarity condition
                        min_weight = min(base_weight_map[b] for b in candidate_bases)
                        max_weight = max(base_weight_map[b] for b in candidate_bases)
                        similarity_ratio = min_weight / max_weight if max_weight > 0 else -1

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
            edge_origin = graph[node][next_node]['occurances']
            valid_neighbor_bases = {seq[-1] for seq in valid_neighbors}
            iupac_code = IUPAC_CODES[frozenset(valid_neighbor_bases)] if len(valid_neighbor_bases) > 1 else next(iter(valid_neighbor_bases))
            seq += iupac_code  # Add to sequence

            accumulated_weight += graph[node][next_node]['weight']
            traversed_nodes.append(next_node)  # Keep track of traversed nodes
            node = next_node

        # Store the sequence and its weight
        reconstructed_seq.append((seq, accumulated_weight, traversed_nodes))

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