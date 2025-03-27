import networkx as nx
from matplotlib import pyplot as plt
from itertools import combinations
class DeBruijnGraph:
    def __init__(self, sequences, k, allow_gaps=False, kmer_mismatch_length=None):
        """
        Initializes the DeBruijnGraph object with sequences, k-mer length, and a flag for allowing gaps.
        
        Parameters:
        - sequences: A list of DNA sequences (strings).
        - k: The length of the k-mers for constructing the De Bruijn graph.
        - allow_gaps: A boolean flag indicating whether gaps should be allowed.
        - kmer_mismatch_length: The number of mismatches allowed in the gapped k-mers (only used if allow_gaps=True).
        """
        self.sequences = sequences
        self.k = k
        self.allow_gaps = allow_gaps
        self.kmer_mismatch_length = kmer_mismatch_length if kmer_mismatch_length is not None else k // 2
        self.graph = nx.DiGraph()

        self.build_graph()

    def build_graph(self):
        """
        Builds the De Bruijn graph (regular or gapped) based on the provided DNA sequences.
        Adjusts edge weights based on the occurrence of k-mers and applies a gap penalty if needed.
        
        Returns:
        - graph: The NetworkX DiGraph object representing the De Bruijn graph.
        """
        print(self.allow_gaps)
        print(self.kmer_mismatch_length)
        for seq_index, dna_sequence in enumerate(self.sequences):
            for i in range(len(dna_sequence) - self.k + 1):
                k_mer_start = dna_sequence[i:i + self.k - 1]
                k_mer_end = dna_sequence[i + 1:i + self.k]
                if not self.kmer_filter(dna_sequence[i:i + self.k]):
                    continue

                if self.allow_gaps:
                    # Use masked k-mers with allowed mismatches if gaps are allowed
                    kmer = dna_sequence[i:i + self.k]
                    masked_kmers = self.get_masked_kmers(kmer, self.kmer_mismatch_length)
                    for masked_kmer in masked_kmers:
                        k_mer_start = masked_kmer[:self.k - 1]
                        k_mer_end = masked_kmer[1:]
                        self._add_or_increment_edge(k_mer_start, k_mer_end)
                else:
                    # Regular De Bruijn graph construction without gaps
                    self._add_or_increment_edge(k_mer_start, k_mer_end, seq_index, i)

        self._remove_isolated_nodes_and_edges()
        self._apply_gap_penalty()

        return self.graph

    def _add_or_increment_edge(self, k_mer_start, k_mer_end, seq_index, i):
        """
        Helper function to add an edge to the graph or increment its weight.
        
        Parameters:
        - k_mer_start: The start of the k-mer.
        - k_mer_end: The end of the k-mer.
        """
        if self.graph.has_edge(k_mer_start, k_mer_end):
            self.graph[k_mer_start][k_mer_end]['weight'] += 1
            self.graph[k_mer_start][k_mer_end]['occurances'].append((seq_index, i))
        else:
            self.graph.add_edge(k_mer_start, k_mer_end, weight=1)
            self.graph[k_mer_start][k_mer_end]['occurances'] = [(seq_index, i)]


    def kmer_filter(self, kmer):
        """
        Filters k-mers based on a given condition (you can modify this method for custom filtering).
        
        Parameters:
        - kmer: The k-mer string to be filtered.
        
        Returns:
        - bool: Whether the k-mer passes the filter condition.
        """
        # Example filter: Only allow k-mers that consist of A, T, G, C (this could be customized)
        return all(base in 'ATGC' for base in kmer)

    def get_masked_kmers(self, kmer, mismatch_length):
        """
        Generates masked k-mers with up to `mismatch_length` mismatches in the original k-mer.
        
        Parameters:
        - kmer: The k-mer string to be masked.
        - mismatch_length: The number of allowed mismatches in the masked k-mers.
        
        Returns:
        - list: A list of masked k-mers.
        """
        masked_kmers = set()
        kmer_length = len(kmer)

        # Generate all possible bitmasks with at most kmer_mask_length wildcards
        for num_wildcards in range(0, mismatch_length + 1):
            for wildcard_positions in combinations(range(kmer_length), num_wildcards):
                masked_kmer = list(kmer)
                for pos in wildcard_positions:
                    masked_kmer[pos] = '*'
                masked_kmers.add(''.join(masked_kmer))

        return masked_kmers

    def _remove_isolated_nodes_and_edges(self):
        """
        Removes all nodes with degree 0 and edges with weight 0 from the graph.
        """
        nodes_to_remove = [node for node in self.graph.nodes() if self.graph.degree(node) == 0]
        self.graph.remove_nodes_from(nodes_to_remove)
        edges_to_remove = [(u, v) for u, v in self.graph.edges() if self.graph[u][v]['weight'] == 0]
        self.graph.remove_edges_from(edges_to_remove)

    def _apply_gap_penalty(self):
        """
        Applies a gap penalty to all edges with wildcards in them by reducing their weight by half.
        """
        if self.allow_gaps:
            for u, v, data in self.graph.edges(data=True):
                if '*' in u or '*' in v:
                    data['weight'] *= 0.5


    def draw_graph(self):
        """
        Draws the De Bruijn graph with edge weights using matplotlib.
        
        Parameters:
        - graph: A NetworkX DiGraph object representing the De Bruijn graph.
        """
        pos = nx.circular_layout(self.graph)
        plt.figure(figsize=(15, 15))
        
        # Draw the nodes
        nx.draw_networkx_nodes(self.graph, pos, node_size=700)
        
        # Draw the edges with labels indicating edge weights
        edge_labels = nx.get_edge_attributes(self.graph, 'weight')
        nx.draw_networkx_edges(self.graph, pos, edgelist=self.graph.edges(), arrowstyle='->', arrowsize=20)
        nx.draw_networkx_labels(self.graph, pos, font_size=10, font_family="sans-serif")
        
        # Add edge labels (weights)
        nx.draw_networkx_edge_labels(self.graph, pos, edge_labels=edge_labels)
        
        plt.title("De Bruijn Graph with Edge Weights")
        plt.show()
