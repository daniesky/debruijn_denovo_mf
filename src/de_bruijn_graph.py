import networkx as nx
from matplotlib import pyplot as plt
class DeBruijnGraph:
    def __init__(self, sequences, k):
        """
        Initializes the DeBruijnGraph object with sequences, k-mer length
        
        Parameters:
        - sequences: A list of DNA sequences (strings).
        - k: The length of the k-mers for constructing the De Bruijn graph.
        """
        self.sequences = sequences
        self.k = k
        self.graph = nx.DiGraph()

        self.build_graph()

    def build_graph(self):
        """
        Builds the De Bruijn graph based on the provided DNA sequences.
        
        Returns:
        - graph: The NetworkX DiGraph object representing the De Bruijn graph.
        """
        for seq_index, dna_sequence in enumerate(self.sequences):
            for i in range(len(dna_sequence) - self.k + 1):
                k_mer_start = dna_sequence[i:i + self.k - 1]
                k_mer_end = dna_sequence[i + 1:i + self.k]

                self._add_or_increment_edge(k_mer_start, k_mer_end, seq_index, i)

        self._remove_isolated_nodes_and_edges()

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

    def _remove_isolated_nodes_and_edges(self):
        """
        Removes all nodes with degree 0 and edges with weight 0 from the graph.
        """
        nodes_to_remove = [node for node in self.graph.nodes() if self.graph.degree(node) == 0]
        self.graph.remove_nodes_from(nodes_to_remove)
        edges_to_remove = [(u, v) for u, v in self.graph.edges() if self.graph[u][v]['weight'] == 0]
        self.graph.remove_edges_from(edges_to_remove)


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
