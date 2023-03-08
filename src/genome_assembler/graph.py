"""Main file of assembly module with de Bruijn Graph and Node classes."""
import sys
import copy
from typing import List, Dict,Set, Tuple, TypeVar, Type

import graphviz
import numpy as np

from genome_assembler import utils
from genome_assembler import error_correction
from genome_assembler import scs
from genome_assembler.analyze_training_data import parse_fasta

sys.setrecursionlimit(100000)

S = TypeVar("S", bound="DeBruijn")


class NotEulerianGraphError(AttributeError):
    """Raised when graph is expected to be eulerian but isn't"""


class Node:
    """Represents a node in de Gruijn graphself.

    Each node in de Bruijn graph represents a k-1 mer of a string
    """

    def __init__(self, k_minus_1_mer: str):
        self.k_minus_1_mer: str = k_minus_1_mer
        self._indegree: int = 0
        self._outdegree: int = 0
        self.predecessors: Set[Node] = set()

    def __hash__(self) -> int:
        return hash(self.k_minus_1_mer)

    def __str__(self) -> str:
        return self.k_minus_1_mer

    @property
    def indegree(self) -> int:
        """Number of input egdes."""
        return self._indegree

    @indegree.setter
    def indegree(self, value: int) -> None:
        self._indegree = value

    @property
    def outdegree(self) -> int:
        """Number of output edges"""
        return self._outdegree

    @outdegree.setter
    def outdegree(self, value: int) -> None:
        self._outdegree = value

    def is_balanced(self) -> bool:
        """Returs true if Node is is balances."""
        return self._indegree == self._outdegree

    def is_semi_balanced(self) -> bool:
        """Returns true if Node is semi-balaned."""
        return abs(self._indegree - self._outdegree) == 1

    def is_tip(self) -> bool:
        """Returns True if node is the end of a ,,tip'' on the graph.

        In other words returns True if node's outdegree or indegree is 0."""
        return self.indegree == 0 or self.outdegree == 0

    def copy(self):
        return copy.copy(self)


AdjacencyList = Dict[Node, Dict[Node, int]]


class DeBruijn:
    """Represent de Bruijn graph"""

    def _update_nodes(self, k_minus_1_mer: str) -> Node:
        if k_minus_1_mer not in self.nodes:
            self.nodes[k_minus_1_mer] = Node(k_minus_1_mer)
        return self.nodes[k_minus_1_mer]

    def _balance_check(self) -> None:
        self.balanced, self.semi_balanced, self.unbalanced = 0, 0, 0
        for node in self.nodes.values():
            if node.is_balanced():
                self.balanced += 1
            elif node.is_semi_balanced():
                self.semi_balanced += 1
            else:
                self.unbalanced += 1

    def _update_tips(self) -> None:
        self.tips = set()
        for node in self.nodes.values():
            if node.indegree == 0:
                self.heads.add(node)
            if node.outdegree == 0:
                self.tails.add(node)
            if node.is_tip():
                self.tips.add(node)

    def _to_dot(self, weights: bool = False) -> graphviz.Digraph:
        graph = graphviz.Digraph(comment="DeBruijn graph")
        for node in self.graph:
            graph.node(node.k_minus_1_mer, node.k_minus_1_mer)
        for source, destinations in self.graph.items():
            if weights:
                for destination, weight in destinations.items():
                    graph.edge(
                        source.k_minus_1_mer,
                        destination.k_minus_1_mer,
                        label=str(weight),
                    )
            else:
                for destination in destinations:
                    graph.edge(source.k_minus_1_mer, destination.k_minus_1_mer)
        return graph

    def __init__(self) -> None:
        self.nodes: Dict[str, Node] = {}
        self.graph: AdjacencyList = {}
        self.balanced: int = 0
        self.semi_balanced: int = 0
        self.unbalanced: int = 0
        self.heads: Set[Node] = set()
        self.tails: Set[Node] = set()
        self.tips: Set[Node] = set()
        self.kmer_hist = {}

    @classmethod
    def from_reads(cls: Type[S], reads: List[str], k: int) -> S:
        """Create a de Bruijn Graph from list of strings."""
        new_graph = cls()
        for read in reads:
            for k_mer, _ in utils.chop(read, k):
                new_graph.kmer_hist.setdefault(k_mer, 0)
                new_graph.kmer_hist[k_mer] += 1
                k_mer_L, k_mer_R = k_mer[:-1], k_mer[1:]
                nodeL = new_graph._update_nodes(k_mer_L)
                nodeR = new_graph._update_nodes(k_mer_R)
                new_graph.graph.setdefault(nodeL, {}).setdefault(nodeR, 0)
                new_graph.graph[nodeL][nodeR] += 1
                new_graph.graph.setdefault(nodeR, {})
                nodeR.predecessors.add(nodeL)
                nodeL.outdegree += 1
                nodeR.indegree += 1
        new_graph._balance_check()
        new_graph._update_tips()
        return new_graph

    @classmethod
    def from_graph(cls: Type[S], graph: AdjacencyList) -> S:
        """Creates a de Bruijn graph from weighted adjacency list representation."""
        new_graph = cls()
        new_graph.graph = graph.copy()
        for node in new_graph.graph:
            new_graph.nodes[node.k_minus_1_mer] = node
        new_graph._balance_check()
        new_graph._update_tips()
        return new_graph

    @classmethod
    def from_kmer_hist(cls: Type[S], kmer_hist: Dict[str, int]) -> S:
        """Creates a de Bruijn graph from k_mer spectrum hashmap."""
        new_graph = cls()
        new_graph.kmer_hist = kmer_hist
        for k_mer, count in kmer_hist.items():
            k_mer_L, k_mer_R = k_mer[:-1], k_mer[1:]
            nodeL = new_graph._update_nodes(k_mer_L)
            nodeR = new_graph._update_nodes(k_mer_R)
            new_graph.graph.setdefault(nodeL, {}).setdefault(nodeR, count)
            new_graph.graph.setdefault(nodeR, {})
            nodeR.predecessors.add(nodeL)
            nodeL.outdegree += count
            nodeR.indegree += count
        new_graph._balance_check()
        new_graph._update_tips()
        return new_graph

    @property
    def node_count(self) -> int:
        """Number of nodes in the graph."""
        return len(self.nodes)

    @property
    def edge_count(self) -> int:
        """Number of edges in the graph."""
        edges: int = 0
        for node_neighbours in self.graph.values():
            edges += len(node_neighbours)
        return edges

    def has_eulerian_walk(self) -> bool:
        """Returns True if the graph has eulerian walk."""
        return self.unbalanced == 0 and self.semi_balanced == 2

    def has_eulerian_cycle(self) -> bool:
        """Returns True if the graph has eulerian walk."""
        return self.unbalanced == 0 and self.semi_balanced == 0

    def is_eulerian(self) -> bool:
        """Returns true if the graph is eulerian."""
        return self.has_eulerian_walk() or self.has_eulerian_cycle()

    def show(self, filename: str, weights: bool = False) -> None:
        """Saves graph visualization to a file and opens it."""
        self._to_dot(weights=weights).view(filename)

    def get_sub_graph(
        self, supergraph: AdjacencyList, node: Node
    ) -> Tuple[AdjacencyList, List[Node]]:
        """Returns a connected components containing the node from a super-graphself.

        Uses recursive dfs to do so."""
        discovered: List[Node] = []

        def _travel(node: Node):
            discovered.append(node)
            for neighbour in supergraph[node]:
                if neighbour not in discovered:
                    _travel(neighbour)

        _travel(node)
        sub_graph: Dict[Node, Dict[Node, int]] = {}
        for node in discovered:
            sub_graph.setdefault(node, self.graph[node])

        return sub_graph.copy(), discovered

    def connected_components(self) -> List[Dict[Node, Dict[Node, int]]]:
        """Returns a list of strongly connected components and saves it as an instance variable."""
        visited: Set[Node] = set()
        self.components: List[Dict[Node, Dict[Node, int]]] = []
        for node in self.heads:
            if node not in visited:
                sub_graph, explored_path = self.get_sub_graph(self.graph, node)
                if explored_path:
                    visited.update(explored_path)
                    self.components.append(sub_graph)
        return self.components

    def largest_component(self) -> AdjacencyList:
        """Returns a strongly connected components with largest number of nodes."""
        return max(self.components, key=lambda x: len(x))

    def to_adjacency_matrix(self) -> np.ndarray:
        """Returns graph as an adejecency matrix representation."""
        adjacency_matrix: np.ndarray = np.zeros([self.node_count, self.node_count])
        nodes = list(self.nodes.values())
        for node in nodes:
            neighbours = self.graph.get(node, [])
            for neighbour in neighbours:
                adjacency_matrix[nodes.index(node), nodes.index(neighbour)] += 1
        return adjacency_matrix

    def branching_nodes(self) -> List[Node]:
        """Returns a list of nodes with more than one unique neighbours."""
        branching_nodes: List[Node] = []
        for node, neighbours in self.graph.items():
            if len(set(neighbours)) > 1:
                branching_nodes.append(node)
        return branching_nodes

    def _travel_to_tip(self, _from: Node) -> List[Node]:
        """Helper function. Returns a path from node _from to closest tip."""
        discovered = []

        def _travel(node: Node):
            discovered.append(node)
            if node.is_tip():
                return
            for neighbour in self.graph[node]:
                if neighbour not in discovered:
                    _travel(neighbour)

        _travel(_from)
        return discovered

    def _travel_to_branch(self, _from: Node) -> List[Node]:
        """Helper funciton. Returns a path from node _from to closest branching node."""
        discovered = []

        def _travel(node: Node):
            discovered.append(node)
            if len(self.graph[node]) > 1:
                return
            for neighbour in self.graph[node]:
                if neighbour not in discovered:
                    _travel(neighbour)

        _travel(_from)
        return discovered

    def convergent_nodes(self) -> List[Node]:
        """Returns a list of nodes with more than 1 unique predecessors."""
        counts: Dict[Node, int] = {}
        for neighbours in self.graph.values():
            for node in neighbours:
                counts.setdefault(node, 0)
                counts[node] += 1
        self.convergent = [node for node, count in counts.items() if count > 1]
        return self.convergent

    def _travel_to_convergent(self, _from: Node) -> List[Node]:
        """Helper function. Returns path from node _from to closest convergent node."""
        discovered = []

        def _travel(node: Node):
            discovered.append(node)
            if len(node.predecessors) > 1:
                return
            for neighbour in self.graph[node]:
                if neighbour not in discovered:
                    _travel(neighbour)

        _travel(_from)
        return discovered

    def tips_lengths(self) -> List[int]:
        """Returns list of lengths of dead-end paths in a graph"""
        lengths: List[int] = []
        for node in self.branching_nodes():
            neighbours = self.graph.get(node, {})
            trails = []
            for neighbour in neighbours:
                trails.append(self._travel_to_tip(neighbour))
            lengths.append(len(min(trails, key=lambda x: len(x))))
        for node in self.heads:
            neighbours = self.graph.get(node, {})
            trails = []
            for neighbour in neighbours:
                trails.append(self._travel_to_convergent(neighbour))
            lengths.append(len(min(trails, key=lambda x: len(x))))
        return lengths

    def get_shortest_tips(self) -> List[List[Node]]:
        """Returns a list of shortest tips in a graph"""
        shortest_tips: List[List[Node]] = []
        for node in self.branching_nodes():
            neighbours = self.graph.get(node, {})
            trails = []
            for neighbour in neighbours:
                trails.append(self._travel_to_tip(neighbour))
            shortest_tips.append(min(trails, key=lambda x: len(x)))
        for node in self.heads:
            neighbours = self.graph.get(node, {})
            trails = []
            for neighbour in neighbours:
                trails.append(self._travel_to_branch(neighbour))
            shortest_tips.append(min(trails, key=lambda x: len(x)))
        return shortest_tips

    def trim_short_tips(self, threshold: int):
        """Returns a new graph, with trimmed tips of lenghts below thresholdself.

        This function does so iteratively. First removing the shortest_tips and then longer
        ones, until it reaches the threshold length."""
        shortest_tips = self.get_shortest_tips()
        k_mers = set()
        for i in range(threshold):
            for tip in [tip for tip in shortest_tips if len(tip) <= i]:
                for node in tip:
                    neighbours = self.graph[node]
                    for neighbour in neighbours:
                        k_mers.add(node.k_minus_1_mer + neighbour.k_minus_1_mer[-1])
        new_kmer_hist = {k_mer: count for k_mer, count in self.kmer_hist.items() if k_mer not in k_mers}
        return DeBruijn.from_kmer_hist(new_kmer_hist)

    def get_bubbles(self) -> List[List[Node]]:
        """Returns a list of paths in a graph that are part of 'bubbles'."""
        weaker_bubbles: List[List[Node]] = []
        for node in self.branching_nodes():
            neighbours = self.graph.get(node, {})
            weaker_bubbles.append(self._travel_to_convergent(min(neighbours, key=lambda x: neighbours[x])))
        return weaker_bubbles
    
    def trim_bubbles(self):
        """Returns a new graph with bubbles removed."""
        bubbles = self.get_bubbles()
        k_mers = set()
        for bubble in bubbles:
            for node in bubble:
                neighbours = self.graph[node]
                for neighbour in neighbours:
                    k_mers.add(node.k_minus_1_mer + neighbour.k_minus_1_mer[-1])
        new_kmer_hist = {k_mer: count for k_mer, count in self.kmer_hist.items() if k_mer not in k_mers}
        return DeBruijn.from_kmer_hist(new_kmer_hist)

    def merge_nodes(self, k: int, max_length) -> List[str]:
        """Merges graph nodes into contings using greedy superstring algorithm."""
        contigs = []
        for head in self.heads:
            discovered = []
            def _travel(node: Node):
                discovered.append(node)
                if node.outdegree == 0:
                    return
                neighbour = max(self.graph[node], key=lambda x: self.graph[node][x], default=None)
                if neighbour not in discovered and neighbour is not None:
                    _travel(neighbour)
            _travel(head)
            k_mers = [node.k_minus_1_mer for node in discovered]
            # contig = k_mers[0]
            # for kmer in k_mers[1:]:
            #     if scs.overlap(contig, kmer):
            #         contig += kmer[-1]
            # contigs.append(contig)
            contigs.append(scs.greedy_scs(k_mers, k, max_length))
        return contigs


def main() -> int:
    """Main function of this module. USED ONLY FOR DEBUGGIN PURPOSES!"""
    data = parse_fasta("../training/reads/reads1.fasta")
    k_mers = utils.draw_k_mers(data, 17)
    data = error_correction.correct_dataset(data, 17, k_mers)
    # data = data[0 : len(data) // 2]
    # data = utils.reads_from_sequence(utils.random_sequence(16000, "actg"), 80, 100)
    # data = ["a_long_", "g_long_time_ago", "long_gttime_ago"]
    # print(len(data))
    supergraph = DeBruijn.from_reads(data, 17)
    print(f"nodes on the supergraph: {supergraph.node_count}")
    print(f"edes on the supergraph: {supergraph.edge_count}")
    print(f"balanced: {supergraph.balanced}")
    print(f"semi-balanced: {supergraph.semi_balanced}")
    print(f"unbalanced: {supergraph.unbalanced}")
    print(f"number of heads: {len(supergraph.heads)}")
    print(f"number of tails: {len(supergraph.tails)}")
    print(f"number of tips: {len(supergraph.tips)}")
    print(f"branching nodes: {len(supergraph.branching_nodes())}")
    print(f"Convergent nodes: {len(supergraph.convergent_nodes())}")
    print(f"lengths of tips in supergraph: {supergraph.tips_lengths()}")
    print(f"Lenghts of bubbles in supergraph: {[len(bubble) for bubble in supergraph.get_bubbles()]}")
    # supergraph.show("graph", weights=true)
    supergraph = supergraph.trim_short_tips(17)
    supergraph = supergraph.trim_bubbles()
    components = [
        DeBruijn.from_graph(subgraph) for subgraph in supergraph.connected_components()
    ]
    # components = [component for component in components if component.node_count > 70]
    print(f"number of connected components/subgraphs: {len(components)}")
    print(
        f"lenghts of components: {[component.node_count for component in components]}"
    )
    print(
        f"sum of lenghts of components: {sum(component.node_count for component in components)}"
    )
    print(
        f"tips in components: {[component.tips_lengths() for component in components]}"
    )
    print(
        f"bubbles in components: {[len(component.get_bubbles()) for component in components]}"
    )
    print(
        f"branches in components: {[len(component.branching_nodes()) for component in components]}"
    )
    print(f"number of tips: {[len(component.tips) for component in components]}")
    # for component in components:
    #     # if any(0 < length < 18 for length in component.tips_lengths()):
    #     # component.trim_short_tips()
    #     print(component.branching_nodes())
    #     # component.show("graph.pdf", weights=True)
    #     # for k_1m_mer, node in component.nodes.items():
    #     #     print(
    #     #         f"k-1 mer: {k_1m_mer}, in: {node.indegree}, out: {node.outdegree}, tip: {node.is_tip()}"
    #     #     )
    #     # input(">>>")

    for component in components:
        component._update_tips()

    print(f"Number of tips: {[len(component.tips) for component in components]}")

    contigs = []
    for i, component in enumerate(components):
        print(f"Merging contig {i}")
        contigs.append(component.merge_nodes(13, 15))
    for contig in contigs:
        print(f"Number of contigs: {len(contig)}, Contigs' lenght: {[len(c) for c in contig]}")
    with open("./output_contings.fasta", "w") as output:
        for i, contig in enumerate(contigs):
            output.write(f">contig_{i}\n")
            output.write(f"{contig[0]}\n")
    return 0

# Run main only if this module is run directly and not imported
if __name__ == "__main__":
    sys.exit(main())
