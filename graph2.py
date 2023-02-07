import sys
import copy
from typing import List, Dict, Optional, Set, Tuple, Union, TypeVar, Type

import graphviz
import numpy as np

import utils
import progress
import error_correction
import greedy_scs
from analyze_training_data import parse_fasta

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
        # self.distance: Union[float, int]
        self.predecessor: Optional[Node]

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


AdjacencyList = Dict[Node, List[Node]]


def get_sub_graph(
    super_graph: AdjacencyList, node: Node
) -> Tuple[AdjacencyList, List[Node]]:
    """Returns a connected components containing the node from a super-graph"""
    discovered: List[Node] = []
    sub_graph: AdjacencyList = {}

    def _travel(node: Node):
        discovered.append(node)
        sub_graph.setdefault(node, []).extend(super_graph.get(node, []))
        for neighbour in super_graph[node]:
            if neighbour not in discovered:
                _travel(neighbour)

    _travel(node)
    return sub_graph, discovered


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
                if node.indegree == node.outdegree + 1:
                    self.heads.add(node)
                if node.outdegree == node.indegree + 1:
                    self.tails.add(node)
            else:
                self.unbalanced += 1

    def _update_tips(self) -> None:
        for node in self.nodes.values():
            if node.is_tip():
                self.tips.add(node)

    def _to_dot(self, weights: bool = False) -> graphviz.Digraph:
        graph = graphviz.Digraph(comment="DeBruijn graph")
        for node in self.graph:
            graph.node(node.k_minus_1_mer, node.k_minus_1_mer)
        for source, destinations in self.graph.items():
            if weights:
                weightmap = {}
                for destination in destinations:
                    weightmap[destination] = weightmap.get(destination, 0) + 1
                for destination, weight in weightmap.items():
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

    @classmethod
    def from_reads(cls: Type[S], reads: List[str], k: int) -> S:
        """Create a de Bruijn Graph from list of strings."""
        new_graph = cls()
        for read in reads:
            for k_mer, _ in utils.chop(read, k):
                k_mer_L, k_mer_R = k_mer[:-1], k_mer[1:]
                nodeL = new_graph._update_nodes(k_mer_L)
                nodeR = new_graph._update_nodes(k_mer_R)
                new_graph.graph.setdefault(nodeL, []).append(nodeR)
                new_graph.graph.setdefault(nodeR, [])
                nodeL.outdegree += 1
                nodeR.indegree += 1
        new_graph._balance_check()
        new_graph._update_tips()
        return new_graph

    @classmethod
    def from_graph(cls: Type[S], graph: AdjacencyList) -> S:
        """Creates a de Bruijn graph from weighted adjacency list representation."""
        new_graph = cls()
        for node in graph:
            new_graph.nodes[node.k_minus_1_mer] = node.copy()
        new_graph.graph = graph.copy()
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

    # Raises NotEulerianGraphError
    def eulerian_walk_or_cycle(self) -> List[str]:
        """Returns sequences of nodes corresponding to eulerian walk or cycle on the graph.

        Returned nodes are represented by their k-1mersself.
        Raises NotEulerianGraphError if graph is not eulerian."""
        # if not self.is_eulerian():
        #     raise NotEulerianGraphError
        g = self.graph
        if walk := self.has_eulerian_walk():
            g = g.copy()
            # If graph has eulerian walk there should be exactly one head nad exactly one tail
            # That's why we should be able to simply get their values from heads and tails sets
            head, tail = self.heads.pop(), self.tails.pop()
            g.setdefault(tail, []).append(head)

        tour: List[Node] = []

        def _visit(node: Node) -> None:
            while len(g[node]) > 0:
                destination = g[node].pop()
                _visit(destination)
            tour.append(node)

        source = next(iter(g.keys()))
        _visit(source)
        tour.reverse()
        # tour.pop()
        try:
            head_index = tour.index(head)  # type: ignore
            tour = tour[head_index:] + tour[:head_index]
        except (NameError, ValueError):
            pass
        return list(map(str, tour))

    def show(self, filename: str, weights: bool = False) -> None:
        self._to_dot(weights=weights).view(filename)

    def connected_components(self) -> List[AdjacencyList]:
        visited: Set[Node] = set()
        self.components: List[AdjacencyList] = []
        for node in self.graph:
            if node not in visited:
                sub_graph, explored_path = get_sub_graph(self.graph, node)
                if explored_path:
                    visited.update(explored_path)
                    self.components.append(sub_graph)
        return self.components

    def largest_component(self) -> AdjacencyList:
        return max(self.components, key=lambda x: len(x))

    def to_adjacency_matrix(self) -> np.ndarray:
        adjacency_matrix: np.ndarray = np.zeros([self.node_count, self.node_count])
        nodes = list(self.nodes.values())
        for node in nodes:
            neighbours = self.graph.get(node, [])
            for neighbour in neighbours:
                adjacency_matrix[nodes.index(node), nodes.index(neighbour)] += 1
        return adjacency_matrix

    def branching_nodes(self) -> int:
        branching: int = 0
        for neighbours in self.graph.values():
            if len(set(neighbours)) > 1:
                branching += 1
        return branching


def main() -> int:
    # data = parse_fasta("../training/reads/reads1.fasta")
    # k_mers = utils.draw_k_mers(data, 17)
    # data = error_correction.correct_dataset(data, 17, k_mers)
    # data = data[0: len(data) // 5]
    # data = utils.reads_from_sequence(utils.random_sequence(16000, "ACTG"), 80, 100)
    data = ["A_long_long_long_", "g_long_time_ago"]
    # print(len(data))
    supergraph = DeBruijn.from_reads(data, 5)
    print(f"Nodes on the supergraph: {supergraph.node_count}")
    print(f"Edes on the supergraph: {supergraph.edge_count}")
    print(f"Balanced: {supergraph.balanced}")
    print(f"Semi-balanced: {supergraph.semi_balanced}")
    print(f"Unbalanced: {supergraph.unbalanced}")
    print(f"Number of Heads: {len(supergraph.heads)}")
    print(f"Number of tails: {len(supergraph.tails)}")
    print(f"Number of tips: {len(supergraph.tips)}")
    print(f"Branching nodes: {supergraph.branching_nodes()}")
    components = [
        DeBruijn.from_graph(subgraph) for subgraph in supergraph.connected_components()
    ]
    # components = [component for component in components if component.node_count > 70]
    print(f"Number of connected components/subgraphs: {len(components)}")
    print(
        f"Lenghts of components: {[component.node_count for component in components]}"
    )
    components[0].show("graph")
    # with open("./output_contings.fasta", "w") as output:
    #     for i, component in enumerate(components):
    #         # walk = component.eulerian_walk_or_cycle()
    #         # walk = walk[0] + "".join(map(lambda x: x[-1], walk[1:]))
    #         string = greedy_scs.greedy_scs(component.nodes, 40)
    #         output.write(f">contig_{i}\n")
    #         output.write(f"{string}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
