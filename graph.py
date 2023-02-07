import sys
from typing import Callable, Iterable, List, Dict, Optional, Set, Tuple, Union

import graphviz

import utils
import progress
import error_correction
from analyze_training_data import parse_fasta

import greedy_scs

sys.setrecursionlimit(10000)


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
        self.distance: Union[float, int]
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


class DeBruijn:
    """Represents a de Beuijn graph constructed from k-mers"""

    def _update_nodes(self, k_mins_1_mer: str) -> Node:
        if k_mins_1_mer not in self.nodes:
            self.nodes[k_mins_1_mer] = Node(k_mins_1_mer)
        return self.nodes[k_mins_1_mer]

    def _set_head_and_tail(self) -> None:
        for node in self.nodes.values():
            if node.is_balanced():
                self.balenced += 1
            elif node.is_semi_balanced():
                if node.indegree == node.outdegree + 1:
                    self.head = node
                if node.outdegree == node.indegree + 1:
                    self.tail = node
                self.semi_balanced += 1
            else:
                self.unbalanced += 1


    def __init__(self) -> None:
        self.nodes: Dict[str, Node] = {}
        self.graph: Dict[Node, List[Node]] = {}
        self.balenced: int = 0
        self.semi_balanced: int = 0
        self.unbalanced: int = 0
        self.head: Node
        self.tail: Node

    @classmethod
    def from_reads(cls, strings: Iterable[str], k: int):
        new_graph = cls()
        for string in strings:
            for k_mer, _ in utils.chop(string, k):
                k_mer_L, k_mer_R = k_mer[:-1], k_mer[1:]
                nodeL = new_graph._update_nodes(k_mer_L)
                nodeR = new_graph._update_nodes(k_mer_R)
                new_graph.graph.setdefault(nodeL, []).append(nodeR)
                new_graph.graph.setdefault(nodeR, [])
                nodeL.outdegree += 1
                nodeR.indegree += 1
        new_graph._set_head_and_tail()
        new_graph.graph.setdefault(new_graph.tail, [])
        return new_graph

    @classmethod
    def from_graph(cls, graph: Dict[Node, List[Node]]):
        new_graph = cls()
        for node in graph:
            new_graph.nodes[node.k_minus_1_mer] = node
        new_graph.graph = graph
        new_graph._set_head_and_tail()
        return new_graph

    def heads(self) -> List[Node]:
        possible_heads: List[Node] = []
        for node in self.nodes.values():
            if node.outdegree > node.indegree and node.indegree == 0:
                possible_heads.append(node)
        return possible_heads

    def best_head(self) -> Optional[Node]:
        for node in self.nodes.values():
            if node.outdegree > node.indegree and node.outdegree == 0:
                return node
        return self.heads()[0]

    def tails(self) -> List[Node]:
        possible_tails: List[Node] = []
        for node in self.nodes.values():
            if node.indegree > node.outdegree and node.outdegree == 0: 
                possible_tails.append(node)
        return possible_tails

    @property
    def node_count(self) -> int:
        """Number of nodes in a graph."""
        return len(self.nodes)

    @property
    def edge_count(self) -> int:
        """Number of egdes in a graph"""
        all_egdes: List[Node] = []
        for values in self.graph.values():
            all_egdes.extend(values)
        return len(all_egdes)


    def has_eulerian_walk(self) -> bool:
        """Returns True if graph has eulerian walk."""
        return self.unbalanced == 0 and self.semi_balanced == 2

    def has_eulerian_cycle(self) -> bool:
        """Returns True if graph has eulerian cycle."""
        return self.unbalanced == 0 and self.semi_balanced == 0

    def is_eulerian(self) -> bool:
        """Returns True if graph has either eulerian walk or cycle."""
        return self.has_eulerian_walk() or self.has_eulerian_cycle()

    def eulerian_walk_or_cycle(self) -> list[str]:
        """Returns sequence of nodes correspond to eulerian walk or cycle.

        Returned nodes are represented by their k-1 mers."""
        if self.is_eulerian is False:
            raise NotEulerianGraphError
        g = self.graph
        if walk := self.has_eulerian_walk():
            g = g.copy()
            g.setdefault(self.tail, []).append(self.head)

        tour: list[Node] = []

        def _visit(node: Node) -> None:
            while len(g[node]) > 0:
                destination = g[node].pop()
                _visit(destination)
            tour.append(node)

        source = next(iter(g.keys()))
        _visit(source)
        tour.reverse()
        tour.pop()

        if walk:
            head_index = tour.index(self.head)
            tour = tour[head_index:] + tour[:head_index]

        return list(map(str, tour))

    def to_dot(self, weights: bool = False) -> graphviz.Digraph:
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

    def show(self, filename: str, weights: bool = False) -> None:
        self.to_dot(weights=weights).view(filename)

    def connected_components(self) -> List[Dict[Node, List[Node]]]:
        bar = progress.ProgressBar(len(self.graph))
        count = 0

        def sub_graph(G: Dict[Node, List[Node]], node: Node):
            discovered: List[Node] = []
            sub_graph: Dict[Node, List[Node]] = {}

            def _travel(node: Node):
                discovered.append(node)
                sub_graph.setdefault(node, []).extend(G.get(node, []))
                for neighbour in G[node]:
                    if neighbour not in discovered:
                        _travel(neighbour)

            _travel(node)
            return sub_graph, discovered

        visited: Set[Node] = set()
        self.components: List[Dict[Node, List[Node]]] = []
        for node in self.graph:
            count += 1
            if node not in visited:
                component, explord_path = sub_graph(self.graph, node)
                if explord_path:
                    visited.update(explord_path)
                    self.components.append(component)
            bar.update(len(visited))
        bar.finish()
        return self.components

    def count_components(self) -> int:
        bar = progress.ProgressBar(len(self.graph))
        count = 0

        visited: Set[Node] = set()
        components: int = 0

        def _explore(graph, current):
            if current in visited:
                return False
            visited.add(current)
            for neighbour in graph[current]:
                _explore(graph, neighbour)
            return True

        for node in self.graph:
            count += 1
            if _explore(self.graph, node):
                components += 1
        return components

    def largest_component(self) -> Dict[Node, List[Node]]:
        return max(self.components, key=lambda x: len(x))

    def _init_single_source(self, start: Node) -> None:
        for node in self.nodes.values():
            node.distance = float("inf")
            node.predecessor = None
        start.distance = 0

    def _weight(self, _from: Node, _to: Node) -> int:
        neighbours = self.graph.get(_from, [])
        return neighbours.count(_to)

    def _relax(self, _from: Node, _to: Node, weight: Callable[[Node, Node], int]):
        if _to.distance > _from.distance + weight(_from, _to):
            _to.distance = _from.distance + weight(_from, _to)
            _to.predecessor = _from

    def bellman_ford(self, start) -> bool:
        self._init_single_source(start)
        for i in range(len(self.nodes) - 1):
            for node, neighbours in self.graph.items():
                for neighbour in neighbours:
                    self._relax(node, neighbour, self._weight)
        for node, neighbours in self.graph.items():
            for neighbour in neighbours:
                if neighbour.distance > node.distance + self._weight(node, neighbour):
                    return False
        return True


def main() -> int:
    # s1 = "to_every_thing_turn_turn_turn_there_is_a_season_turn_turn_turn"
    # s2 = "ABCDFED"
    data = parse_fasta("../training/reads/reads1.fasta")
    k_mers = utils.draw_k_mers(data, 45)
    data = error_correction.correct_dataset(data, 45, k_mers)
    data = data[: len(data) // 5]
    graph = DeBruijn.from_reads(data, k=45)
    print(graph.balenced)
    print(graph.semi_balanced)
    print(graph.unbalanced)
    print(len(graph.nodes))
    print(graph.edge_count)
    print(graph.has_eulerian_walk())
    print(graph.has_eulerian_cycle())
    components = graph.connected_components()
    graphs = [DeBruijn.from_graph(component) for component in components]
    print(len(graphs))
    print([g.node_count for g in graphs])
    # cyclic: int = 0
    # for g in graphs:
    #     try:
    #         g.bellman_ford(g.heads()[0])
    #     except IndexError:
    #         pass
    # 
    # walks = []
    # for tail in graphs[0].tails():
    #     node = tail
    #     walk = [node]
    #     while node.predecessor is not None:
    #         walk = [node.predecessor] + walk
    #         node = node.predecessor
    #     walks.append(walk)
    #
    # print(walks)
    # for node in graph.nodes.values():
    #     print(f"{node}: in: {node.indegree}, out: {node.outdegree}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
