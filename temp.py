
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
