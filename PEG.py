class Graph:
    def __init__(self):
        self.nodes = {}
        self.edges = []

    def add_node(self, node, node_type):
        self.nodes[node] = {'node_type': node_type}

    def add_edge(self, node1, node2):
        if not self.has_edge(node1, node2):
            self.edges.append((node1, node2))

    def neighbors(self, node):
        neighbors = []
        for edge in self.edges:
            if edge[0] == node:
                neighbors.append(edge[1])
            elif edge[1] == node:
                neighbors.append(edge[0])
        return neighbors

    def has_edge(self, node1, node2):
        return (node1, node2) in self.edges or (node2, node1) in self.edges

def create_graph(n, m, row_weight):
    graph = Graph()
    s_nodes = ['s{}'.format(i) for i in range(n)]
    for s_i in s_nodes:
        graph.add_node(s_i, 's')
    c_nodes = ['c{}'.format(i) for i in range(m)]
    for c_i in c_nodes:
        graph.add_node(c_i, 'c')

    for s_i in s_nodes:
        available_c_nodes = [c for c in c_nodes if not graph.has_edge(s_i, c)]
        while len(graph.neighbors(s_i)) < row_weight:
            min_degree_c_nodes = [c for c in available_c_nodes if len(graph.neighbors(c)) == min([len(graph.neighbors(c1)) for c1 in available_c_nodes])]
            c_node = min_degree_c_nodes[0]
            graph.add_edge(s_i, c_node)
            available_c_nodes.remove(c_node)

    return graph

def create_parity_matrix(graph, n, m):
    parity_matrix = [[0] * n for _ in range(m)]
    for s_i in range(n):
        for c_i in range(m):
            if graph.has_edge('s{}'.format(s_i), 'c{}'.format(c_i)):
                parity_matrix[c_i][s_i] = 1
    return parity_matrix


# Example usage:
# n = 5
# m = 4
# row_weight = 3
# graph = create_graph(n, m, row_weight)
# parity_matrix = create_parity_matrix(graph, n, m)
#
# for row in parity_matrix:
#     print(row)


