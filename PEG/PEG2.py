# class Graph:
#     def __init__(self):
#         self.nodes = {}
#         self.edges = []
#
#     def add_node(self, node, node_type):
#         self.nodes[node] = {'node_type': node_type}
#
#     def add_edge(self, node1, node2):
#         if not self.has_edge(node1, node2):
#             self.edges.append((node1, node2))
#
#     def neighbors(self, node):
#         neighbors = []
#         for edge in self.edges:
#             if edge[0] == node:
#                 neighbors.append(edge[1])
#             elif edge[1] == node:
#                 neighbors.append(edge[0])
#         return neighbors
#
#     def has_edge(self, node1, node2):
#         return (node1, node2) in self.edges or (node2, node1) in self.edges
#
# def create_graph(n, m, row_weight):
#     graph = Graph()
#     s_nodes = ['s{}'.format(i) for i in range(n)]
#     for s_i in s_nodes:
#         graph.add_node(s_i, 's')
#     c_nodes = ['c{}'.format(i) for i in range(m)]
#     for c_i in c_nodes:
#         graph.add_node(c_i, 'c')
#
#     for s_i in s_nodes:
#         available_c_nodes = [c for c in c_nodes if not graph.has_edge(s_i, c)]
#         while len(graph.neighbors(s_i)) < row_weight:
#             min_degree_c_nodes = [c for c in available_c_nodes if len(graph.neighbors(c)) == min([len(graph.neighbors(c1)) for c1 in available_c_nodes])]
#             c_node = min_degree_c_nodes[0]
#             graph.add_edge(s_i, c_node)
#             available_c_nodes.remove(c_node)
#
#     return graph
#
# def create_parity_matrix(graph, n, m):
#     parity_matrix = [[0] * n for _ in range(m)]
#     for s_i in range(n):
#         for c_i in range(m):
#             if graph.has_edge('s{}'.format(s_i), 'c{}'.format(c_i)):
#                 parity_matrix[c_i][s_i] = 1
#     return parity_matrix


# Example usage:
# n = 1000
# m = 1500
# row_weight = 9
# graph = create_graph(n, m, row_weight)
# parity_matrix = create_parity_matrix(graph, n, m)
#
# for row in parity_matrix:
#     print(row)


import numpy as np

class Graph:
    def __init__(self, n, m):
        self.n = n
        self.m = m
        self.adj_matrix = np.zeros((n, m), dtype=int)

    def add_edge(self, node1, node2):
        self.adj_matrix[node1, node2] = 1

    def neighbors(self, node):
        if node < self.n:
            return np.where(self.adj_matrix[node, :] == 1)[0]
        else:
            return np.where(self.adj_matrix[:, node - self.n] == 1)[0]

    def has_edge(self, node1, node2):
        return self.adj_matrix[node1, node2] == 1

def create_graph(n, m, row_weight):
    graph = Graph(n, m)

    for s_i in range(n):
        available_c_nodes = np.where(graph.adj_matrix[s_i, :] == 0)[0]
        while len(graph.neighbors(s_i)) < row_weight:
            if len(available_c_nodes) == 0:
                break
            min_degree_c_nodes = available_c_nodes[np.argmin([len(graph.neighbors(c)) for c in available_c_nodes])]
            graph.add_edge(s_i, min_degree_c_nodes)
            available_c_nodes = np.delete(available_c_nodes, np.where(available_c_nodes == min_degree_c_nodes))

    return graph

def create_parity_matrix(graph):
    return graph.adj_matrix

# Example usage:
n = 4
m = 5
row_weight = 3

graph = create_graph(n, m, row_weight)
parity_matrix = create_parity_matrix(graph)
print(parity_matrix)
