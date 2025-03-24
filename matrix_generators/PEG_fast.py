import numpy as np
import heapq

class PEG_fast:
    def __init__(self, m, n, l):
        self.n = n
        self.m = m
        self.l = l
        self.H = np.zeros((m, n), dtype=np.int8)
        self.D_s = np.zeros(n, dtype=np.int16)
        self.D_c = np.zeros(m, dtype=np.int16)
        self.D_s_future = np.zeros(n, dtype=np.int16) + l

    def find_min_degree_node(self, j):
        available_nodes = np.where(self.H[:, j] == 0)[0]
        if available_nodes.size == 0:
            return np.argmin(self.D_c)
        degrees = self.D_c[available_nodes]
        return available_nodes[np.argmin(degrees)]

    def peg(self):
        for j in range(self.n):
            for k in range(self.D_s_future[j]):
                if k == 0:
                    c = np.argmin(self.D_c)
                else:
                    c = self.find_min_degree_node(j)
                self.H[c, j] = 1
                self.D_s[j] += 1
                self.D_c[c] += 1

            np.savetxt(f'data/H_matrix/H{self.l}_{self.m}x{self.n}_fastPEG.txt', self.H, fmt='%d', delimiter=',')

# n = 1500
# m = 1000
# for l in range(5, 9):
#     PEG_fast(m, n, l).peg()