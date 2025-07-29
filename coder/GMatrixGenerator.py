import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from common.LinearCodes import LinearCodes
import numpy as np

class GMatrixGenerator(LinearCodes):
    def __init__(self, H):
        self.H = H
        self.operations = []

    def __to_sistematic(self):
        sistematic_H = self.H.copy()
        rows, cols = sistematic_H.shape
        n = cols
        k = cols - rows

        for i in range(n - k):
            diag_col = i + k
            if sistematic_H[i, diag_col] != 1:
                for j in range(n):
                    if sistematic_H[i, j] == 1:
                        self.operations.append(np.array((diag_col, j)))
                        sistematic_H[:, diag_col] ^= sistematic_H[:, j]
                        break

            for j_k in range(n):
                if sistematic_H[i, j_k] == 1 and j_k != diag_col:
                    self.operations.append((j_k, diag_col))
                    sistematic_H[:, j_k] ^= sistematic_H[:, diag_col]

    def save_matrix(self):
        rows, cols = self.H.shape
        n = cols
        k = cols - rows

        self.__to_sistematic()
        sistematic_G = np.hstack((np.eye(k, dtype=int), np.zeros((k, n - k), dtype=int)))
        for operation in self.operations[::-1]:
            first, second = operation[0], operation[1]
            sistematic_G[:, second] ^= sistematic_G[:, first]
        np.savetxt("data/g_matrix.txt", sistematic_G, delimiter=',', fmt='%d')
    
    def get_matrix(self) -> np.array:
        rows, cols = self.H.shape
        n = cols
        k = cols - rows

        self.__to_sistematic()
        sistematic_G = np.hstack((np.eye(k, dtype=int), np.zeros((k, n - k), dtype=int)))
        for operation in self.operations[::-1]:
            first, second = operation[0], operation[1]
            sistematic_G[:, second] ^= sistematic_G[:, first]
        return sistematic_G
