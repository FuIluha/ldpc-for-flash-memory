import numpy as np
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np
from common.LinearCodes import LinearCodes

class GallagerMatrixGenerator(LinearCodes):
    @staticmethod
    def generate(m, n, l):
        if m % l != 0:
            raise ValueError("m должно быть кратно l")
        n0 = (n * l) // m
        print(n0)
        
        layers = l
        layer_height = m // layers
        base_layer = np.zeros((layer_height, n), dtype=int)
        for i in range(layer_height):
            base_layer[i, i * n0:(i + 1) * n0] = 1
        gallager_matrix = np.vstack([base_layer] * layers)
        for layer in range(1, layers):
            permuted_columns = np.random.permutation(n)
            new_layer = base_layer[:, permuted_columns]
            gallager_matrix[layer * layer_height:(layer + 1) * layer_height, :] = new_layer

        return gallager_matrix

# n = 1000
# m = 1500
# l = 6
# n0 = 4
# matrix = create_gallager_matrix(n, m, l, n0)
# text_output_file = "gallager_matrix.txt"
# np.savetxt(text_output_file, matrix, fmt="%d")