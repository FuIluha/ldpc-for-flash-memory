import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from common.LinearCodes import LinearCodes
from coder.GMatrixGenerator import GMatrixGenerator
import numpy as np

class Coder(LinearCodes):
    @staticmethod
    def get_message(H: np.array) -> np.array:
        G = GMatrixGenerator(H).get_matrix()
        u = np.random.randint(0, 2, G.shape[0], dtype=np.int8)
        return Coder._dot_for_codes(u, G)
