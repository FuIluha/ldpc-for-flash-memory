import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from common.LinearCodes import LinearCodes
import numpy as np

class ChanelWithErrors(LinearCodes):
    @staticmethod
    def get_messages(c: np.array, ps: np.array, n: int = 100) -> np.array:
        e_list = np.random.choice(2, (n, len(c)), p=[1. - ps[0], ps[0]])
        e_list = np.array([
            np.random.choice(2, size=(n, len(c)), p=[1 - p, p]) for p in ps
        ])
        c_list = np.tile(c, (ps.shape[0], n, 1))
        r_list = e_list ^ c_list
        return r_list
    
