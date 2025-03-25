import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from common.LinearCodes import LinearCodes
from coder.ChanelWithErrors import ChanelWithErrors
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import shared_memory
import numpy as np

class Decoder(LinearCodes):
    @staticmethod
    def __bitflipping_decoding(r, H, equation_matrix, max_deep, deep=0, with_threshold: bool=False):
        if deep >= max_deep:
            return r

        S = (H @ r.T) % 2

        if not np.any(S):
            return r
        
        strange_bits = np.zeros(H.shape[1])
        error_S = np.argwhere(S == 1)
        np.add.at(strange_bits, equation_matrix[error_S].ravel(), 1)

        if with_threshold:
            l = np.sum(H[:, 0])
            t = l / 2 + l % 2
            r[strange_bits >= t] ^= 1
        else:
            r[np.argmax(strange_bits)] ^= 1
        return Decoder.__bitflipping_decoding(r, H, equation_matrix, max_deep, deep=deep + 1, with_threshold=with_threshold)

    @staticmethod
    def _decode_single_element(args):
        i, j, H, equation_matrix, deep, c, r_ij, with_threshold = args
        decoded = Decoder.__bitflipping_decoding(r_ij, H, equation_matrix, max_deep=deep, with_threshold=with_threshold)
        error = not np.array_equal(c, decoded)
        return i, j, error, decoded

    @staticmethod
    def bitflipping_decoding_for_many(
        H: np.ndarray, 
        c: np.ndarray, 
        p: np.ndarray, 
        repeat: int, 
        deep: int, 
        ber: bool = False, 
        with_threshold: bool = False
    ) -> np.ndarray:
        n = repeat
        r = ChanelWithErrors.get_messages(c, p, n=n)

        equation_matrix = np.array([np.where(H[i] == 1)[0] for i in range(H.shape[0])], dtype=np.int16)
        
        tasks = [
            (i, j, H, equation_matrix, deep, c, r[i, j], with_threshold)
            for i in range(r.shape[0])
            for j in range(r.shape[1])
        ]

        count_errors = np.zeros(r.shape[0])

        with ProcessPoolExecutor() as executor:
            results = executor.map(Decoder._decode_single_element, tasks)

        if ber:
            for i, j, error, decoded in results:
                if error:
                    count_errors[i] += np.sum(decoded)
            
            return count_errors / (n * c.shape[0])
        else:
            for i, j, error, decoded in results:
                if error:
                    count_errors[i] += 1

            return count_errors / n
