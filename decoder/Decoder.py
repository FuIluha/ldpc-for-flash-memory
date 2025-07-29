import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from common.LinearCodes import LinearCodes
from coder.ChanelWithErrors import ChanelWithErrors
from concurrent.futures import ProcessPoolExecutor
import numpy as np

class Decoder(LinearCodes):
    @staticmethod
    def __bitflipping_decoding(r, H, equation_matrix, bit_in_sindrom, max_deep, deep=0):
        if deep >= max_deep:
            return r

        S = Decoder._dot_for_codes(H, r.T)

        if not np.any(S):
            return r
        
        strange_bits = {}
        for i in range(len(S)):
            if S[i] == 1:
                for bit in equation_matrix[i]:
                    strange_bits[bit] = strange_bits.get(bit, 0) + 1
        strangest_bit = min(strange_bits.items(), key=lambda x: bit_in_sindrom[x[0]] - x[1])[0]
        r[strangest_bit] = (r[strangest_bit] + 1) % 2
        return Decoder.__bitflipping_decoding(r, H, equation_matrix, bit_in_sindrom, max_deep, deep=deep + 1)

    @staticmethod
    def _decode_single_element(args):
        i, j, H, equation_matrix, bit_in_sindrom, deep, c, r_ij = args
        decoded = Decoder.__bitflipping_decoding(r_ij, H, equation_matrix, bit_in_sindrom, max_deep=deep)
        error = not np.array_equal(c, decoded)
        return i, j, error, decoded

    @staticmethod
    def bitflipping_decoding_for_many(H: np.array, c: np.array, p: np.array, repeat: int, deep: int) -> np.array:
        n = repeat
        r = ChanelWithErrors.get_messages(c, p, n=n)

        equation_matrix = np.array([np.where(H[i] == 1)[0] for i in range(H.shape[0])], dtype=object)
        bit_in_sindrom = np.sum(H, axis=0)

        count_errors = np.zeros(r.shape[0])
        
        tasks = [
            (i, j, H, equation_matrix, bit_in_sindrom, deep, c, r[i, j])
            for i in range(r.shape[0])
            for j in range(r.shape[1])
        ]

        decoded_r = np.empty_like(r)

        with ProcessPoolExecutor() as executor:
            results = executor.map(Decoder._decode_single_element, tasks)

        for i, j, error, decoded in results:
            if error:
                count_errors[i] += 1
            decoded_r[i, j] = decoded

        return count_errors / n