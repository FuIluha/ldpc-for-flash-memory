import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from common.LinearCodes import LinearCodes
from coder.ChanelWithErrors import ChanelWithErrors
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import shared_memory
from functools import cache
from numba import jit, prange
import numpy as np

class Decoder(LinearCodes):
    @staticmethod
    def _decode_single_element(args):
        i, j, H, equation_matrix, deep, r_ij, with_threshold = args
        
        @cache
        def bf_decoder(deep, r: tuple, with_threshold, current_deep = 0) -> tuple:
            r = np.array(r, dtype=np.int8)
            S = (H @ r.T) % 2

            if current_deep >= deep:
                return tuple(r)

            if not np.any(S):
                return tuple(r)
            
            strange_bits = np.zeros(H.shape[1])
            error_S = np.argwhere(S == 1).flatten()
            
            indices = np.concatenate([equation_matrix[idx] for idx in error_S]).astype(int)
            np.add.at(strange_bits, indices, 1)

            if with_threshold:
                l = np.sum(H[:, 0])
                t = l / 2 + l % 2
                r[strange_bits >= t] ^= 1
            else:
                r[np.argmax(strange_bits)] ^= 1
            
            return bf_decoder(deep, tuple(r), with_threshold, current_deep + 1)

        r = bf_decoder(deep, tuple(r_ij), with_threshold)

        return i, j, np.any(r), np.sum(r)
    
    @staticmethod
    def _tasks_generator(H, equation_matrix, deep, r, with_threshold):
        for i in prange(r.shape[0]):
            for j in prange(r.shape[1]):
                yield i, j, H, equation_matrix, deep, r[i, j], with_threshold

    @staticmethod
    @jit(nopython=True)
    def _count_error(ber, results, r, c, n):
        count_errors = np.zeros(r.shape[0])

        if ber:
            for i, j, error, decoded in results:
                if error:
                    count_errors[i] += decoded
            
            return count_errors / (n * c.shape[0])
        else:
            for i, j, error, decoded in results:
                if error:
                    count_errors[i] += 1

            return count_errors / n
        

    @staticmethod
    def _create_shared_array(arr):
        shm = shared_memory.SharedMemory(create=True, size=arr.nbytes)
        shared_arr = np.ndarray(arr.shape, dtype=arr.dtype, buffer=shm.buf)
        shared_arr[:] = arr[:]
        return shm, shared_arr

    @staticmethod
    def bitflipping_decoding_for_many(H, c, p, repeat, deep, ber=False, with_threshold=False):
        n = repeat
        r = ChanelWithErrors.get_messages(c, p, n=n)
        equation_matrix = np.array([np.where(H[i] == 1)[0] for i in range(H.shape[0])], dtype=object)

        shm_H, shared_H = Decoder._create_shared_array(H)
        shm_eq, shared_eq = Decoder._create_shared_array(equation_matrix)
        shm_r, shared_r = Decoder._create_shared_array(r)

        with ProcessPoolExecutor() as executor:
            results = executor.map(Decoder._decode_single_element, 
                                Decoder._tasks_generator(shared_H, shared_eq, deep, shared_r, with_threshold))
        
        decoded_errors = Decoder._count_error(ber, list(results), shared_r, c, n)

        shm_H.close()
        shm_H.unlink()
        shm_eq.close()
        shm_eq.unlink()
        shm_r.close()
        shm_r.unlink()

        return decoded_errors
