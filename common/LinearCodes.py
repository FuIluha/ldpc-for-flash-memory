import numpy as np

class LinearCodes:
    @staticmethod
    def _dot_for_codes(a: np.array, b: np.array) -> np.array:
        return np.array(a @ b % 2)