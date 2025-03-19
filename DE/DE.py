import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
from scipy.special import comb
from coder.ChanelWithErrors import ChanelWithErrors
import matplotlib.pyplot as plt

def density_evolution(H: np.ndarray):
    l = np.sum(H[:, 0])
    n0 = np.sum(H[0, :])
    t = l / 2

    p = np.linspace(0, 0.35, 2000)

    pcn = (1 - np.power((1 - 2 * p), 2)) / 2
    
    i_values = np.arange(t, l + 1)
    binom_coeffs = comb(l, i_values)
    p_cn_powers = np.power(pcn[:, np.newaxis], i_values)
    one_minus_p_cn_powers = np.power(1 - pcn[:, np.newaxis], l - i_values)

    pvn = np.sum(binom_coeffs * p_cn_powers * one_minus_p_cn_powers, axis=1)

    f = pvn - p

    plt.plot(p, f, label=f'$n_0 = {n0:.0f}$, $l = {l:.0f}$')
