import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
from scipy.special import comb
from coder.ChanelWithErrors import ChanelWithErrors
import matplotlib.pyplot as plt

def density_evolution_f(H: np.ndarray):
    l = np.sum(H[:, 0])
    n0 = np.sum(H[0, :])

    R = np.unique(np.sum(H, axis=1), return_counts=True)
    degres = R[0] * R[1]

    t = int(l // 2 + l % 2)

    p = np.linspace(0, 0.1, 2000)

    rho = np.sum(degres * np.power((1 - 2 * p[:, np.newaxis]), R[0] - 1), axis=1) / np.sum(degres)

    pcn = (1 - rho) / 2
    
    i_values = np.arange(t, l + 1)
    binom_coeffs = comb(l, i_values)
    p_cn_powers = np.power(pcn[:, np.newaxis], i_values)
    one_minus_p_cn_powers = np.power(1 - pcn[:, np.newaxis], l - i_values)

    pvn = np.sum(binom_coeffs * p_cn_powers * one_minus_p_cn_powers, axis=1)

    f = pvn - p

    rho_x = ''
    for i in range(len(degres)):
        rho_x += f'\\frac{{{degres[i]:.0f}}}{{{np.sum(degres):.0f}}}x^{{{R[0][i] - 1:.0f}}}+'

    plt.plot(p, f, label=f'$\\rho(x) = {rho_x[:-1]}$, $l = {l:.0f}$, $t = {t}$')

def density_evolution(H: np.ndarray):
    l = np.sum(H[:, 0])
    R = np.unique(np.sum(H, axis=1), return_counts=True)
    degres = R[0] * R[1]

    t = int(l // 2 + l % 2)

    for p in np.linspace(0, 0.02, 6):
        last_pvn = p

        f = []
        all_pvn = []

        for _ in range(100):
            rho = np.sum((degres * np.power((1 - 2 * last_pvn), R[0] - 1))) / np.sum(degres)
            pcn = (1 - rho) / 2

            i_values = np.arange(t, l + 1)
            binom_coeffs = comb(l, i_values)
            p_cn_powers = np.power(pcn, i_values)
            one_minus_p_cn_powers = np.power(1 - pcn, l - i_values)

            pvn = np.sum(binom_coeffs * p_cn_powers * one_minus_p_cn_powers)

            all_pvn.append(last_pvn)
            f.append(pvn - last_pvn)
            last_pvn = pvn

        rho_x = ''
        for i in range(len(degres)):
            rho_x += f'\\frac{{{degres[i]:.0f}}}{{{np.sum(degres):.0f}}}x^{{{R[0][i] - 1:.0f}}}+'

        plt.scatter(all_pvn, f, label=f'$\\rho(x) = {rho_x[:-1]}$, $l = {l:.0f}$, $p = {p:.3f}$, $t = {t}$')

def density_evolution_p(H: np.ndarray, ps: float, repeat: int = 10) -> np.ndarray:
    l = np.sum(H[:, 0])
    R = np.unique(np.sum(H, axis=1), return_counts=True)
    degres = R[0] * R[1]

    t = int(l // 2 + l % 2)

    list_pvn = []

    for p in ps:
        last_pvn = p

        for _ in range(repeat):
            rho = np.sum((degres * np.power((1 - 2 * last_pvn), R[0] - 1))) / np.sum(degres)
            pcn = (1 - rho) / 2

            i_values = np.arange(t, l + 1)
            binom_coeffs = comb(l, i_values)
            p_cn_powers = np.power(pcn, i_values)
            one_minus_p_cn_powers = np.power(1 - pcn, l - i_values)

            pvn = np.sum(binom_coeffs * p_cn_powers * one_minus_p_cn_powers)
            last_pvn = pvn 
        list_pvn.append(last_pvn)
    
    rho_x = ''
    for i in range(len(degres)):
        rho_x += f'\\frac{{{degres[i]:.0f}}}{{{np.sum(degres):.0f}}}x^{{{R[0][i] - 1:.0f}}}+'
    
    return np.array(list_pvn)
