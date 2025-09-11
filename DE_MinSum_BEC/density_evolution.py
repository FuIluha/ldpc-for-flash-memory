import numpy as np
from scipy.special import factorial

def density_evolution_bec(epsilon, delta, l, k, num_iters):
    """
    Моделирует эволюцию плотностей p0 и p_minus для LDPC декодера.
    
    Аргументы:
        l (int): степень variable node (d_v)
        k (int): степень check node (d_c)
        epsilon (float): вероятность стирания в канале
        delta (float): вероятность ошибки (бит-флипа) в канале
        num_iterations (int): количество итераций

    Возвращает:
        p0_list, p_minus_list — массивы длинной num_iterations+1
    """
    epsilon = np.float64(epsilon)
    delta = np.float64(delta)

    # Начальные вероятности из канала
    p0 = epsilon
    p_minus = delta
    p_plus = 1 - p0 - p_minus

    p0_list = [p0]
    p_minus_list = [p_minus]
    p_plus_list = [1 - p0 - p_minus]

    for it in range(num_iters):
        # 1. Вычисляем вероятности на check nodes
        p0_cn = 1 - np.power(1 - p0, k - 1)
        p_minus_cn = ((1 - p0)**(k - 1) - (1 - p0 - 2 * p_minus)**(k - 1)) / 2
        p_plus_cn = ((1 - p0)**(k - 1) + (1 - p0 - 2 * p_minus)**(k - 1)) / 2

        # 2. Вычисляем вероятности на variable nodes
        new_p0 = 0
        new_p_minus = 0
        new_p_plus = 0

        if it == num_iters - 1:
            l += 1

        for k1 in range(l):
            for k2 in range(l - k1):
                k3 = l - 1 - k1 - k2

                # Вероятность конкретной конфигурации
                if (p_minus_cn > 0 or k1 == 0) and (p0_cn > 0 or k2 == 0) and (p_plus_cn > 0 or k3 == 0):
                    prob = 0
                    if k1 > 0:
                        prob += k1 * np.log(p_minus_cn)
                    if k2 > 0:
                        prob += k2 * np.log(p0_cn)
                    if k3 > 0:
                        prob += k3 * np.log(1 - p0_cn - p_minus_cn)
                    prob = np.exp(prob)
                else:
                    prob = 0

                for j in [-1, 0, +1]:
                    sum_j = -k1 + k3 + j
                    base = delta * (j == -1) + epsilon * (j == 0) + (1 - epsilon - delta) * (j == +1)

                    if sum_j == 0:
                        new_p0 += prob * base * factorial(l - 1) / (factorial(k1) * factorial(k2) * factorial(k3))
                    elif sum_j > 0:
                        new_p_plus += prob * base * factorial(l - 1) / (factorial(k1) * factorial(k2) * factorial(k3))
                    elif sum_j < 0:
                        new_p_minus += prob * base * factorial(l - 1) / (factorial(k1) * factorial(k2) * factorial(k3))

        # Сохраняем
        p0_list.append(new_p0)
        p_minus_list.append(new_p_minus)
        p_plus_list.append(new_p_plus)

        # Обновляем
        p0 = new_p0
        p_minus = new_p_minus

    return np.array(p0_list), np.array(p_minus_list), np.array(p_plus_list)