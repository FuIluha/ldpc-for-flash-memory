import numpy as np

def dot_for_codes(a: np.array, b: np.array) -> np.array:
    return np.array(a @ b % 2)
# мои примеры
H = np.loadtxt('verification_matrix.txt', delimiter=',', dtype=int)

r_messages = np.loadtxt('error_messages.txt', delimiter=',', dtype=int)

# Вычисление синдромов S другим способом, потому что иначе размерность не совпадет
S = []
for r in r_messages:
    S.append(np.mod(np.dot(H, r.T), 2).T )

corrected_r = []
for current_r, current_S in zip(r, S):
    if not np.any(current_S):  
        corrected_r.append(current_r)
        continue

    for i, h_column in enumerate(H.T):
        if np.array_equal(h_column, current_S):  
            print(f"Ошибка найдена в бите {i}.")
            current_r[i] = (current_r[i] + 1) % 2  
            break
    corrected_r.append(current_r)
corrected_r = np.array(corrected_r)
print("Исправленные кодовые слова:")
print(corrected_r)

