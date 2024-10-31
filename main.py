def read_matrix_from_file(filename):
    matrix = []
    try:
        with open(filename, 'r') as file:
            for line in file:
                elements = list(map(int, line.strip().split()))
                if elements:
                    matrix.append(elements)
    except FileNotFoundError:
        print(f"Файл '{filename}' не найден.")
    return matrix

def print_file(filename):
    try:
        with open(filename, 'r') as file:
            for line in file:
                print(line, end='')
    except FileNotFoundError:
        print(f"Файл '{filename}' не найден.")

def print_matrix(matrix):
    for row in matrix:
        print(" ".join(map(str, row)))

def process_matrix(matrix):
    def print_matrix(matrix):
        for row in matrix:
            print(' '.join(map(str, row)))
        print()
    
    print("Исходная матрица:")
    print_matrix(matrix)
    
    error_rows = []

    for row_index, row in enumerate(matrix):
        positions = [f"r{index + 1}" for index, value in enumerate(row) if value == 1]
        positions_str = ' + '.join(positions)

        xor_sum = 0
        for value in row:
            xor_sum ^= value

        print(f"S{row_index + 1} = {positions_str} = {xor_sum}")
        
        if xor_sum == 1:
            error_rows.append(row_index)
            print(f"Ошибка в S{row_index + 1}")

    corrected = False
   
    for i in range(len(error_rows)):
        for j in range(i + 1, len(error_rows)):
            row_i = error_rows[i]
            row_j = error_rows[j]
            
            common_bits = [index for index in range(len(matrix[row_i])) if matrix[row_i][index] == 1 and matrix[row_j][index] == 1]
            
            for common_bit in common_bits:
                
                bit_in_error_free_row = False
                for row_index, row in enumerate(matrix):
                    if row_index not in error_rows and row[common_bit] == 1:
                        bit_in_error_free_row = True
                        break
                
                if not bit_in_error_free_row:
                    matrix[row_i][common_bit] = 0
                    matrix[row_j][common_bit] = 0
                    corrected = True
                    print(f"Исправляем r{common_bit + 1} в строках S{row_i + 1} и S{row_j + 1}")
                    print("Исправленная матрица:")
                    print_matrix(matrix)

    if corrected:
        remaining_errors = []
        for row_index, row in enumerate(matrix):
            xor_sum = 0
            for value in row:
                xor_sum ^= value
            
            if xor_sum == 1:
                remaining_errors.append(row_index)
                print(f"Оставшаяся ошибка в S{row_index + 1}")

    return matrix  # Вернуть исправленную матрицу

def swap_rows(matrix, row1, row2):
    matrix[row1], matrix[row2] = matrix[row2], matrix[row1]

def xor_rows(matrix, target_row, source_row):
    matrix[target_row] = [a ^ b for a, b in zip(matrix[target_row], matrix[source_row])]

def gaussian_elimination(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    operations = []

    for col in range(min(rows, cols)):
        pivot_row = None
        for row in range(col, rows):
            if matrix[row][col] == 1:
                pivot_row = row
                break
        
        if pivot_row is None:
            continue

        if pivot_row != col:
            swap_rows(matrix, col, pivot_row)
            operations.append(('swap_rows', col, pivot_row))
        
        for row in range(col + 1, rows):
            if matrix[row][col] == 1:
                xor_rows(matrix, row, col)
                operations.append(('xor', row, col))

    for i in range(min(rows, cols)):
        if matrix[i][i] != 1:
            for j in range(i + 1, cols):
                if matrix[i][j] == 1:
                    for k in range(rows):
                        matrix[k][i], matrix[k][j] = matrix[k][j], matrix[k][i]
                    operations.append(('swap_cols', i, j))
                    break

    for col in range(min(rows, cols) - 1, -1, -1):
        for row in range(col):
            if matrix[row][col] == 1:
                xor_rows(matrix, row, col)
                operations.append(('xor', row, col))

    return matrix, operations

def make_identity(matrix):
    matrix, operations = gaussian_elimination(matrix)

    print("\nМатрица после приведения к диагональному виду:")
    print_matrix(matrix)

    return matrix, operations

def zero_main_diagonal(matrix):
    for i in range(min(len(matrix), len(matrix[0]))):
        matrix[i][i] = 0

def delete_zero_columns(matrix):
    non_zero_columns = [i for i in range(len(matrix[0])) if any(row[i] != 0 for row in matrix)]
    short_matrix = [[row[i] for i in non_zero_columns] for row in matrix]
    return short_matrix

def transpose(matrix):
    return list(map(list, zip(*matrix)))

def add_identity_matrix(transposed_matrix):
    num_rows = len(transposed_matrix)
    identity_matrix = [[int(i == j) for i in range(num_rows)] for j in range(num_rows)]
    augmented_matrix = [transposed_matrix[i] + identity_matrix[i] for i in range(num_rows)]
    return augmented_matrix

def make_correct_matrix(augmented_matrix, operations):
    for op in reversed(operations):
        if op[0] == 'swap_rows':
            swap_rows(augmented_matrix, op[1], op[2])
        elif op[0] == 'xor':
            xor_rows(augmented_matrix, op[1], op[2])
        elif op[0] == 'swap_cols':
            for row in augmented_matrix:
                row[op[1]], row[op[2]] = row[op[2]], row[op[1]]
    return augmented_matrix

def matrix_multiplication(matrix1, matrix2):
    result = []
    for i in range(len(matrix1)):
        row = []
        for j in range(len(matrix2[0])):
            total = 0
            for k in range(len(matrix1[0])):
                total += matrix1[i][k] * matrix2[k][j]
            row.append(total)
        result.append(row)
    return result

filename = "matrix.txt"
matrix = read_matrix_from_file(filename)
print()

print_file(filename)
print()

if matrix:
    original_matrix = process_matrix(matrix)  # Обновлено

    result_matrix, operations = make_identity(matrix.copy())
    print("\nВыполненные операции:")
    for operation in operations:
        print(operation)

    zero_main_diagonal(result_matrix)
    print("\nОбнуляем:")
    for row in result_matrix:
        print(row)
    
    result_matrix = delete_zero_columns(result_matrix)
    print("\nУдаляем:")
    for row in result_matrix:
        print(row)
    
    transposed_matrix = transpose(result_matrix)
    print("\nТранспонируем:")
    for row in transposed_matrix:
        print(row)
    
    augmented_matrix = add_identity_matrix(transposed_matrix)
    print("\nДобавляем единичную:")
    for row in augmented_matrix:
        print(row)

    corrected_matrix = make_correct_matrix(augmented_matrix, operations)
    print("\nКорректируем матрицу:")
    for row in corrected_matrix:
        print(row)

    transposed_matrix1 = transpose(corrected_matrix)
    print("\nТранспонируем:")
    for row in transposed_matrix1:
        print(row)

    print("\nМатрица 1 (исходная матрица):")
    print_matrix(original_matrix)

    print("\nМатрица 2 (транспонированная и корректированная):")
    print_matrix(transposed_matrix1)

    should_be_zero_matrix = matrix_multiplication(original_matrix, transposed_matrix1)
    print("\nУмножаем:")
    for row in should_be_zero_matrix:
        print