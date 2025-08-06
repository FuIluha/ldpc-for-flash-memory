#include "differential_evolution_poly_find.h"
#include "density_evolution_poly.h"
#include "progress_bar.h"
#include <math.h>
#include <stdlib.h>

typedef struct Bounds {
    long double p_0;
    long double p_minus;
} Bounds;

typedef struct PolyItem {
    int d;
    long double poly;
} PolyItem;

typedef struct Poly {
    PolyItem* items;
    int size;
} Poly;

static long double random_double() {
    return rand() / (RAND_MAX + 1.0);
}

static void shuffle(int* arr, int size) {
    for (int i = size - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }
}

static long double h2(long double p) {
    // Расчитывает энтропию Шеннона
    long double answer = 0;
    if (p != 0)
        answer += p * log2l(p);
    if (p != 1)
        answer += (1 - p) * log2l(1 - p);
    return answer;
}

static Bounds calc_bounds(double r) {
    // Расчитывает граничные вероятности декодирования
    Bounds bounds = {1 - r, 0};
    long double a = 0, b = 0.5;

    for (int i = 0; i < 300; i++) {
        long double cur = (a + b) / 2;
        long double cur_p_minus_bound = 1 - h2(cur);
        if (fabsl(b - a) < 1e-20) {
            bounds.p_minus = cur;
            return bounds;
        }
        if (cur_p_minus_bound < r)
            a = cur;
        else
            b = cur;
    }
    bounds.p_minus = (a + b) / 2;
    return bounds;
}

static void calc_rho(
    double r, 
    Poly lambda,
    Poly* rho
) {
    // Подбирает полином rho
    long double sum_lmbd = 0;
    for (int i = 0; i < lambda.size; i++) {
        sum_lmbd += lambda.items[i].poly / lambda.items[i].d;
    }
    long double k = 1 / ((1 - r) * sum_lmbd);
    long double x = floor(k), y = ceil(k);
    if (x == y) {
        rho->items = (PolyItem*)malloc(sizeof(PolyItem));
        rho->items[0].poly = k;
        rho->items[0].d = 1;
        rho->size = 1;
        return;
    }
    rho->items = (PolyItem*)malloc(2 * sizeof(PolyItem));
    rho->size = 2;
    rho->items[0].d = x;
    rho->items[1].d = y;
    rho->items[0].poly = x * (y * (1 - r) * sum_lmbd - 1) / (y - x);
    rho->items[1].poly = 1 - rho->items[0].poly;
    return;
}

static long double calc_f(
    int num_iters,
    int dv_max,
    double r,
    long double* p_0_vals,
    long double* p_minus_vals,
    int p_size,
    long double threshold,
    Poly lambda
) {
    // Расчитывает площадь под графиком со знаком минус, функция, которую мы будем минимизировать
    Poly rho;
    calc_rho(r, lambda, &rho);
    int* lambda_d_values = (int*)malloc((lambda.size) * sizeof(int));
    long double* lambda_poly_values = (long double*)malloc((lambda.size) * sizeof(long double));
    for (int i = 0; i < lambda.size; i++) {
        lambda_d_values[i] = lambda.items[i].d;
        lambda_poly_values[i] = lambda.items[i].poly;
    }
    int* rho_d_values = (int*)malloc((rho.size) * sizeof(int));
    long double* rho_poly_values = (long double*)malloc((rho.size) * sizeof(long double));
    for (int i = 0; i < rho.size; i++) {
        rho_d_values[i] = rho.items[i].d;
        rho_poly_values[i] = rho.items[i].poly;
    }
    free(rho.items);

    long double* border = (long double*)malloc((p_size) * sizeof(long double));
    density_evolution_grid_irregular(
        p_0_vals, 
        p_size, 
        p_minus_vals, 
        p_size, 
        lambda_d_values, 
        rho_d_values, 
        lambda_poly_values, 
        rho_poly_values, 
        lambda.size, 
        rho.size, 
        num_iters, 
        threshold,
        0,
        border
    );
    free(rho_d_values);
    free(rho_poly_values);
    free(lambda_d_values);
    free(lambda_poly_values);

    long double square = 0;
    for (int i = 1; i < p_size && border[i] >= 0; i++) {
        square += (p_0_vals[i] - p_0_vals[i - 1]) * (border[i] + border[i - 1]) / 2;
    }
    free(border);

    return -square;
}


__attribute__((visibility("default")))
void differential_evolution_poly_find(
    int dv_max, 
    double r, 
    long double threshold, 
    int np, 
    double cr, 
    double f, 
    int max_iters, 
    long double* lambda,
    long double* rho,
    int* dc,
    int* dc_size
) {
    srand(666);
    int p_size = 600;
    int decoding_iters = 100;
    Bounds bounds = calc_bounds(r);

    long double* p_0_vals = (long double*)malloc((p_size) * sizeof(long double));
    long double* p_minus_vals = (long double*)malloc((p_size) * sizeof(long double));
    long double delta_p_minus = bounds.p_minus / (p_size - 1), delta_p_0 = bounds.p_0 / (p_size - 1);
    for (int i = 0; i < p_size; i++) {
        p_0_vals[i] = i * delta_p_0;
        p_minus_vals[i] = i * delta_p_minus;
    }

    Poly* lambda_values = (Poly*)malloc((np) * sizeof(Poly));
    for (int i = 0; i < np; i++) {
        lambda_values[i].items = (PolyItem*)malloc((dv_max) * sizeof(PolyItem));
        lambda_values[i].size = dv_max;

        long double sum = 0.0;
        for (int j = 0; j < dv_max; j++) {
            lambda_values[i].items[j].d = j + 1;
            lambda_values[i].items[j].poly = random_double();
            sum += lambda_values[i].items[j].poly;
        }

        for (int j = 0; j < dv_max; j++) {
            lambda_values[i].items[j].poly /= sum;
        }
    }

    for (int i = 0; i < max_iters; i++) {
        print_progress_bar(i, max_iters - 1);
        for (int x_idx = 0; x_idx < np; x_idx++) {
            int indices[np - 1], pos = 0;
            for (int j = 0; j < np; j++) {
                if (j != x_idx) indices[pos++] = j;
            }

            shuffle(indices, np - 1);

            int a_idx = indices[0];
            int b_idx = indices[1];
            int c_idx = indices[2];

            // Случайный индекс R (обязательно включить этот индекс в дифференциацию)
            int R = rand() % dv_max;

            Poly trial_vector;
            trial_vector.items = (PolyItem*)malloc(dv_max * sizeof(PolyItem));
            trial_vector.size = dv_max;

            for (int j = 0; j < dv_max; j++) {
                double r_j = (double)rand() / RAND_MAX;

                if (r_j < cr || j == R) {
                    // y_j = a_j + f * (b_j - c_j)
                    trial_vector.items[j].d = j + 1;
                    trial_vector.items[j].poly =
                        lambda_values[a_idx].items[j].poly +
                        f * (lambda_values[b_idx].items[j].poly - lambda_values[c_idx].items[j].poly);
                } else {
                    // y_j = x_j
                    trial_vector.items[j].d = j + 1;
                    trial_vector.items[j].poly = lambda_values[x_idx].items[j].poly;
                }
            }

            // (нормировка trial_vector, если нужно, например, сумма = 1)
            long double sum = 0.0;
            for (int j = 0; j < dv_max; j++) {
                if (trial_vector.items[j].poly < 0.0)
                    trial_vector.items[j].poly = 0.0; // clamp отрицательные

                sum += trial_vector.items[j].poly;
            }
            if (sum > 0.0) {
                for (int j = 0; j < dv_max; j++) {
                    trial_vector.items[j].poly /= sum;
                }
            }
            
            long double fy = calc_f(decoding_iters, dv_max, r, p_0_vals, p_minus_vals, p_size, threshold, trial_vector);
            long double fx = calc_f(decoding_iters, dv_max, r, p_0_vals, p_minus_vals, p_size, threshold, lambda_values[x_idx]);

            // Если новый полином не хуже предыдущего — принимаем его
            if (fy <= fx) {
                for (int j = 0; j < dv_max; j++) {
                    lambda_values[x_idx].items[j].poly = trial_vector.items[j].poly;
                }
            }
            {
                for (int j = 0; j < dv_max; j++) {
                    lambda_values[x_idx].items[j].poly = trial_vector.items[j].poly;
                }
            }

            free(trial_vector.items);
        }
    }

    // Найдём индекс лучшего агента
    int best_idx = 0;
    long double best_val = calc_f(decoding_iters, dv_max, r, p_0_vals, p_minus_vals, p_size, threshold, lambda_values[0]);

    for (int i = 1; i < np; i++) {
        long double val = calc_f(decoding_iters, dv_max, r, p_0_vals, p_minus_vals, p_size, threshold, lambda_values[i]);
        if (val < best_val) {
            best_val = val;
            best_idx = i;
        }
    }

    // Сохраняем лучшие коэффициенты в results
    for (int j = 0; j < dv_max; j++) {
        lambda[j] = lambda_values[best_idx].items[j].poly;
    }

    Poly poly_rho;
    calc_rho(r, lambda_values[best_idx], &poly_rho);
    for (int i = 0; i < poly_rho.size; i++) {
        rho[i] = poly_rho.items[i].poly;
        dc[i] = poly_rho.items[i].d;
    }
    *dc_size = poly_rho.size;

    for (int i = 0; i < np; i++) {
        free(lambda_values[i].items);
    }
    free(lambda_values);
    free(p_0_vals);
    free(p_minus_vals);
}
