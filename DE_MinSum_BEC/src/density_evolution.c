#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "density_evolution.h"
#include "progress_bar.h"

static unsigned long long factorial(int n) {
    unsigned long long res = 1;
    for (int i = 2; i <= n; ++i)
        res *= i;
    return res;
}

static unsigned long long binomial(int n, int k) {
    if (k < 0 || k > n) return 0;
    return factorial(n) / (factorial(k) * factorial(n - k));
}

__attribute__((visibility("default")))
void density_evolution_bec(
    long double epsilon, long double delta,
    int l, int k, int num_iters,
    long double* p0_out, long double* p_minus_out, long double* p_plus_out
) {
    long double p0 = epsilon;
    long double p_minus = delta;
    long double p_plus = 1.0 - p0 - p_minus;

    p0_out[0] = p0;
    p_minus_out[0] = p_minus;
    p_plus_out[0] = p_plus;

    for (int it = 0; it < num_iters; it++) {
        long double p0_cn = 1.0 - pow(1.0 - p0, k - 1);
        long double tmp = pow(1.0 - p0 - 2.0 * p_minus, k - 1);
        long double pow_1mp0 = pow(1.0 - p0, k - 1);
        long double p_minus_cn = (pow_1mp0 - tmp) / 2.0;
        long double p_plus_cn = (pow_1mp0 + tmp) / 2.0;

        long double new_p0 = 0.0;
        long double new_p_minus = 0.0;

        int l_adj = (it == num_iters - 1) ? l + 1 : l;

        for (int i = -1; i <= 1; i++) {
            long double cur_val_by_i = (l_adj - 1) * log2l(p_minus_cn);
            long double cur_p0_by_k1= 0.0;
            long double cur_p_minus_by_k1 = 0.0;
            long double theta = 0.0;
            switch (i) {
            case -1:
                theta = delta;
                break;
            case 0:
                theta = epsilon;
                break;
            default:
                theta = 1 - delta - epsilon;
                break;
            }

            unsigned short is_p0_changed = 0;
            unsigned short is_p_minus_changed = 0;
            for (int k1 = 0; k1 <= l_adj - 1; k1++) {
                long double cur_val_by_k1 = (l_adj - 1 - k1) * log2l(p_plus_cn / p_minus_cn) + log2l(binomial(l_adj - 1, k1));
                long double cur_p0_by_k2 = 0.0;
                long double cur_p_minus_by_k2 = 0.0;
                if (k1 == 0 && p_minus_cn <= 0) {
                    cur_val_by_i = 0;
                    cur_val_by_k1 = (l_adj - 1 - k1) * log2l(p_plus_cn) + log2l(binomial(l_adj - 1, k1));
                } else if (p_minus_cn <= 0) {
                    break;
                }

                if (theta > 0) {
                    unsigned short is_p0_changed_local = 0;
                    unsigned short is_p_minus_changed_local = 0;
                    for (int k2 = 0; k2 <= l_adj - 1 - k1; k2++) {
                        long double cur_val_by_k2 = log2l(binomial(l_adj - 1 - k1, k2)) + k2 * log2l(p0_cn / p_plus_cn);
                        if (l_adj - 1 - k1 - k2 == 0 && p_plus_cn <= 0) {
                            cur_val_by_k1 = (1 + k1 - l_adj) * log2l(p_minus_cn) + log2l(binomial(l_adj - 1, k1));
                            cur_val_by_k2 = log2l(binomial(l_adj - 1 - k1, k2)) + k2 * log2l(p0_cn);
                        } else if (p_plus_cn <= 0) {
                            break;
                        }
                        if (k2 == 0 && p0_cn <= 0) {
                            cur_val_by_k2 = log2l(binomial(l_adj - 1 - k1, k2));
                        } else if (p0_cn <= 0) {
                            break;
                        }
                        if (l_adj - 1 - 2 * k1 - k2 + i == 0) {
                            is_p0_changed_local = 1;
                            cur_p0_by_k2 += exp2l(cur_val_by_k2 + log2l(theta));
                        } else if (l_adj - 1 - 2 * k1 - k2 + i < 0) {
                            is_p_minus_changed_local = 1;
                            cur_p_minus_by_k2 += exp2l(cur_val_by_k2 + log2l(theta));
                        }
                    }
                    if (is_p0_changed_local) {
                        is_p0_changed = 1;
                        cur_p0_by_k1 += exp2l(cur_val_by_k1 + log2l(cur_p0_by_k2));
                    }
                    if (is_p_minus_changed_local) {
                        is_p_minus_changed = 1;
                        cur_p_minus_by_k1 += exp2l(cur_val_by_k1 + log2l(cur_p_minus_by_k2));
                    }
                }
            }
            if (is_p0_changed)
                new_p0 += exp2l(cur_val_by_i + log2l(cur_p0_by_k1));
            if (is_p_minus_changed)
                new_p_minus += exp2l(cur_val_by_i + log2l(cur_p_minus_by_k1));
        }

        p0_out[it + 1] = new_p0;
        p_minus_out[it + 1] = new_p_minus;
        p_plus_out[it + 1] = 1 - new_p0 - new_p_minus;

        p0 = new_p0;
        p_minus = new_p_minus;
    }
}

__attribute__((visibility("default")))
void density_evolution_grid(
    long double* p0_vals, int n_p0,
    long double* p_minus_vals, int n_p_minus,
    int l, int k, int num_iters, long double threshold,
    int* results
) {
    int total = n_p0 - 1;
    long double last_p_minus_idx = 0;

    for (int i = 0; i < n_p0; ++i) {
        results[i] = -1;
        print_progress_bar(i, total);
        
        int step = -1;
        for (int j = last_p_minus_idx; j < n_p_minus && j >= 0; j += step) {
            long double p0 = p0_vals[i];
            long double p_minus = p_minus_vals[j];

            long double* p0_out = (long double*)malloc((num_iters + 1) * sizeof(long double));
            long double* p_minus_out = (long double*)malloc((num_iters + 1) * sizeof(long double));
            long double* p_plus_out = (long double*)malloc((num_iters + 1) * sizeof(long double));

            if (!p0_out || !p_minus_out || !p_plus_out) {
                if (p0_out) free(p0_out);
                if (p_minus_out) free(p_minus_out);
                if (p_plus_out) free(p_plus_out);
                continue;
            }

            density_evolution_bec(
                p0, p_minus, l, k, num_iters,
                p0_out, p_minus_out, p_plus_out
            );

            long double cur_ber = p0_out[num_iters] + p_minus_out[num_iters];

            free(p0_out);
            free(p_minus_out);
            free(p_plus_out);

            if (threshold > cur_ber) {
                results[i] = j;
                step = 1;
            } else if (step > 0) {
                break;
            }
        }
        last_p_minus_idx = results[i];
    }
}