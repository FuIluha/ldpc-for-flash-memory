#ifndef DENSITY_EVOLUTION_H
#define DENSITY_EVOLUTION_H

void density_evolution_bec(
    long double epsilon, long double delta,
    int l, int k, int num_iters,
    long double* p0_out, long double* p_minus_out, long double* p_plus_out
);

void density_evolution_grid(
    long double* p0_vals, int n_p0,
    long double* p_minus_vals, int n_p_minus,
    int l, int k, int num_iters, long double threshold,
    int* results
);

#endif