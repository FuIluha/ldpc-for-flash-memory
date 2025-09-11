#ifndef DENSITY_EVOLUTION_POLY_H
#define DENSITY_EVOLUTION_POLY_H

void density_evolution_bec_irregular(
    long double epsilon, 
    long double delta,
    int* dv,
    int* dc,
    long double* lambda,
    long double* rho,
    int dv_size,
    int dc_size,
    int num_iters,
    long double* p0_out, 
    long double* p_minus_out, 
    long double* p_plus_out
);

void density_evolution_grid_irregular(
    long double* p0_vals, 
    int n_p0,
    long double p_minus_limit,
    int p_minus_max_iters,
    int* dv,
    int* dc,
    long double* lambda,
    long double* rho,
    int dv_size,
    int dc_size,
    int num_iters, 
    long double threshold,
    int is_progress_visible,
    long double* results
);

#endif