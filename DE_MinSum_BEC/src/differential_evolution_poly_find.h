#ifndef DIFFERENTIAL_EVOLUTION_POLY_FIND_H
#define DIFFERENTIAL_EVOLUTION_POLY_FIND_H

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
    int* dc_size,
    int type_of_metric,
    int p_size,
    int decoding_iters,
    int p_minus_split_iters
);

void differential_evolution_poly_find_for_p0_prob(
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
    int* dc_size,
    long double p_0,
    int decoding_iters,
    int p_minus_split_iters
);

#endif