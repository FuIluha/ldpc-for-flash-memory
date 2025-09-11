#ifndef DIFFERENTIAL_EVOLUTION_PERFECT_BORDER_FIND_H
#define DIFFERENTIAL_EVOLUTION_PERFECT_BORDER_FIND_H

void differential_evolution_perfect_border_find(
    int dv_max,
    double r,
    long double threshold,
    int np,
    double cr,
    double f,
    int max_iters,
    int p_size,
    int decoding_iters,
    int p_minus_split_iters,
    long double* p_minus
);

#endif