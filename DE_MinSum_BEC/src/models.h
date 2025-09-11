#ifndef MODELS_H
#define MODELS_H

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

#endif