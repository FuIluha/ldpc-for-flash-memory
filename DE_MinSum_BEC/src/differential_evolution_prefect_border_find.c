#include "models.h"
#include "progress_bar.h"
#include "density_evolution_poly.h"
#include "differential_evolution_poly_find.h"
#include "differential_evolution_prefect_border_find.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdatomic.h>
#include <unistd.h>
#include <pthread.h>
#include <omp.h>

static long double h2(long double p) {
    // Расчитывает энтропию Шеннона
    long double answer = 0;
    if (p != 0)
        answer += p * log2l(p);
    if (p != 1)
        answer += (1 - p) * log2l(1 - p);
    return -answer;
}

static Bounds calc_bounds(double r) {
    // Расчитывает граничные вероятности декодирования
    Bounds bounds = {1 - r, 0};
    long double a = 0, b = 0.5;

    for (int i = 0; i < 300; i++) {
        long double cur = (a + b) / 2;
        long double cur_p_minus_bound = 1 - h2(cur);
        if (b - a < 1e-20) {
            bounds.p_minus = cur;
            return bounds;
        }
        if (cur_p_minus_bound > r)
            a = cur;
        else
            b = cur;
    }
    bounds.p_minus = (a + b) / 2;
    return bounds;
}

// Глобальный атомарный счётчик
static atomic_int progress_counter;
static int total_tasks;

// Флаг завершения
static atomic_int progress_done;

// Поток, который рисует прогресс
void* progress_thread_func(void* arg) {
    while (!atomic_load(&progress_done)) {
        int done = atomic_load(&progress_counter);
        print_progress_bar(done, total_tasks);
        usleep(100000); // обновление каждые 0.1 сек
    }
    // Финальный вывод
    print_progress_bar(total_tasks, total_tasks);
    return NULL;
}

__attribute__((visibility("default")))
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
) {
    Bounds bounds = calc_bounds(r);

    long double* p_0_vals = (long double*)malloc((p_size) * sizeof(long double));
    long double delta_p_0 = bounds.p_0 / (p_size - 1);
    for (int i = 0; i < p_size; i++) {
        p_0_vals[i] = i * delta_p_0;
    }

    atomic_store(&progress_counter, 0);
    atomic_store(&progress_done, 0);
    total_tasks = p_size;

    // Запускаем поток прогресса
    pthread_t progress_thread;
    pthread_create(&progress_thread, NULL, progress_thread_func, NULL);

    // Запускаем вычисления параллельно
    #pragma omp parallel for num_threads(4) schedule(dynamic)
    for (int i = 0; i < p_size; i++) {
        long double* lambda = (long double*)malloc((dv_max) * sizeof(long double));
        long double* rho = (long double*)malloc((2) * sizeof(long double));
        int* dc = (int*)malloc((2) * sizeof(int));
        int dc_size;

        differential_evolution_poly_find_for_p0_prob(
            dv_max, r, threshold, np,
            cr, f, max_iters, lambda, rho,
            dc, &dc_size, p_0_vals[i], decoding_iters,
            p_minus_split_iters
        );

        int* dv = (int*)malloc((dv_max) * sizeof(int));
        for (int j = 0; j < dv_max; j++) {
            dv[j] = j + 1;
        }

        density_evolution_grid_irregular(
            (p_0_vals + i), 1, bounds.p_minus, p_minus_split_iters,
            dv, dc, lambda, rho, dv_max,
            dc_size, decoding_iters, threshold, 0,
            (p_minus + i)
        );

        free(lambda);
        free(rho);
        free(dc);
        free(dv);

        // отметим задачу
        atomic_fetch_add(&progress_counter, 1);
    }

    // Всё готово → останавливаем поток прогресса
    atomic_store(&progress_done, 1);
    pthread_join(progress_thread, NULL);

    free(p_0_vals);
}