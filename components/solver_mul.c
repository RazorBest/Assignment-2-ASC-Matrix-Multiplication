/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"

// Solver that performs matrix multiplication A*B
double* my_solver(int N, register double *A, register double* B) {
    printf("MUL SOLVER\n");
    register double *C = malloc(N*N * sizeof(*C));
    register size_t n = N;
    register size_t i, j, k;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double sum = 0;
            for (k = 0; k < n; k++) {
                sum += A[i*n + k] * B[k*n + j];
            }
            C[i*n + j] = sum;
        }
    }

    return C;
}
