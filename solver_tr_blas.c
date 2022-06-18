/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include <string.h>
#include "utils.h"
#include "cblas.h"

#define BLOCK_SIZE 40

/* 
 * Add your BLAS implementation here
 */

double* my_solver(int N, double *A, double *B) {
    printf("BLAS SOLVER\n");

    register double *C = calloc(N*N, sizeof(*C));

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0f, A, N, B, N, 1.0f, C, N);

    return C;
}
