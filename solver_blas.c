/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include <string.h>
#include "utils.h"
#include <cblas.h>

#define BLOCK_SIZE 40

/* 
 * Add your BLAS implementation here
 */
double* my_solver(int N, double *A, double *B) {
    printf("BLAS SOLVER\n");

    register double *C = calloc(N*N, sizeof(*C));
    register double *D = calloc(N*N, sizeof(*D));

    // C = A
    memcpy(C, A, N*N*sizeof(A));

    // C = C*A^t
    // A is upper triangular
    cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasTrans, CblasNonUnit, N, N, 1.0, A, N, C, N);

    // D = B^t*B + D
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, N, N, 1.0, B, N, B, N, 1.0, D, N);

    // D = B*C + D
    // C is symmetric
    cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, N, N, 1.0, C, N, B, N, 1.0, D, N);

    free(C);

    return D;
}
