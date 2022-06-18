/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"
#include <string.h>

#define max(a, b) (a) > (b) ? a : b;

// Compute B = A * A^t
// Knowing A is upper triangular
void mult_mat_trans_tri(register double * B, register double * A, register size_t n) {
    for (register size_t i = 0; i < n; i++) {
        for (register size_t j = i; j < n; j++) {
            size_t start = max(i, j);
            B[i*n + j] = 0;
            for (register size_t k = start; k < n; k++) {
                B[i*n + j] += A[i*n + k] * A[j*n + k];
            }

            // The output is symetrical
            B[j*n + i] = B[i*n + j];
        }
    }
}

void mul_sym(register double * C, register double * A, register double * B, register size_t n) {

    for (register size_t i = 0; i < n; i++) {
        for (register size_t j = 0; j < n; j++) {
            C[i*n + j] = 0;
            for (register size_t k = 0; k < n; k++) {
                C[i*n + j] += A[i*n + k] * B[j*n + k];
            }
        }
    }
}

// Set A = A^t
void transpose(register double *A, register int n) {
    for (register size_t i = 0; i < n; i++) {
        for (register size_t j = i + 1; j < n; j++) {
            register double aux = A[i*n + j];
            A[i*n + j] = A[j*n + i];
            A[j*n + i] = aux;
        }
    }
}

void mult_mat_trans(register double * restrict B, register double * restrict A, register size_t n) {
    for (register size_t i = 0; i < n; i++) {
        for (register size_t j = i; j < n; j++) {
            B[i*n + j] = 0;
            for (register size_t k = 0; k < n; k++) {
                B[i*n + j] += A[i*n + k] * A[j*n + k];
            }

            // The output is symetrical
            B[j*n + i] = B[i*n + j];
        }
    }
}

// Compute A = A + B
void add_mat(register double *A, register double *B, const register size_t n) {
    for (register size_t i = 0; i < n; i++) {
        for (register size_t j = 0; j < n; j++) {
            A[i*n + j] += B[i*n + j];
        }
    }
}

/*
 * Unoptimized implementation of B*A*A^t + B^t*B
 */
double* my_solver(int N, double *A, double* B) {
    printf("NEOPT SOLVER\n");
    double *C = malloc(N*N * sizeof(*C));
    double *D = malloc(N*N * sizeof(*D));
    size_t n = N;

    // Compute C = A * A^t
    mult_mat_trans_tri(C, A, n);

    // Compute D = B * C
    // C is symetrical
    mul_sym(D, B, C, n);

    // Tranpose: B = B^t
    transpose(B, n);

    // Compute C = B * B^t
    mult_mat_trans(C, B, n);

    // Compute C = C + D
    add_mat(C, D, n);

    free(D);
    
    return C;
}
