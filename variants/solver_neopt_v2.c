/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"
#include <string.h>

#define max(a, b) (a) > (b) ? a : b;

// Compute B = A * A^t
// Knowing A is upper triangular
void mul_tr_uptri(double * B, double * A, size_t n) {
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n; j++) {
            size_t start = max(i, j);
            B[i*n + j] = 0;
            register size_t in = i*n;
            register size_t jn = j*n;
            for (size_t k = start; k < n; k++) {
                B[in + j] += A[in + k] * A[jn + k];
            }

            // The output is symetrical
            B[jn + i] = B[in + j];
        }
    }
}

void mul_tr(double * C, double * A, double * B, size_t n) {

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            C[i*n + j] = 0;
            register size_t in = i*n;
            register size_t jn = j*n;
            for (size_t k = 0; k < n; k++) {
                C[in + j] += A[in + k] * B[jn + k];
            }
        }
    }
}

// Set A = A^t
void transpose(double *A, int n) {
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            double aux = A[i*n + j];
            A[i*n + j] = A[j*n + i];
            A[j*n + i] = aux;
        }
    }
}

void mult_mat_trans(double * B, double * A, size_t n) {
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n; j++) {
            B[i*n + j] = 0;
            register size_t in = i*n;
            register size_t jn = j*n;
            for (size_t k = 0; k < n; k++) {
                B[in + j] += A[in + k] * A[jn + k];
            }

            // The output is symetrical
            B[jn + i] = B[in + j];
        }
    }
}

// Compute A = A + B
void add_mat(double *A, double *B, const size_t n) {
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
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
    mul_tr_uptri(C, A, n);

    // Compute D = B * C
    // C is symetrical
    mul_tr(D, B, C, n);

    // Tranpose: B = B^t
    transpose(B, n);

    // Compute C = B * B^t
    mult_mat_trans(C, B, n);

    // Compute C = C + D
    add_mat(C, D, n);

    free(D);
    
    return C;
}
