/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"
#include <string.h>

int n;

// Perform matrix multiplication C = A x B
static void matrix_mul(double *C, double *A, double *B) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i*n + j] = 0;
            for (int k = 0; k < n; k++) {
                C[i*n + j] += A[i*n + k] * B[k*n + j];
            }
        }
    }
}

// Perform matrix addition C = A + B
static int matrix_add(double *C, double *A, double *B) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i*n + j] = A[i*n + j] + B[i*n + j];
        }
    }

    return 1;
}

// Perform matrix tranposition A = A^t
static int matrix_sq_transpose_inplace(double *A) {
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            // Swap a[i][j] with a[j][i]
            double aux = A[i*n + j]; 
            A[i*n + j] = A[j*n + i];
            A[j*n + i] = aux;
        }
    }

    return 1;
}

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	printf("BASIC SOLVER\n");

    // Initialisation
    n = N;
    double *C = calloc(n*n, sizeof(*C));
    double *D = calloc(n*n, sizeof(*D));
    double *E = calloc(n*n, sizeof(*D));

    // D = B * A * A^t
    matrix_mul(C, B, A);
    matrix_sq_transpose_inplace(A);
    matrix_mul(D, C, A);

    // C = B^t*B
    memcpy(A, B, n*n*sizeof(*A));
    matrix_sq_transpose_inplace(B);
    matrix_mul(C, B, A);

    // E = C + D
    matrix_add(E, C, D);

	return E;
}
