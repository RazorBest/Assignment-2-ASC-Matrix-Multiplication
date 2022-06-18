/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"
#include <string.h>

#define max(a, b) (a) > (b) ? a : b;

#define BLOCK_SIZE 40

// Set A = A^t
void transpose(double *A, int n) {
    register size_t i, j, k;
    // Transpose the blocks that are on the main diagonal
    register size_t i_start = 0;
    register size_t n_blocks = n / BLOCK_SIZE;
    for (k = 0; k < n_blocks; k++) {
        for (i = i_start; i < i_start + BLOCK_SIZE; i++) {
            double *A_lin_1 = &A[i*n + (i + 1)];
            double *A_col_1 = &A[(i + 1)*n + i];
            for (j = i + 1; j < i_start + BLOCK_SIZE; j++) {
                double aux = *A_lin_1;
                *A_lin_1 = *A_col_1;
                *A_col_1 = aux;

                A_lin_1++;
                A_col_1 += n;
            }
        }
    }

    for (size_t lb = 0; lb < n_blocks; lb++) {
        register size_t j_start = (lb + 1) *  BLOCK_SIZE;
        for (register size_t rb = lb + 1; rb < n_blocks; rb++) {
            register double *A_lin_2 = &A[(lb*BLOCK_SIZE)*n + j_start];
            for (i = lb*BLOCK_SIZE; i < lb*BLOCK_SIZE + BLOCK_SIZE; i++) {
                register double *A_lin_1 = A_lin_2;//&A[i*n + j_start];
                register double *A_col_1 = &A[j_start*n + i];
                for (j = j_start; j < j_start + BLOCK_SIZE; j++) {
                    double aux = *A_lin_1;
                    *A_lin_1 = *A_col_1;
                    *A_col_1 = aux;

                    A_lin_1++;
                    A_col_1 += n;
                }
                A_lin_2 += n;
            }
            j_start += BLOCK_SIZE;
        }
    }
}

/*
 * Add your optimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
    printf("OPT SOLVER\n");
    double *C = malloc(N*N * sizeof(*C));
    double *D = malloc(N*N * sizeof(*D));
    size_t n = N;
    register size_t i;
    register size_t j;
    register size_t k;

    // Compute C = A * A^t
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            size_t start = max(i, j);
            C[i*n + j] = 0;
            register double sum = 0;
            register double *A_lin_1 = &A[i*n + start];
            register double *A_lin_2 = &A[j*n + start];
            for (k = start; k < n; k++) {
                sum += *A_lin_1 * (*A_lin_2);
                A_lin_1++;
                A_lin_2++;
            }
            C[i*n + j] = sum;

            // The output is symetrical
            C[j*n + i] = C[i*n + j];
        }
    }

    /*
    // Compute D = C * B^t
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            D[i*n + j] = 0;
            for (k = 0; k < n; k++) {
                D[i*n + j] += C[i*n + k] * B[j*n + k];
            }
        }
    }
    */

    // Compute D = B * C
    // C is symetrical
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            D[i*n + j] = 0;
            register double sum1 = 0;
            register double sum2 = 0;
            for (k = 0; k < n; k += 2) {
                sum1 += B[i*n + k] * C[j*n + k];
                sum2 += B[i*n + k + 1] * C[j*n + k + 1];
            }
            D[i*n + j] = sum1 + sum2;
        }
    }

    // Tranpose: B = B^t
    transpose(B, n);

    // Compute C = B * B^t
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            register double sum1 = 0;
            register double sum2 = 0;
            for (k = 0; k < n; k += 2) {
                sum1 += B[i*n + k] * B[j*n + k];
                sum2 += B[i*n + k] * B[j*n + k];
            }
            C[i*n + j] = sum1 + sum2;
            // The output is symmetrical
            C[j*n + i] = C[i*n + j];
        }
    }


    // Compute C = C + D
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j += 2) {
            C[i*n + j] += D[i*n + j];
            C[i*n + j + 1] += D[i*n + j + 1];
        }
    }

    free(D);
    
    return C;
}
