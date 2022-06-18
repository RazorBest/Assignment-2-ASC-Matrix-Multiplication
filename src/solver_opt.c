/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include <string.h>
#include <emmintrin.h>
#include "utils.h"
#include "./mat_mul.h"

#define BLOCK_SIZE 40

// Set A = A^t
void transpose(double *A, size_t n) {
    register size_t i, j, k;
    register size_t i_start = 0;
    register size_t n_blocks = n / BLOCK_SIZE;

    // Transpose the blocks that are on the main diagonal
    // Iterates though the blocks
    for (k = 0; k < n_blocks; k++) {
        for (i = i_start; i < i_start + BLOCK_SIZE; i++) {
            register double *A_lin_1 = &A[i*n + (i + 1)];
            register double *A_col_1 = &A[(i + 1)*n + i];
            for (j = i + 1; j < i_start + BLOCK_SIZE; j++) {
                register double aux = *A_lin_1;
                *A_lin_1 = *A_col_1;
                *A_col_1 = aux;

                A_lin_1++;
                A_col_1 += n;
            }
        }
        i_start += BLOCK_SIZE;
    }

    // Transpose the rest of the blocks
    for (size_t lb = 0; lb < n_blocks; lb++) {
        register size_t j_start = (lb + 1) *  BLOCK_SIZE;
        for (register size_t rb = lb + 1; rb < n_blocks; rb++) {
            register double *A_lin_2 = &A[(lb*BLOCK_SIZE)*n + j_start];
            // Swap+transpose the blocks
            for (i = lb*BLOCK_SIZE; i < lb*BLOCK_SIZE + BLOCK_SIZE; i++) {
                register double *A_lin_1 = A_lin_2;
                register double *A_col_1 = &A[j_start*n + i];
                for (j = j_start; j < j_start + BLOCK_SIZE; j++) {
                    register double aux = *A_lin_1;
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

// Compute A = A + B
void add_mat(register double * restrict A, 
        register double * restrict B, const register size_t n) {
    register size_t i;
    register size_t j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j += 2) {
            A[i*n + j] += B[i*n + j];
            A[i*n + j + 1] += B[i*n + j + 1];
        }
    }
}

/*
 * Optimized implementation of B*A*A^t + B^t*B
 */
double* my_solver(int N, register double *A, double* B) {
    printf("OPT SOLVER\n");

    double *C = malloc(N*N * sizeof(*C));
    double *D = malloc(N*N * sizeof(*D));
    register size_t n = N;

    // Compute C = A * A^t
    // Knowing A is upper triangular
    mul_tr_uptri(C, A, n);

    // Compute D = B * C
    // C is symetrical, so B * C = B * C^t
    mul_tr(D, B, C, n);

    // Tranpose: B = B^t
    transpose(B, n);

    // Compute C = B * B^t
    mul_tr_self(C, B, n);

    // Compute C = C + D
    add_mat(C, D, n);

    free(D);
    
    return C;
}
