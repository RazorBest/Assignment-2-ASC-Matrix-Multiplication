/*
 * Source code that contains functions for matrix multiplication
 * It only contains multiplication by matrix transpose.
 * */
#include <emmintrin.h>
#include <string.h>
#include "./asm_mul.h"
#include "./mat_mul.h"

#define min(a, b) ((a) < (b) ? (a) : (b))

// Compute B = A * A^t
// Knowing A is upper triangular
void mul_tr_uptri(double * restrict B, double * restrict A, register size_t n) {
    size_t i;
    register size_t j;

    memset(B, 0, n*n*sizeof(B));
    for (i = 0; i < n; i++) {
        // Take only the elements above the diagonal
        for (j = i; j < n; j++) {
            // Set start such that is divisible by 8
            size_t start = j / 8 * 8;
            register __m128d sum1 = _mm_set_pd(0.0, 0.0);
            register __m128d sum2 = _mm_set_pd(0.0, 0.0);

            register double *A_lin_1 = &A[i*n + start];
            register double *A_lin_2 = &A[j*n + start];
            register double *A_fin = &A[i*n + n];

            while (A_lin_1 < A_fin) {
                ASM_MUL_ADD_8(sum1, sum2, A_lin_1, A_lin_2);

                A_lin_1 += 8;
                A_lin_2 += 8;
            }

            // Perform B[i*n + j] += sum1 + sum2;
            B += i*n + j;
            ASM_ACCUMULATE_S1_S2(sum1, sum2, B);
            B -= i*n + j;

            // The output is symetrical
            B[j*n + i] = B[i*n + j];
        }
    }
    
}


// Compute A*B^t, where A is a N1 x M block and B is a M x N2 block
// All A, B and C are subblocks of matrices that are N x N
// The output block will be a N1 x N2 matrix and will be added
//  to the block that starts at C
// In this case, M must be 64
static inline void mul_tr_block_64(register double * restrict C, register double * restrict A, 
        register double * restrict B, register size_t N1, size_t M,
        register size_t N2, register size_t N) {
    register size_t i;
    register size_t j;

    register double *A_lin = A;
    // For each line in A
    for (i = 0; i < N1; i++) {
        register double *B_lin = B;
        // For each line in B
        for (j = 0; j < N2; j++) {
            // Perform dot product for 64 pairs (8x8)
            register __m128d sum1 = _mm_set_pd(0.0, 0.0);
            register __m128d sum2 = _mm_set_pd(0.0, 0.0);

            ASM_MUL_ADD_8(sum1, sum2, A_lin, B_lin);
            A_lin += 8;
            B_lin += 8;

            ASM_MUL_ADD_8(sum1, sum2, A_lin, B_lin);
            A_lin += 8;
            B_lin += 8;

            ASM_MUL_ADD_8(sum1, sum2, A_lin, B_lin);
            A_lin += 8;
            B_lin += 8;

            ASM_MUL_ADD_8(sum1, sum2, A_lin, B_lin);
            A_lin += 8;
            B_lin += 8;

            ASM_MUL_ADD_8(sum1, sum2, A_lin, B_lin);
            A_lin += 8;
            B_lin += 8;

            ASM_MUL_ADD_8(sum1, sum2, A_lin, B_lin);
            A_lin += 8;
            B_lin += 8;

            ASM_MUL_ADD_8(sum1, sum2, A_lin, B_lin);
            A_lin += 8;
            B_lin += 8;

            ASM_MUL_ADD_8(sum1, sum2, A_lin, B_lin);
            A_lin += 8;
            B_lin += 8;

            A_lin -= M;
            B_lin -= M;
            B_lin += N;

            // Perform C[i*N + j] = sum1 + sum2
            C += i*N + j;
            ASM_ACCUMULATE_S1_S2(sum1, sum2, C);
            C -= i*N + j;
        }
        A_lin += N;
    }
}

// Compute A*B^t, where A is a N1 x M block and B is a M x N2 block
// All A, B and C are subblocks of matrices that are N x N
// The output block will be a N1 x N2 matrix and will be added
//  to the block that starts at C
// Use this function for M < 64
static inline void mul_tr_block(register double * restrict C, register double * restrict A, 
        register double * restrict B, register size_t N1, register size_t M,
        register size_t N2, register size_t N) {
    size_t i;
    register size_t j;
    register size_t k;

    register double *A_lin = A;
    // For each line in A
    for (i = 0; i < N1; i++) {
        register double *B_lin = B;
        // For each line in B
        for (j = 0; j < N2; j++) {
            register double sum1 = 0;
            register double sum2 = 0;
            register double sum3 = 0;
            register double sum4 = 0;

            for (k = 0; k + 4 <= M; k += 4) {
                sum1 += *A_lin * (*B_lin);
                sum2 += *(A_lin + 1) * (*(B_lin + 1));
                sum3 += *(A_lin + 2) * (*(B_lin + 2));
                sum4 += *(A_lin + 3) * (*(B_lin + 3));

                A_lin += 4;
                B_lin += 4;
            }

            while (k < M) {
                sum1 += *A_lin * (*B_lin);
                A_lin++;
                B_lin++;
                k++;
            }

            A_lin -= M;
            B_lin -= M;
            B_lin += N;

            C[i*N + j] += sum1 + sum2 + sum3 + sum4;
        }
        A_lin += N;
    }
}

#ifndef BL_Y1
    #define BL_Y1 32//12//12/*64*/
#endif
#ifndef BL_X
    #define BL_X 64//8 /*64*/
#endif
#ifndef BL_Y2
    #define BL_Y2 32//8/*32*/
#endif

#if BL_X == 64
    #define mul_tr_block_BL_X mul_tr_block_64
#else
    #define mul_tr_block_BL_X mul_tr_block
#endif

void mul_tr(double * restrict C, register double * A, double * B, 
        register size_t n) {
    size_t i;
    size_t j;
    size_t k;

    memset(C, 0, n*n*sizeof(*C));
    
    for (i = 0; i < n; i += BL_Y1) {
        for (j = 0; j < n; j += BL_Y2) {
            for (k = 0; k < n; k += BL_X) {
                // Multiply block (i, k) from A with block (j, k) from B
                // And store in block (i, j) from C
                // Where block A is BL_Y1 x BL_X
                // Block B is BL_Y2 x BL_X
                // And block C is BL_Y1 x BL_Y2

                if (n - k >= BL_X) {
                    // Use block multiplication with BL_X line size
                    mul_tr_block_BL_X(C + i*n + j, A + i*n + k, B + j*n + k, 
                            min(BL_Y1, n - i), 
                            BL_X, 
                            min(BL_Y2, n -j), n);
                // If the block length is smaller than BL_X
                } else {
                    // Use generic block multiplication
                    mul_tr_block(C + i*n + j, A + i*n + k, B + j*n + k, 
                            min(BL_Y1, n - i), 
                            min(BL_X, n - k), 
                            min(BL_Y2, n -j), n);
                }
            }
        }
    }
}

// Compute B = A*A^t
void mul_tr_self(double * restrict B, double * restrict A, size_t n) {
    size_t i;
    size_t j;
    size_t k;
    memset(B, 0, n*n*sizeof(*B));
    
    for (i = 0; i < n; i += BL_Y1) {
        // Take only the blocks that contain elements above the diagonal
        // The others will be deduced by symmetry
        for (j = (i / BL_Y2) * BL_Y2; j < n; j += BL_Y2) {
            for (k = 0; k < n; k += BL_X) {
                // Multiply block (i, k) from A with block (j, k) from B
                // And store in block (i, j) from C
                // Where block A is BL_Y1 x BL_X
                // Block B is BL_Y2 x BL_X
                // And block C is BL_Y1 x BL_Y2
                if (n - k >= BL_X) {
                    mul_tr_block_BL_X(B + i*n + j, A + i*n + k, A + j*n + k, 
                            min(BL_Y1, n - i), 
                            min(BL_X, n - k), 
                            min(BL_Y2, n -j), n);
                } else {
                    mul_tr_block(B + i*n + j, A + i*n + k, A + j*n + k, 
                            min(BL_Y1, n - i), 
                            min(BL_X, n - k), 
                            min(BL_Y2, n -j), n);
                }
            }
        }
    }

    // B is symmetric
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            B[j*n + i] = B[i*n + j];
        }
    }
}
