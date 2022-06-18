/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"
#include <string.h>
#include <emmintrin.h>

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

#define BLOCK_SIZE 40

// Compute B = A * A^t
// Knowing A is upper triangular
void mul_tr_uptri(double * restrict B, double * restrict A, register size_t n) {
    register size_t i;
    register size_t j;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            register size_t start = max(i, j);
            B[i*n + j] = 0;
            // register __m128d sum1 = _mm_set_pd(0.0, 0.0);
            register double sum1 = 0.0;
            register double sum2 = 0.0;

            register double *A_lin_1 = &A[i*n + start];
            register double *A_lin_2 = &A[j*n + start];
            register double *A_fin = &A[i*n + n];
            while (A_lin_1 < A_fin) {
                sum1 += *A_lin_1 * (*A_lin_2);
                sum2 += *(A_lin_1 + 1) * (*(A_lin_2 + 1));

                A_lin_1 += 2;
                A_lin_2 += 2;
            }
            B[i*n + j] = sum1 + sum2;

            // The output is symetrical
            B[j*n + i] = B[i*n + j];
        }
    }
    
}

// Assembly routine for vector-vector multiplication of size 8
#define ASM_MUL_ADD_8(sum1, sum2, A, B) __asm__ (   \
                            "movapd (%2), %%xmm1\n\t"\
                            "movapd (%3), %%xmm2\n\t"\
                            "mulpd %%xmm2, %%xmm1\n\t"\
                            "addpd %0, %%xmm1\n\t"\
                            "movapd %%xmm1, %0\n\t"\
                            \
                            "movapd 16(%2), %%xmm1\n\t"\
                            "movapd 16(%3), %%xmm2\n\t"\
                            "mulpd %%xmm2, %%xmm1\n\t"\
                            "addpd %1, %%xmm1\n\t"\
                            "movapd %%xmm1, %1\n\t"\
                            \
                            "movapd 32(%2), %%xmm1\n\t"\
                            "movapd 32(%3), %%xmm2\n\t"\
                            "mulpd %%xmm2, %%xmm1\n\t"\
                            "addpd %0, %%xmm1\n\t"\
                            "movapd %%xmm1, %0\n\t"\
                            \
                            "movapd 48(%2), %%xmm1\n\t"\
                            "movapd 48(%3), %%xmm2\n\t"\
                            "mulpd %%xmm2, %%xmm1\n\t"\
                            "addpd %1, %%xmm1\n\t"\
                            "movapd %%xmm1, %1\n\t"\
                            \
                            : "+x"(sum1), "+x"(sum2)\
                            : "rm"(A), "rm"(B)\
                            : "%xmm1", "%xmm2")
                            /*
                            */

// Compute A*B^t, where A is a N1 x M block and B is a M x N2 block
// All A, B and C are subblocks of matrices that are N x N
// The output block will be a N1 x N2 matrix and will be added
//  to the block that starts at C
static inline void mul_tr_block_64(register double * restrict C, register double * restrict A, 
        register double * restrict B, register size_t N1, register size_t M,
        register size_t N2, register size_t N) {
    register size_t i;
    register size_t j;

    register double *A_lin = A;
    for (i = 0; i < N1; i++) {
        register double *B_lin = B;
        for (j = 0; j < N2; j++) {
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

            C += i*N + j;
            __asm__ (   
                        "addpd %1, %0\n\t"
                        //"addpd %3, %2\n\t"
                        // "addpd %2, %0\n\t"
                        "pxor %1, %1\n\t"
                        "haddpd %1, %0\n\t"
                        "movlpd (%2), %1\n\t"
                        "addpd %1, %0\n\t"
                        "movlpd %0, (%2)\n\t"

                        : "+x"(sum1), "+x"(sum2), "+rm"(C)
                        :
                        :);
            C -= i*N + j;
        }
        A_lin += N;
    }
}

static inline void mul_tr_block(register double * restrict C, register double * restrict A, 
        register double * restrict B, register size_t N1, register size_t M,
        register size_t N2, register size_t N) {
    register size_t i;
    register size_t j;
    register size_t k;

    register double *A_lin = A;
    for (i = 0; i < N1; i++) {
        register double *B_lin = B;
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

void mul_tr(double * restrict C, register double * A, double * B, int N) {
    size_t i;
    size_t j;
    size_t k;
    const register size_t n = N;

    memset(C, 0, n*n*sizeof(*C));
    
    for (i = 0; i < n; i += BL_Y1) {
        for (j = 0; j < n; j += BL_Y2) {
            for (k = 0; k < n; k += BL_X) {
                // Multiply block (i, k) from A with block (j, k) from B
                // And store in block (i, j) from C
                // Where block A is BL_Y1 * BL_X
                // Block B is BL_Y2 * BL_X
                // And block C is BL_Y1 * BL_Y2

                if (n - k >= BL_X) {
                    mul_tr_block_BL_X(C + i*n + j, A + i*n + k, B + j*n + k, 
                            min(BL_Y1, n - i), 
                            min(BL_X, n - k), 
                            min(BL_Y2, n -j), N);
                } else {
                    mul_tr_block(C + i*n + j, A + i*n + k, B + j*n + k, 
                            min(BL_Y1, n - i), 
                            min(BL_X, n - k), 
                            min(BL_Y2, n -j), N);
                }
                // mul_tr_m8_a4_b4(C + i*n + j, A + i*n + k, B + j*n + k, N);

            }
        }
    }
}

// Compute B = A*A^t + C
void mult_mat_trans_add(double *B, double *A, double *C, size_t n) {
    register size_t i;
    register size_t j;
    const register size_t nsq = n * n;
    register size_t ni = 0;

    memcpy(B, C, n*n*sizeof(*B));

    register double *A_lin_1 = A;
    register double *A_lin_2 = A + n;
    register double *B_lin = B;
    register double *B_col = B;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            register double sum1 = 0;
            register double sum2 = 0;
            register double *A_fin = A_lin_1 + n;
            while (A_lin_1 < A_fin) {
                sum1 += *A_lin_1 * (*A_lin_2);
                sum2 += *(A_lin_1 + 1) * (*(A_lin_2 + 1));
                A_lin_1 += 2;
                A_lin_2 += 2;
            }
            A_lin_1 -= n;
            *B_lin += sum1 + sum2;
            // The output is symmetrical
            *B_col += sum1 + sum2;

            B_lin++;
            B_col += n;
            // C_lin++;
            // C_col += n;
        }
        B_lin += i + 2;
        B_col += 2*n + 1 - nsq + ni;
        // C_lin += i + 2;
        // C_col += 2*n + 1 - nsq + ni;

        A_lin_1 += n;
        A_lin_2 -= nsq - ni;
        ni += n;
        A_lin_2 += 2*n;
    }

    A_lin_1 = A;
    A_lin_2 = A;
    B_lin = B;
    B_col = B;
    // C_lin = C;
    // C_col = C;
    ni = 0;

    for (i = 0; i < n; i++) {
        register double sum1 = 0;
        register double sum2 = 0;
        register double *A_fin = A_lin_1 + n;

        while (A_lin_1 < A_fin) {
            sum1 += *A_lin_1 * (*A_lin_2);
            sum2 += *(A_lin_1 + 1) * (*(A_lin_2 + 1));
            A_lin_1 += 2;
            A_lin_2 += 2;
        }
        // Assign the result
        *B_lin += sum1 + sum2;//+ *C_lin;

        B_lin += n + 1;
        B_col += n + 1;
        // C_lin += n + 1;
        // C_col += n + 1;

        ni += n;
    }
}

// Set A = A^t
void transpose(double *A, int n) {
    register size_t i, j, k;
    // Transpose the blocks that are on the main diagonal
    register size_t i_start = 0;
    register size_t n_blocks = n / BLOCK_SIZE;
    // Iterates though the blocks
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
        i_start += BLOCK_SIZE;
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

// Compute B = A*A^t
void mult_mat_trans_2(double * restrict B, double * restrict A, size_t n) {
    size_t i;
    size_t j;
    size_t k;
    memset(B, 0, n*n*sizeof(*B));
    
    for (i = 0; i < n; i += BL_Y1) {
        for (j = (i / BL_Y2) * BL_Y2; j < n; j += BL_Y2) {
            for (k = 0; k < n; k += BL_X) {
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

// Compute B = A*A^t
void mult_mat_trans(double * restrict B, double * restrict A, size_t n) {
    register size_t i;
    register size_t j;
    const register size_t nsq = n * n;
    register size_t ni = 0;

    register double *A_lin_1 = A;
    register double *A_lin_2 = A;
    register double *B_lin = B;
    register double *B_col = B;
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            register __m128d sum = _mm_set_pd(0.0, 0.0);
            register double *A_fin = A_lin_1 + n;
            while (A_lin_1 < A_fin) {
                __asm__ (
                            "movapd (%1), %%xmm1\n\t"
                            "movapd (%2), %%xmm2\n\t"
                            "mulpd %%xmm2, %%xmm1\n\t"
                            "addpd %0, %%xmm1\n\t"
                            "movapd %%xmm1, %0\n\t"

                            : "+x"(sum)
                            : "r"(A_lin_1), "r"(A_lin_2)
                            : "%xmm1", "%xmm2");
                A_lin_1 += 2;
                A_lin_2 += 2;
            }
            A_lin_1 -= n;
            __asm__ volatile (   
                        "pxor %%xmm0, %%xmm0\n\t"
                        "haddpd %2, %%xmm0\n\t"
                        "movhpd %%xmm0, (%0)\n\t"
                        "movhpd %%xmm0, (%1)\n\t"
                        "add $8, %0\n\t"
                        "mov %3, %%r15\n\t"
                        "sal $3, %%r15\n\t"
                        "add %%r15, %1\n\t"

                        : "+rm"(B_lin), "+rm"(B_col)
                        : "x"(sum), "r"(n)
                        : "%xmm0", "%r15");
            // *B_lin = sum1 + sum2;
            // The output is symmetrical
            // B_lin++;
            // B_col += n;
        }
        B_lin += i + 1;
        B_col += n + 1 - nsq + ni;

        A_lin_1 += n;
        A_lin_2 -= nsq - ni;
        ni += n;
        A_lin_2 += n;
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
    // C is symetrical
    mul_tr(D, B, C, n);

    // Tranpose: B = B^t
    transpose(B, n);

    // Compute C = B * B^t
    mult_mat_trans_2(C, B, n);

    // Compute C = C + D
    add_mat(C, D, n);

    free(D);
    
    return C;
}
