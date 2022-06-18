/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"
#include <string.h>

#define min(a, b) ((a) < (b) ? (a) : (b))

#define BLOCK_SIZE 40

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

void mul_tr_block_2(double *C, double *A, double *B, int N) {
    register size_t i;
    register size_t j;
    register size_t k;
    register size_t n = N;
    const register size_t block_size = 40;

    const register size_t nb = n * block_size;

    memset(C, 0, n*n*sizeof(*C));
    register double *C_lin = C;
    for (register size_t bl = 0; bl < n; bl += block_size) {
        for (register size_t br = 0; br < n; br += block_size) {
            register double *A_lin = A + bl * n + br;
            // C_lin = C + bl*n;
            for (register size_t bl2 = 0; bl2 < n; bl2 += block_size) {
                // register double *A_lin = A + bl * n + br;
                register double *B_lin = B + bl2 * n + br;
                // C_lin = C + bl*n + bl2;
                for (i = 0; i < block_size; i++) {
                    // register double *B_lin = B + bl2 * n + br;
                    for (j = 0; j < block_size; j++) {
                        register double sum1 = 0;
                        register double sum2 = 0;
                        register double sum3 = 0;
                        register double sum4 = 0;

                        // A_lin will be A + (i + bl)*n + br
                        // B_lin will be B + (j + bl2)*n + br
                        for (k = 0; k < block_size; k += 4) {
                            sum1 += *A_lin * (*B_lin);
                            sum2 += *(A_lin + 1) * (*(B_lin + 1));
                            sum3 += *(A_lin + 2) * (*(B_lin + 2));
                            sum4 += *(A_lin + 3) * (*(B_lin + 3));

                            A_lin += 4;
                            B_lin += 4;
                        }
                        A_lin -= block_size;//n;
                        B_lin += n - block_size;//Weird
                        *(C_lin) += sum1 + sum2 + sum3 + sum4;
                        C_lin++;
                    }
                    C_lin += n - block_size;
                    B_lin -= nb;
                    A_lin += n;
                }
                C_lin -= nb;
                C_lin += block_size;
                A_lin -= nb;
            }
            C_lin -= n;
        }
        C_lin += nb;
    }
}

void mul_tr_block_1(double *C, double *A, double *B, int N) {
    register size_t i;
    register size_t j;
    register size_t k;
    register size_t n = N;
    const register size_t block_size = n / 2;

    memset(C, 0, n*n*sizeof(*C));
    for (i = 0; i < n; i++) {
        for (register size_t bl = 0; bl < n; bl += block_size) {
            for (j = 0; j < n; j++) {
                register double sum1 = 0;
                register double sum2 = 0;
                register double sum3 = 0;
                register double sum4 = 0;

                for (k = 0; k < block_size; k += 4) {
                    sum1 += A[i*n + bl + k] * B[j*n + bl + k];
                    sum2 += A[i*n + bl + k + 1] * B[j*n + bl + k + 1];
                    sum3 += A[i*n + bl + k + 2] * B[j*n + bl + k + 2];
                    sum4 += A[i*n + bl + k + 3] * B[j*n + bl + k + 3];
                }
                C[i*n + j] += sum1 + sum2 + sum3 + sum4;
            }
        }
    }
}

// Copmute C = A * B^t
void mul_tr(register double *C, register double *A, register double *B, int N) {
    register size_t n = N;
    register size_t i, j, k;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            register double sum = 0;
            for (k = 0; k < n; k++) {
                sum += A[i*n + k] * B[j*n + k];
            }
            C[i*n + j] = sum;
        }
    }
}

void mul_tr_2(double * restrict C, register double * restrict A, double * restrict B, int N) {
    register size_t i;
    const register size_t n = N;
    const register size_t nsq = n*n;

    register double *A_lin = A;
    register double *C_lin = C;

    for (i = 0; i < n; i++) {
        register double *B_lin = B;
        register double *B_fin = B + nsq;
        while (B_lin < B_fin) {
            register double sum1 = 0;
            register double sum2 = 0;
            register double sum3 = 0;
            register double sum4 = 0;

            register double *A_fin = A_lin + n;
            while (A_lin < A_fin) {
                sum1 += *A_lin * (*B_lin);
                sum2 += *(A_lin + 1) * (*(B_lin + 1));
                sum3 += *(A_lin + 2) * (*(B_lin + 2));
                sum4 += *(A_lin + 3) * (*(B_lin + 3));

                A_lin += 4;
                B_lin += 4;
            }
            A_lin -= n;
            *(C_lin) = sum1 + sum2 + sum3 + sum4;
            C_lin++;
        }
        A_lin += n;
    }
}

#define BL_SIZE 16

void mul_tr_3(double * restrict C, register double * restrict A, double * restrict B, int N) {
    register size_t i;
    register size_t j;
    register size_t k;
    const register size_t n = N;

    memset(C, 0, n*n*sizeof(*C));

    transpose(B, n);

    register size_t i2, k2;
    register double *restrict rres;
    register double *restrict rmul1;
    register double *restrict rmul2;

    register size_t in = 0;

    for (i = 0; i < n; i += 2*BL_SIZE) {
        for (j = 0; j < n; j += 12*BL_SIZE) {
            register size_t kn = 0;
            for (k = 0; k < n; k += BL_SIZE) {
                rres = &C[in + j];
                rmul1 = &A[in + k];
                for (i2 = 0; i2 < 2*BL_SIZE; i2 += 2, rres += 2*n, rmul1 += 2*n) { rmul2 = &B[kn + j];
                    for (k2 = 0; k2 < BL_SIZE; k2 += 4, rmul2 += 4*n) {
                        register double rm = rmul1[k2];
                        register double rm2 = rmul1[k2 + 1];
                        register double rm3 = rmul1[k2 + 2];
                        register double rm4 = rmul1[k2 + 3];

                        const register double rm5 = rmul1[k2 + n];
                        const register double rm6 = rmul1[k2 + n + 1];
                        const register double rm7 = rmul1[k2 + n + 2];
                        const register double rm8 = rmul1[k2 + n + 3];

                        register size_t j2n = n;
                        register size_t j3n = j2n + n;
                        register size_t j4n = j3n + n;
                        /*
                        register size_t j5n = j4n + n;
                        register size_t j6n = j5n + n;
                        register size_t j7n = j6n + n;
                        register size_t j8n = j7n + n;
                        */

                        for (register size_t j2 = 0; j2 < 12*BL_SIZE; j2 += 4, j2n += 4, j3n += 4, j4n += 4/*, j5n += 4, j6n += 4, j7n += 4, j8n += 4*/) {
                            rres[j2] += rm * rmul2[j2] + rm2 * rmul2[j2n] + rm3 * rmul2[j3n] + rm4 * rmul2[j4n];// + rm5 * rmul2[j5n] + rm6 * rmul2[j6n] + rm7 * rmul2[j7n] + rm8 * rmul2[j8n];
                            rres[j2 + 1] += rm * rmul2[j2 + 1] + rm2 * rmul2[j2n + 1] + rm3 * rmul2[j3n + 1] + rm4 * rmul2[j4n + 1];// + rm5 * rmul2[j5n + 1] + rm6 * rmul2[j6n + 1] + rm7 * rmul2[j7n + 1] + rm8 * rmul2[j8n + 1];
                            rres[j2 + 2] += rm * rmul2[j2 + 2] + rm2 * rmul2[j2n + 2] + rm3 * rmul2[j3n + 2] + rm4 * rmul2[j4n + 2];// + rm5 * rmul2[j5n + 2] + rm6 * rmul2[j6n + 2] + rm7 * rmul2[j7n + 2] + rm8 * rmul2[j8n + 2];
                            rres[j2 + 3] += rm * rmul2[j2 + 3] + rm2 * rmul2[j2n + 3] + rm3 * rmul2[j3n + 3] + rm4 * rmul2[j4n + 3];// + rm5 * rmul2[j5n + 3] + rm6 * rmul2[j6n + 3] + rm7 * rmul2[j7n + 3] + rm8 * rmul2[j8n + 3];

                            rres[j2n] += rm5 * rmul2[j2] + rm6 * rmul2[j2n] + rm7 * rmul2[j3n] + rm8 * rmul2[j4n];// + rm5 * rmul2[j5n] + rm6 * rmul2[j6n] + rm7 * rmul2[j7n] + rm8 * rmul2[j8n];
                            rres[j2n + 1] += rm5 * rmul2[j2 + 1] + rm6 * rmul2[j2n + 1] + rm7 * rmul2[j3n + 1] + rm8 * rmul2[j4n + 1];// + rm5 * rmul2[j5n + 1] + rm6 * rmul2[j6n + 1] + rm7 * rmul2[j7n + 1] + rm8 * rmul2[j8n + 1];
                            rres[j2n + 2] += rm5 * rmul2[j2 + 2] + rm6 * rmul2[j2n + 2] + rm7 * rmul2[j3n + 2] + rm8 * rmul2[j4n + 2];// + rm5 * rmul2[j5n + 2] + rm6 * rmul2[j6n + 2] + rm7 * rmul2[j7n + 2] + rm8 * rmul2[j8n + 2];
                            rres[j2n + 3] += rm5 * rmul2[j2 + 3] + rm6 * rmul2[j2n + 3] + rm7 * rmul2[j3n + 3] + rm8 * rmul2[j4n + 3];// + rm5 * rmul2[j5n + 3] + rm6 * rmul2[j6n + 3] + rm7 * rmul2[j7n + 3] + rm8 * rmul2[j8n + 3];
                        }
                    }
                }
                kn += n*BL_SIZE;
            }
        }
        in += 2*n*BL_SIZE;
    }
}

void mul_tr_m4_a4_b4(register double * restrict C, register double * restrict A, 
        register double * restrict B, register size_t N) {
    register double a00 = A[0];

    register double a10 = A[N];
    register double a20 = A[2*N];
    register double a30 = A[3*N];

    register size_t N3 = N*3;

#define b11 (B[N + 1])
#define b22 (B[2*N + 2])
#define b33 (B[3*N + 3])

    register double b00 = B[0];
    register double b10 = B[N];
    register double b20 = B[2*N];
    register double b30 = B[N3];

    C[0] += a00*b00 + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    C[1] += a00*b10 + A[1]*b11 + A[2]*B[1*N + 2] + A[3]*B[1*N + 3];
    C[2] += a00*b20 + A[1]*B[2*N + 1] + A[2]*b22 + A[3]*B[2*N + 3];
    C[3] += a00*b30 + A[1]*B[N3 + 1] + A[2]*B[N3 + 2] + A[3]*b33;

    A += N;
    C += N;
    C[0] += a10*b00 + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    C[1] += a10*b10 + A[1]*b11 + A[2]*B[1*N + 2] + A[3]*B[1*N + 3];
    C[2] += a10*b20 + A[1]*B[2*N + 1] + A[2]*b22 + A[3]*B[2*N + 3];
    C[3] += a10*b30 + A[1]*B[N3 + 1] + A[2]*B[N3 + 2] + A[3]*b33;

    A += N;
    C += N;
    C[0] += a20*b00 + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    C[1] += a20*b10 + A[1]*b11 + A[2]*B[1*N + 2] + A[3]*B[1*N + 3];
    C[2] += a20*b20 + A[1]*B[2*N + 1] + A[2]*b22 + A[3]*B[2*N + 3];
    C[3] += a20*b30 + A[1]*B[N3 + 1] + A[2]*B[N3 + 2] + A[3]*b33;

    A += N;
    C += N;
    C[0] += a30*b00 + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    C[1] += a30*b10 + A[1]*b11 + A[2]*B[1*N + 2] + A[3]*B[1*N + 3];
    C[2] += a30*b20 + A[1]*B[2*N + 1] + A[2]*b22 + A[3]*B[2*N + 3];
    C[3] += a30*b30 + A[1]*B[N3 + 1] + A[2]*B[N3 + 2] + A[3]*b33;
}

#include <emmintrin.h>
void mul_tr_m8_a4_b4(register double * restrict C, register double * restrict A, 
        register double * restrict B, register size_t N) {
    register double a00 = A[0];

    register double a10 = A[N];
    register double a20 = A[2*N];
    register double a30 = A[3*N];

    register size_t N3 = N*3;

#define b11 (B[N + 1])
#define b22 (B[2*N + 2])
#define b33 (B[3*N + 3])

    register double b00 = B[0];
    register double b10 = B[N];
    register double b20 = B[2*N];
    register double b30 = B[N3];

    register __m128d sum1 = _mm_set_pd(0.0, 0.0);
    // register __m128d sum2;
    __asm__ (   
                "movq %2, %%xmm0\n\t"
                "movq %3, %%xmm1\n\t"
                "movhpd 8(%4), %%xmm0\n\t"
                "movhpd 8(%5), %%xmm1\n\t"

                "movapd 16(%4), %%xmm2\n\t"
                "movapd 16(%5), %%xmm3\n\t"

                "mulpd %%xmm1, %%xmm0\n\t"
                "mulpd %%xmm3, %%xmm2\n\t"
                "haddpd %%xmm2, %%xmm0\n\t"

                "movapd 32(%4), %%xmm1\n\t"
                "movapd 32(%5), %%xmm2\n\t"
                "mulpd %%xmm2, %%xmm1\n\t"
                "haddpd %%xmm1, %%xmm0\n\t"

                "movapd 48(%4), %%xmm1\n\t"
                "movapd 48(%5), %%xmm2\n\t"
                "mulpd %%xmm2, %%xmm1\n\t"
                "haddpd %%xmm1, %%xmm0\n\t"

                "pxor %%xmm1, %%xmm1\n\t"
                "vhaddpd %%xmm1, %%xmm0, %0\n\t"

                "movhpd (%1), %0\n\t"
                "pxor %%xmm1, %%xmm1\n\t"
                "haddpd %%xmm1, %0\n\t"
                "movlpd %0, (%1)\n\t"

                : "+x"(sum1), "+rm"(C)
                : "r"(a00), "r"(b00), "r"(A), "r"(B)
                : "%xmm0", "%xmm1", "%xmm2", "%xmm3");

    B += N;
    __asm__ (   
                "movq %2, %%xmm0\n\t"
                "movq %3, %%xmm1\n\t"
                "movhpd 8(%4), %%xmm0\n\t"
                "movhpd 8(%5), %%xmm1\n\t"

                "movapd 16(%4), %%xmm2\n\t"
                "movapd 16(%5), %%xmm3\n\t"

                "mulpd %%xmm1, %%xmm0\n\t"
                "mulpd %%xmm3, %%xmm2\n\t"
                "haddpd %%xmm2, %%xmm0\n\t"

                "movapd 32(%4), %%xmm1\n\t"
                "movapd 32(%5), %%xmm2\n\t"
                "mulpd %%xmm2, %%xmm1\n\t"
                "haddpd %%xmm1, %%xmm0\n\t"

                "movapd 48(%4), %%xmm1\n\t"
                "movapd 48(%5), %%xmm2\n\t"
                "mulpd %%xmm2, %%xmm1\n\t"
                "haddpd %%xmm1, %%xmm0\n\t"

                "pxor %%xmm1, %%xmm1\n\t"
                "vhaddpd %%xmm1, %%xmm0, %0\n\t"

                "movhpd 8(%1), %0\n\t"
                "pxor %%xmm1, %%xmm1\n\t"
                "haddpd %%xmm1, %0\n\t"
                // Load in C[1]
                "movlpd %0, 8(%1)\n\t"

                : "+x"(sum1), "+rm"(C)
                : "r"(a00), "r"(b10), "r"(A), "r"(B)
                : "%xmm0", "%xmm1", "%xmm2", "%xmm3");

    B += N;
    __asm__ (   
                "movq %2, %%xmm0\n\t"
                "movq %3, %%xmm1\n\t"
                "movhpd 8(%4), %%xmm0\n\t"
                "movhpd 8(%5), %%xmm1\n\t"

                "movapd 16(%4), %%xmm2\n\t"
                "movapd 16(%5), %%xmm3\n\t"

                "mulpd %%xmm1, %%xmm0\n\t"
                "mulpd %%xmm3, %%xmm2\n\t"
                "haddpd %%xmm2, %%xmm0\n\t"

                "movapd 32(%4), %%xmm1\n\t"
                "movapd 32(%5), %%xmm2\n\t"
                "mulpd %%xmm2, %%xmm1\n\t"
                "haddpd %%xmm1, %%xmm0\n\t"

                "movapd 48(%4), %%xmm1\n\t"
                "movapd 48(%5), %%xmm2\n\t"
                "mulpd %%xmm2, %%xmm1\n\t"
                "haddpd %%xmm1, %%xmm0\n\t"

                "pxor %%xmm1, %%xmm1\n\t"
                "vhaddpd %%xmm1, %%xmm0, %0\n\t"

                // Get previous value from C[2]
                "movhpd 16(%1), %0\n\t"
                "pxor %%xmm1, %%xmm1\n\t"
                "haddpd %%xmm1, %0\n\t"
                // Load in C[2]
                "movlpd %0, 16(%1)\n\t"

                : "+x"(sum1), "+rm"(C)
                : "r"(a00), "r"(b20), "r"(A), "r"(B)
                : "%xmm0", "%xmm1", "%xmm2", "%xmm3");

    B += N;
    __asm__ (   
                "movq %2, %%xmm0\n\t"
                "movq %3, %%xmm1\n\t"
                "movhpd 8(%4), %%xmm0\n\t"
                "movhpd 8(%5), %%xmm1\n\t"

                "movapd 16(%4), %%xmm2\n\t"
                "movapd 16(%5), %%xmm3\n\t"

                "mulpd %%xmm1, %%xmm0\n\t"
                "mulpd %%xmm3, %%xmm2\n\t"
                "haddpd %%xmm2, %%xmm0\n\t"

                "movapd 32(%4), %%xmm1\n\t"
                "movapd 32(%5), %%xmm2\n\t"
                "mulpd %%xmm2, %%xmm1\n\t"
                "haddpd %%xmm1, %%xmm0\n\t"

                "movapd 48(%4), %%xmm1\n\t"
                "movapd 48(%5), %%xmm2\n\t"
                "mulpd %%xmm2, %%xmm1\n\t"
                "haddpd %%xmm1, %%xmm0\n\t"

                "pxor %%xmm1, %%xmm1\n\t"
                "vhaddpd %%xmm1, %%xmm0, %0\n\t"

                // Get previous value from C[2]
                "movhpd 24(%1), %0\n\t"
                "pxor %%xmm1, %%xmm1\n\t"
                "haddpd %%xmm1, %0\n\t"
                // Load in C[2]
                "movlpd %0, 24(%1)\n\t"

                : "+x"(sum1), "+rm"(C)
                : "r"(a00), "r"(b30), "r"(A), "r"(B)
                : "%xmm0", "%xmm1", "%xmm2", "%xmm3");


    // C[0] += a00*b00 + A[1]*B[1] + A[2]*B[2] + A[3]*B[3] + A[4]*B[4] + A[5]*B[5] + A[6]*B[6] + A[7]*B[7];
    // C[1] += a00*b10 + A[1]*b11 + A[2]*B[1*N + 2] + A[3]*B[1*N + 3] + A[4]*B[1*N + 4] + A[5]*B[1*N + 5] + A[6]*B[1*N + 6] + A[7]*B[1*N + 7];
    // C[2] += a00*b20 + A[1]*B[2*N + 1] + A[2]*b22 + A[3]*B[2*N + 3] + A[4]*B[2*N + 4] + A[5]*B[2*N + 5] + A[6]*B[2*N + 6] + A[7]*B[2*N + 7];
    // C[3] += a00*b30 + A[1]*B[N3 + 1] + A[2]*B[N3 + 2] + A[3]*b33 + A[4]*B[N3 + 4] + A[5]*B[N3 + 5] + A[6]*B[N3 + 6] + A[7]*B[N3 + 7];

    B -= N3;

    A += N;
    C += N;
    __asm__ (   
                "movq %2, %%xmm0\n\t"
                "movq %3, %%xmm1\n\t"
                "movhpd 8(%4), %%xmm0\n\t"
                "movhpd 8(%5), %%xmm1\n\t"

                "movapd 16(%4), %%xmm2\n\t"
                "movapd 16(%5), %%xmm3\n\t"

                "mulpd %%xmm1, %%xmm0\n\t"
                "mulpd %%xmm3, %%xmm2\n\t"
                "haddpd %%xmm2, %%xmm0\n\t"

                "movapd 32(%4), %%xmm1\n\t"
                "movapd 32(%5), %%xmm2\n\t"
                "mulpd %%xmm2, %%xmm1\n\t"
                "haddpd %%xmm1, %%xmm0\n\t"

                "movapd 48(%4), %%xmm1\n\t"
                "movapd 48(%5), %%xmm2\n\t"
                "mulpd %%xmm2, %%xmm1\n\t"
                "haddpd %%xmm1, %%xmm0\n\t"

                "pxor %%xmm1, %%xmm1\n\t"
                "vhaddpd %%xmm1, %%xmm0, %0\n\t"

                // Get previous value from C[0]
                "movhpd (%1), %0\n\t"
                "pxor %%xmm1, %%xmm1\n\t"
                "haddpd %%xmm1, %0\n\t"
                // Load in C[0]
                "movlpd %0, (%1)\n\t"

                : "+x"(sum1), "+rm"(C)
                : "r"(a10), "r"(b00), "r"(A), "r"(B)
                : "%xmm0", "%xmm1", "%xmm2", "%xmm3");

    B += N;
    __asm__ (   
                "movq %2, %%xmm0\n\t"
                "movq %3, %%xmm1\n\t"
                "movhpd 8(%4), %%xmm0\n\t"
                "movhpd 8(%5), %%xmm1\n\t"

                "movapd 16(%4), %%xmm2\n\t"
                "movapd 16(%5), %%xmm3\n\t"

                "mulpd %%xmm1, %%xmm0\n\t"
                "mulpd %%xmm3, %%xmm2\n\t"
                "haddpd %%xmm2, %%xmm0\n\t"

                "movapd 32(%4), %%xmm1\n\t"
                "movapd 32(%5), %%xmm2\n\t"
                "mulpd %%xmm2, %%xmm1\n\t"
                "haddpd %%xmm1, %%xmm0\n\t"

                "movapd 48(%4), %%xmm1\n\t"
                "movapd 48(%5), %%xmm2\n\t"
                "mulpd %%xmm2, %%xmm1\n\t"
                "haddpd %%xmm1, %%xmm0\n\t"

                "pxor %%xmm1, %%xmm1\n\t"
                "vhaddpd %%xmm1, %%xmm0, %0\n\t"

                // Get previous value from C[1]
                "movhpd 8(%1), %0\n\t"
                "pxor %%xmm1, %%xmm1\n\t"
                "haddpd %%xmm1, %0\n\t"
                // Load in C[1]
                "movlpd %0, 8(%1)\n\t"

                : "+x"(sum1), "+rm"(C)
                : "r"(a10), "r"(b10), "r"(A), "r"(B)
                : "%xmm0", "%xmm1", "%xmm2", "%xmm3");

    B -= N;

    //C[0] += a10*b00 + A[1]*B[1] + A[2]*B[2] + A[3]*B[3] + A[4]*B[4] + A[5]*B[5] + A[6]*B[6] + A[7]*B[7];
    // C[1] += a10*b10 + A[1]*b11 + A[2]*B[1*N + 2] + A[3]*B[1*N + 3] + A[4]*B[1*N + 4] + A[5]*B[1*N + 5] + A[6]*B[1*N + 6] + A[7]*B[1*N + 7];
    C[2] += a10*b20 + A[1]*B[2*N + 1] + A[2]*b22 + A[3]*B[2*N + 3] + A[4]*B[2*N + 4] + A[5]*B[2*N + 5] + A[6]*B[2*N + 6] + A[7]*B[2*N + 7];
    C[3] += a10*b30 + A[1]*B[N3 + 1] + A[2]*B[N3 + 2] + A[3]*b33 + A[4]*B[N3 + 4] + A[5]*B[N3 + 5] + A[6]*B[N3 + 6] + A[7]*B[N3 + 7];

    A += N;
    C += N;
    C[0] += a20*b00 + A[1]*B[1] + A[2]*B[2] + A[3]*B[3] + A[4]*B[4] + A[5]*B[5] + A[6]*B[6] + A[7]*B[7];
    C[1] += a20*b10 + A[1]*b11 + A[2]*B[1*N + 2] + A[3]*B[1*N + 3] + A[4]*B[1*N + 4] + A[5]*B[1*N + 5] + A[6]*B[1*N + 6] + A[7]*B[1*N + 7];
    C[2] += a20*b20 + A[1]*B[2*N + 1] + A[2]*b22 + A[3]*B[2*N + 3] + A[4]*B[2*N + 4] + A[5]*B[2*N + 5] + A[6]*B[2*N + 6] + A[7]*B[2*N + 7];
    C[3] += a20*b30 + A[1]*B[N3 + 1] + A[2]*B[N3 + 2] + A[3]*b33 + A[4]*B[N3 + 4] + A[5]*B[N3 + 5] + A[6]*B[N3 + 6] + A[7]*B[N3 + 7];

    A += N;
    C += N;
    C[0] += a30*b00 + A[1]*B[1] + A[2]*B[2] + A[3]*B[3] + A[4]*B[4] + A[5]*B[5] + A[6]*B[6] + A[7]*B[7];
    C[1] += a30*b10 + A[1]*b11 + A[2]*B[1*N + 2] + A[3]*B[1*N + 3] + A[4]*B[1*N + 4] + A[5]*B[1*N + 5] + A[6]*B[1*N + 6] + A[7]*B[1*N + 7];
    C[2] += a30*b20 + A[1]*B[2*N + 1] + A[2]*b22 + A[3]*B[2*N + 3] + A[4]*B[2*N + 4] + A[5]*B[2*N + 5] + A[6]*B[2*N + 6] + A[7]*B[2*N + 7];
    C[3] += a30*b30 + A[1]*B[N3 + 1] + A[2]*B[N3 + 2] + A[3]*b33 + A[4]*B[N3 + 4] + A[5]*B[N3 + 5] + A[6]*B[N3 + 6] + A[7]*B[N3 + 7];
}

// Compute A*B^t, where A is a N1 x M block and B is a M x N2 block
// All A, B and C are subblocks of matrices that are N x N
// The output block will be a N1 x N2 matrix and will be added
//  to the block that starts at C
static inline void mul_tr_block_reg(double * restrict C, register double * restrict A, 
        register double * restrict B, size_t N1, register size_t M,
        register size_t N2, register size_t N) {
    register size_t i;
    register size_t j;
    register size_t k;

    for (j = 0; j < N2; j += 4) {
        for (i = 0; i < N1; i += 4) {
            for (k = 0; k < M; k += 8) {
                mul_tr_m8_a4_b4(C + i*N + j, A + i*N + k, B + j*N + k, N);
                // mul_tr_m4_a4_b4(C + i*N + j, A + i*N + k, B + j*N + k, N);
            }
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
    // register size_t k;

    register double *A_lin = A;
    for (i = 0; i < N1; i++) {
        register double *B_lin = B;
        for (j = 0; j < N2; j++) {
            register __m128d sum1 = _mm_set_pd(0.0, 0.0);
            register __m128d sum2 = _mm_set_pd(0.0, 0.0);
            // register double sum1 = 0.0;
            // register double sum2 = 0.0;
            // register double sum3 = 0;
            // register double sum4 = 0;
#define sum3 (sum1)
#define sum4 (sum2)

            // for (k = 0; k < M; k += 8) {
            /*
                sum1 += *A_lin * (*B_lin);
                sum2 += *(A_lin + 1) * (*(B_lin + 1));
                sum3 += *(A_lin + 2) * (*(B_lin + 2));
                sum4 += *(A_lin + 3) * (*(B_lin + 3));
                */

            /*
                sum1 += *(A_lin + 4) * (*(B_lin + 4));
                sum2 += *(A_lin + 5) * (*(B_lin + 5));
                sum3 += *(A_lin + 6) * (*(B_lin + 6));
                sum4 += *(A_lin + 7) * (*(B_lin + 7));
                */
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

                /*
                sum1 += *A_lin * (*B_lin);
                sum2 += *(A_lin + 1) * (*(B_lin + 1));
                sum3 += *(A_lin + 2) * (*(B_lin + 2));
                sum4 += *(A_lin + 3) * (*(B_lin + 3));

                sum1 += *(A_lin + 4) * (*(B_lin + 4));
                sum2 += *(A_lin + 5) * (*(B_lin + 5));
                sum3 += *(A_lin + 6) * (*(B_lin + 6));
                sum4 += *(A_lin + 7) * (*(B_lin + 7));

                A_lin += 8;
                B_lin += 8;

                sum1 += *A_lin * (*B_lin);
                sum2 += *(A_lin + 1) * (*(B_lin + 1));
                sum3 += *(A_lin + 2) * (*(B_lin + 2));
                sum4 += *(A_lin + 3) * (*(B_lin + 3));

                sum1 += *(A_lin + 4) * (*(B_lin + 4));
                sum2 += *(A_lin + 5) * (*(B_lin + 5));
                sum3 += *(A_lin + 6) * (*(B_lin + 6));
                sum4 += *(A_lin + 7) * (*(B_lin + 7));

                A_lin += 8;
                B_lin += 8;

                sum1 += *A_lin * (*B_lin);
                sum2 += *(A_lin + 1) * (*(B_lin + 1));
                sum3 += *(A_lin + 2) * (*(B_lin + 2));
                sum4 += *(A_lin + 3) * (*(B_lin + 3));

                sum1 += *(A_lin + 4) * (*(B_lin + 4));
                sum2 += *(A_lin + 5) * (*(B_lin + 5));
                sum3 += *(A_lin + 6) * (*(B_lin + 6));
                sum4 += *(A_lin + 7) * (*(B_lin + 7));

                A_lin += 8;
                B_lin += 8;

                sum1 += *A_lin * (*B_lin);
                sum2 += *(A_lin + 1) * (*(B_lin + 1));
                sum3 += *(A_lin + 2) * (*(B_lin + 2));
                sum4 += *(A_lin + 3) * (*(B_lin + 3));

                sum1 += *(A_lin + 4) * (*(B_lin + 4));
                sum2 += *(A_lin + 5) * (*(B_lin + 5));
                sum3 += *(A_lin + 6) * (*(B_lin + 6));
                sum4 += *(A_lin + 7) * (*(B_lin + 7));

                A_lin += 8;
                B_lin += 8;

                sum1 += *A_lin * (*B_lin);
                sum2 += *(A_lin + 1) * (*(B_lin + 1));
                sum3 += *(A_lin + 2) * (*(B_lin + 2));
                sum4 += *(A_lin + 3) * (*(B_lin + 3));

                sum1 += *(A_lin + 4) * (*(B_lin + 4));
                sum2 += *(A_lin + 5) * (*(B_lin + 5));
                sum3 += *(A_lin + 6) * (*(B_lin + 6));
                sum4 += *(A_lin + 7) * (*(B_lin + 7));

                A_lin += 8;
                B_lin += 8;

                sum1 += *A_lin * (*B_lin);
                sum2 += *(A_lin + 1) * (*(B_lin + 1));
                sum3 += *(A_lin + 2) * (*(B_lin + 2));
                sum4 += *(A_lin + 3) * (*(B_lin + 3));

                sum1 += *(A_lin + 4) * (*(B_lin + 4));
                sum2 += *(A_lin + 5) * (*(B_lin + 5));
                sum3 += *(A_lin + 6) * (*(B_lin + 6));
                sum4 += *(A_lin + 7) * (*(B_lin + 7));

                A_lin += 8;
                B_lin += 8;

                sum1 += *A_lin * (*B_lin);
                sum2 += *(A_lin + 1) * (*(B_lin + 1));
                sum3 += *(A_lin + 2) * (*(B_lin + 2));
                sum4 += *(A_lin + 3) * (*(B_lin + 3));

                sum1 += *(A_lin + 4) * (*(B_lin + 4));
                sum2 += *(A_lin + 5) * (*(B_lin + 5));
                sum3 += *(A_lin + 6) * (*(B_lin + 6));
                sum4 += *(A_lin + 7) * (*(B_lin + 7));

                A_lin += 8;
                B_lin += 8;
                */

#undef sum3
#undef sum4


            // }
            A_lin -= M;
            B_lin -= M;
            B_lin += N;

            C += i*N + j;
            /*
            double s1, s2, s3, s4;
            _mm_storeh_pd(&s1, sum1);
            _mm_storel_pd(&s2, sum1);
            _mm_storeh_pd(&s3, sum2);
            _mm_storel_pd(&s4, sum2);
            *C += s1 + s2 + s3 + s4;
            */
            // *C += sum1 + sum2;
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
            // C[i*N + j] += sum1 + sum2;// + sum3 + sum4;
            // s_fin += sum1;// + sum2 + sum3 + sum4;
        }
        // C[i*N] += s_fin;
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

// BL_X must be divisible by 8, for the mul_tr_block function to be correct
#ifndef BL_Y1
    #define BL_Y1_SRC
    #define BL_Y1 32//12//12/*64*/
#endif
#ifndef BL_X
    #define BL_X_SRC
    #define BL_X 64//8 /*64*/
#endif
#ifndef BL_Y2
    #define BL_Y2_SRC
    #define BL_Y2 32//8/*32*/
#endif

#if BL_X == 64
    #define mul_tr_block_BL_X mul_tr_block_64
#else
    #define mul_tr_block_BL_X mul_tr_block
#endif

void mul_tr_4(double * restrict C, register double * restrict A, double * restrict B, int N) {
    size_t i;
    size_t j;
    size_t k;
    const register size_t n = N;

    memset(C, 0, n*n*sizeof(*C));
    
    for (i = 0; i < n; i += BL_Y1) {
        for (j = 0; j < n; j += BL_Y2) {
            for (k = 0; k < n; k += BL_X) {
                // Multiply block (i, j) from A with block (k, j) from B
                /*mul_tr_block_reg(C + i*n + j, A + i*n + k, B + j*n + k, 
                        min(BL_Y1, n - i), 
                        min(BL_X, n - k), 
                        min(BL_Y2, n -j), N);
                        */

                // Multiply block A with block B and store in block C
                // Where block A is BL_Y1 * BL_X
                // Block B is BL_Y2 * BL_X
                // And block C is BL_Y1 * BL_Y2
                // _mm_prefetch((char*)(A + row_size),_MM_HINT_T0);

                /*
                register double aux;
                for (register size_t i1 = 0; i1 < BL_Y1; i1++) {
                    for (register size_t j1 = 0; j1 < BL_Y2; j1++) {
                        aux = C[(i + i1)*n + (j + j1)];
                    }
                }
                aux = C[(i + BL_Y1 - 1)*n + (j + BL_Y2 - 1)] = aux;
                */

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

void mul_tr_5(double * restrict C, register double * restrict A,
        double * restrict B, const register size_t n) {
    register size_t i;
    const register size_t nsq = n*n;

    register double *A_lin = A;
    register double *C_lin = C;

    // For every line in A
    for (i = 0; i < n; i++) {
        register double *B_lin = B;
        register double *B_fin = B + nsq;
        // For every line in B
        while (B_lin < B_fin) {
            register __m128d sum1  = _mm_set_pd(0.0, 0.0);
            register __m128d sum2  = _mm_set_pd(0.0, 0.0);
            register __m128d sum3  = _mm_set_pd(0.0, 0.0);
            // register __m128d sum4  = _mm_set_pd(0.0, 0.0);

            register double *A_fin = A_lin + n;
            while (A_lin < A_fin) {
                __asm__ (   
                            "movapd (%4), %%xmm1\n\t"
                            "movapd (%5), %%xmm2\n\t"
                            "mulpd %%xmm2, %%xmm1\n\t"
                            "addpd %0, %%xmm1\n\t"
                            "movapd %%xmm1, %0\n\t"

                            "movapd 16(%4), %%xmm1\n\t"
                            "movapd 16(%5), %%xmm2\n\t"
                            "mulpd %%xmm2, %%xmm1\n\t"
                            "addpd %1, %%xmm1\n\t"
                            "movapd %%xmm1, %1\n\t"

                            /*
                            "movapd 32(%4), %%xmm1\n\t"
                            "movapd 32(%5), %%xmm2\n\t"
                            "mulpd %%xmm2, %%xmm1\n\t"
                            "addpd %2, %%xmm1\n\t"
                            "movapd %%xmm1, %2\n\t"

                            "movapd 48(%4), %%xmm1\n\t"
                            "movapd 48(%5), %%xmm2\n\t"
                            "mulpd %%xmm2, %%xmm1\n\t"
                            "addpd %3, %%xmm1\n\t"
                            "movapd %%xmm1, %3\n\t"
                            */

                            : "+x"(sum1), "+x"(sum2), "=x"(sum3), "=x"(sum3)
                            : "r"(A_lin), "r"(B_lin)
                            : "%xmm1", "%xmm2");

                A_lin += 4;
                B_lin += 4;
            }
            A_lin -= n;

            __asm__ (   
                        "addpd %1, %0\n\t"
                        //"addpd %3, %2\n\t"
                        "addpd %2, %0\n\t"
                        "pxor %1, %1\n\t"
                        "haddpd %1, %0\n\t"
                        "movlpd %0, (%4)\n\t"
                        "add $8, %4\n\t"

                        : "+x"(sum1), "+x"(sum2), "+x"(sum3),"+x"(sum3), "+rm"(C_lin)
                        :
                        :);
        }
        A_lin += n;
    }
}

void mul_tr_6(double * restrict C, register double * A,
        double * B, const register size_t n) {
    size_t i;
    size_t j;

    register double *A_lin = A;
    register double *C_lin = C;
    register double *B_lin;

    // For every row in A
    for (i = 0; i < n; i++) {
        // For every row in B
        for (j = 0; j < n; j++) {
            B_lin = B + j*n;

            register __m128d sum1  = _mm_set_pd(0.0, 0.0);
            register __m128d sum2  = _mm_set_pd(0.0, 0.0);

            // Compute the dot product bwetwen the row in A and the row in B
            // The row in A starts at A_lin
            // The row in B starts at B_lin
            register double *A_fin = A_lin + n;
            while (A_lin < A_fin) {
                __asm__ volatile (   
                            // sum1[0] += (*A_lin) * (*B_lin);
                            // sum1[1] += (*(A_lin + 1)) * (*(B_lin + 1));
                            "movapd (%2), %%xmm1\n\t"
                            "movapd (%3), %%xmm2\n\t"
                            "mulpd %%xmm2, %%xmm1\n\t"
                            "addpd %0, %%xmm1\n\t"
                            "movapd %%xmm1, %0\n\t"

                            // sum2[0] += (*(A_lin + 2)) * (*(B_lin + 2));
                            // sum2[1] += (*(A_lin + 3)) * (*(B_lin + 3));
                            "movapd 16(%2), %%xmm1\n\t"
                            "movapd 16(%3), %%xmm2\n\t"
                            "mulpd %%xmm2, %%xmm1\n\t"
                            "addpd %1, %%xmm1\n\t"
                            "movapd %%xmm1, %1\n\t"

                            : "+x"(sum1), "+x"(sum2)
                            : "r"(A_lin), "r"(B_lin)
                            : "%xmm1", "%xmm2");

                A_lin += 4;
                B_lin += 4;
            }
            A_lin -= n;

            // Accumulate the sums in C[i][j]
            __asm__ volatile (   
                        // *C_lin = sum1[0] + sum1[1] + sum2[0] + sum2[1];
                        "addpd %1, %0\n\t"
                        "pxor %1, %1\n\t"
                        "haddpd %1, %0\n\t"
                        "movlpd %0, (%2)\n\t"
                        // C_lin++;
                        "add $8, %2\n\t"

                        : "+x"(sum1), "+x"(sum2), "+rm"(C_lin)
                        :
                        :);
        }
        A_lin += n;
    }
}

// Compute A*B^t, where A is a N1 x M block and B is a M x N2 block
// All A, B and C are subblocks of matrices that are N x N
// The output block will be a N1 x N2 matrix and will be added
//  to the block that starts at C
static inline void mul_tr_block_64_7(register double * restrict C, register double * restrict A, 
        register double * restrict B, register size_t N1, register size_t M,
        register size_t N2, register size_t N) {
    register size_t i;
    register size_t j;
    // register size_t k;

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

#undef sum3
#undef sum4

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

static inline void mul_tr_block_7(register double * restrict C, register double * restrict A, 
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

#ifdef BL_Y1_SRC
    #undef BL_Y1
    #define BL_Y1 32//12//12/*64*/
#endif

#ifdef BL_X_SRC
    #undef BL_X
    #define BL_X 64//8 /*64*/
#endif

#ifdef BL_Y2_SRC
    #undef BL_Y2
    #define BL_Y2 30//8/*32*/
#endif

#undef mul_tr_block_BL_X
#if BL_X == 64
    #define mul_tr_block_BL_X mul_tr_block_64_7
#else
    #define mul_tr_block_BL_X mul_tr_block_7
#endif

#undef mul_tr_block
#define mul_tr_block mul_tr_block_7

void mul_tr_7(double * restrict C, register double * restrict A, double * restrict B, int N) {
    size_t i;
    size_t j;
    size_t k;
    const register size_t n = N;

    memset(C, 0, n*n*sizeof(*C));
    
    for (i = 0; i < n; i += BL_Y1) {
        for (j = 0; j < n; j += BL_Y2) {
            for (k = 0; k < n; k += BL_X) {
                // Multiply block (i, j) from A with block (k, j) from B

                // Multiply block A with block B and store in block C
                // Where block A is BL_Y1 * BL_X
                // Block B is BL_Y2 * BL_X
                // And block C is BL_Y1 * BL_Y2


                if (n - k >= BL_X) {
                    mul_tr_block_BL_X(C + i*n + j, A + i*n + k, B + j*n + k, 
                            min(BL_Y1, n - i), 
                            min(BL_X, n - k), 
                            min(BL_Y2, n - j), N);
                } else {
                    mul_tr_block(C + i*n + j, A + i*n + k, B + j*n + k, 
                            min(BL_Y1, n - i), 
                            min(BL_X, n - k), 
                            min(BL_Y2, n - j), N);
                }
                // A_ptr += min(BL_Y1, n - i) * min(BL_X, n - k);
                // B_ptr += min(BL_Y2, n - j) * min(BL_X, n - k);
            }
        }
    }
}

double* my_solver(int N, register double *A, register double* B) {
    printf("MUL SOLVER\n");
    register double *C = malloc(N*N * sizeof(*C));

    // mul_tr(C, A, B, N);
    // mul_tr_2(C, A, B, N);
    // mul_tr_3(C, A, B, N);
    // mul_tr_4(C, A, B, N);
    // mul_tr_5(C, A, B, N);
    // mul_tr_6(C, A, B, N);
    mul_tr_7(C, A, B, N);
    // mul_tr_block_2(C, A, B, N);

    return C;
}
