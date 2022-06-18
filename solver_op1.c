/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"
#include <string.h>
#include <time.h>

#define max(a, b) (a) > (b) ? a : b;

#define BLOCK_SIZE 40

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

void mul_sym_block_2(double *C, double *A, double *B, int N) {
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

void mul_sym_block_1(double *C, double *A, double *B, int N) {
    register size_t i;
    register size_t j;
    register size_t k;
    register size_t n = N;
    const register size_t block_size = n / 2;

    memset(C, 0, n*n*sizeof(*C));
    register double *A_lin = A;
    for (i = 0; i < n; i++) {
        for (register size_t bl = 0; bl < n; bl += block_size) {
            register double *B_lin = B + bl;
            for (j = 0; j < n; j++) {
                register double sum1 = 0;
                register double sum2 = 0;
                register double sum3 = 0;
                register double sum4 = 0;

                // A_lin will be A + i*n + bl
                // B_lin will be B + j*n + bl
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
                C[i*n + j] += sum1 + sum2 + sum3 + sum4;
            }
            A_lin += block_size;
        }
    }
}

void mul_sym(double *C, double *A, double *B, int N) {
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

void mult_mat_trans_add(double *B, double *A, double *C, size_t n) {
    register size_t i;
    register size_t j;
    const register size_t nsq = n * n;
    register size_t ni = 0;

    register double *A_lin_1 = A;
    register double *A_lin_2 = A + n;
    register double *B_lin = B;
    register double *B_col = B;
    register double *C_lin = C;
    register double *C_col = C;
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
            *B_lin = sum1 + sum2 + *C_lin;
            // The output is symmetrical
            *B_col = sum1 + sum2 + *C_col;

            B_lin++;
            B_col += n;
            C_lin++;
            C_col += n;
        }
        B_lin += i + 2;
        B_col += 2*n + 1 - nsq + ni;
        C_lin += i + 2;
        C_col += 2*n + 1 - nsq + ni;

        A_lin_1 += n;
        A_lin_2 -= nsq - ni;
        ni += n;
        A_lin_2 += 2*n;
    }

    A_lin_1 = A;
    A_lin_2 = A;
    B_lin = B;
    B_col = B;
    C_lin = C;
    C_col = C;
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
        *B_lin = sum1 + sum2 + *C_lin;

        B_lin += n + 1;
        B_col += n + 1;
        C_lin += n + 1;
        C_col += n + 1;

        ni += n;
    }
}

// Compute B = A*A^t
void mult_mat_trans(double *B, double *A, size_t n) {
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
            *B_lin = sum1 + sum2;// + C[i*n + j];//*C_lin;
            // The output is symmetrical
            *B_col = *B_lin;// + C[j*n + i];//*C_col;
            B_lin++;
            B_col += n;
        }
        B_lin += i + 1;
        B_col += n + 1 - nsq + ni;

        A_lin_1 += n;
        A_lin_2 -= nsq - ni;
        ni += n;
        A_lin_2 += n;
    }

    /*
    A_lin_1 = A;
    A_lin_2 = A;
    B_lin = B;
    B_col = B;
    ni = 0;

    for (i = 0; i < n; i++) {
        for (j = i; j < i + 1; j++) {
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
            *B_lin = sum1 + sum2;// + C[i*n + j];// *C_lin;
            // The output is symmetrical
            *B_col = *B_lin;// + C[j*n + i];// *C_col;
            if (B_col != B_lin) {
                printf("Nooo\n");
                return;
            }
            B_lin++;
            B_col += n;
        }
        A_lin_2 += n * (n - i - 1);
        B_lin += n - i - 1;
        B_col += n*(n - i - 1);

        B_lin += i + 1;
        B_col += n + 1 - nsq + ni;

        A_lin_1 += n;
        A_lin_2 -= nsq - ni;
        ni += n;
        A_lin_2 += n;
    }
    */
}

/*
 * Add your optimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
    printf("OPT SOLVER\n");
    double *C = malloc(N*N * sizeof(*C));
    register size_t n = N;
    register size_t i;
    register size_t j;

    double t1, t2;

    // Compute C = A * A^t
    // Knowing A is upper triangular
    t1 = clock();
    for (i = 0; i < n; i++) {
        // register size_t jn = i*n;
        for (j = i; j < n; j++) {
            size_t start = max(i, j);
            C[i*n + j] = 0;
            register double sum = 0;
            register double *A_lin_1 = &A[i*n + start];
            register double *A_lin_2 = &A[j*n + start];
            register double *A_fin = &A[i*n + n];
            while (A_lin_1 < A_fin) {
                sum += *A_lin_1 * (*A_lin_2);
                A_lin_1++;
                A_lin_2++;
            }
            C[i*n + j] = sum;

            // The output is symetrical
            C[j*n + i] = C[i*n + j];
            // jn += n;
        }
    }
    t2 = clock();
    printf("T1 %f\n", (t2 - t1) / CLOCKS_PER_SEC);
    
    return C;
}
