/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"
#include <string.h>

typedef struct {
    double *mat;
    size_t n, m;
} matrix;

static matrix* new_matrix(size_t n, size_t m) {
    matrix *mat = malloc(sizeof(*mat));
    if (mat == NULL) {
        return NULL;
    }

    mat->mat = calloc(n, m * sizeof(mat->mat));
    if (mat->mat == NULL) {
        free(mat);
        return NULL;
    }

    mat->n = n;
    mat->m = m;

    return mat;
}

static void free_matrix(matrix *mat) {
    free(mat->mat);
    free(mat);
}

static matrix* pack_matrix(double *a, size_t n, size_t m) {
    matrix *mat = malloc(sizeof(*mat));
    if (mat == NULL) {
        return NULL;
    }

    mat->mat = a;
    mat->n = n;
    mat->m = m;

    return mat;
}

static double* unpack_matrix(matrix *mat) {
    double *a = mat->mat;
    free(mat);

    return a;
}

// Perform matrix multiplication m_c = m_a x m_b
static int matrix_mul(matrix m_c, matrix m_a, matrix m_b) {
    size_t i, j, k;
    size_t n1 = m_a.n, m1 = m_a.m;
    size_t n2 = m_b.n, m2 = m_b.m;
    double *a = m_a.mat;
    double *b = m_b.mat;
    double *c = m_c.mat;

    if (m1 != n2) {
        return 0;
    }

    for (i = 0; i < n1; i++) {
        for (j = 0; j < m2; j++) {
            c[i*m2 + j] = 0;
            for (k = 0; k < m1; k++) {
                c[i*m2 + j] += a[i*m1 + k] * b[k*m2 + j];
            }
        }
    }

    return 1;
}

// Perform matrix addition m_c = m_a x m_b
static int matrix_add(matrix m_c, matrix m_a, matrix m_b) {
    size_t i, j;
    size_t n1 = m_a.n, m1 = m_a.m;
    size_t n2 = m_b.n, m2 = m_b.m;
    double *a = m_a.mat;
    double *b = m_b.mat;
    double *c = m_c.mat;

    if (n1 != n2 || m1 != m2) {
        return 0;
    }

    for (i = 0; i < n1; i++) {
        for (j = 0; j < m1; j++) {
            c[i*m1 + j] = a[i*m1 + j] + b[i*m1 + j];
        }
    }

    return 1;
}

// Perform matrix tranposition m_a = m_a^t
static int matrix_sq_transpose_inplace(matrix m_a) {
    size_t i, j;
    size_t n = m_a.n, m = m_a.m;
    double *a = m_a.mat;

    if (n != m) {
        return 0;
    }

    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            // Swap a[i][j] with a[j][i]
            double aux = a[i*m + j]; 
            a[i*m + j] = a[j*m + i];
            a[j*m + i] = aux;
        }
    }

    return 1;
}

// Perform copy m_a = m_b
static int matrix_copy(matrix m_a, matrix m_b) {
    size_t n1 = m_a.n, m1 = m_a.m;
    size_t n2 = m_b.n, m2 = m_b.m;
    double *a = m_a.mat;
    double *b = m_b.mat;

    if (n1 != n2 || m1 != m2) {
        return 0;
    }

    memcpy(a, b, sizeof(*a) * n1 * m1);

    return 1;
}

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	printf("BASIC SOLVER\n");

    matrix *m_a, *m_b, *m_c, *m_d;
    m_a = pack_matrix(A, N, N);
    m_b = pack_matrix(B, N, N);
    m_c = new_matrix(N, N);
    m_d = new_matrix(N, N);

    matrix_mul(*m_c, *m_b, *m_a);
    matrix_sq_transpose_inplace(*m_a);
    matrix_mul(*m_d, *m_c, *m_a);

    matrix_copy(*m_a, *m_b);
    matrix_sq_transpose_inplace(*m_b);
    matrix_mul(*m_c, *m_b, *m_a);

    matrix_add(*m_a, *m_c, *m_d);
    matrix_copy(*m_c, *m_a);

    free_matrix(m_d);
    unpack_matrix(m_b);
    unpack_matrix(m_a);

	return unpack_matrix(m_c);
}
