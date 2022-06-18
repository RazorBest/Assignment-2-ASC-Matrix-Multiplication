#ifndef __MAT_MUL_H__
#define __MAT_MUL_H__
#include <stddef.h>

void mul_tr_uptri(double * restrict B, double * restrict A, register size_t n);
void mul_tr(double * restrict C, register double * A, double * B, size_t n);
void mul_tr_self(double * restrict B, double * restrict A, size_t n);

#endif // __MAT_MUL_H__
