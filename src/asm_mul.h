#pragma once

// Assembly routine for vector-vector multiplication of size 8
// sum1 and sum2 must be xmm* registers
#define ASM_MUL_ADD_8(sum1, sum2, A, B) __asm__ (   \
        /* sum1 += A[0:1]*B[0:1]*/ \
        "movapd (%2), %%xmm1\n\t"\
        "movapd (%3), %%xmm2\n\t"\
        "mulpd %%xmm2, %%xmm1\n\t"\
        "addpd %0, %%xmm1\n\t"\
        "movapd %%xmm1, %0\n\t"\
        \
        /* sum2 += A[2:3]*B[2:3]*/ \
        "movapd 16(%2), %%xmm1\n\t"\
        "movapd 16(%3), %%xmm2\n\t"\
        "mulpd %%xmm2, %%xmm1\n\t"\
        "addpd %1, %%xmm1\n\t"\
        "movapd %%xmm1, %1\n\t"\
        \
        /* sum1 += A[4:5]*B[4:5]*/ \
        "movapd 32(%2), %%xmm1\n\t"\
        "movapd 32(%3), %%xmm2\n\t"\
        "mulpd %%xmm2, %%xmm1\n\t"\
        "addpd %0, %%xmm1\n\t"\
        "movapd %%xmm1, %0\n\t"\
        \
        /* sum2 += A[6:7]*B[6:7]*/ \
        "movapd 48(%2), %%xmm1\n\t"\
        "movapd 48(%3), %%xmm2\n\t"\
        "mulpd %%xmm2, %%xmm1\n\t"\
        "addpd %1, %%xmm1\n\t"\
        "movapd %%xmm1, %1\n\t"\
        \
        : "+x"(sum1), "+x"(sum2)\
        : "rm"(A), "rm"(B)\
        : "%xmm1", "%xmm2")

// Assembly routine for *A = sum1[0] + sum1[1] + sum2[0] + sum2[1]
// sum1 and sum2 must be xmm* registers
#define ASM_ACCUMULATE_S1_S2(sum1, sum2, A) __asm__ (   \
                        "addpd %1, %0\n\t"\
                        "pxor %1, %1\n\t"\
                        "haddpd %1, %0\n\t"\
                        "movlpd (%2), %1\n\t"\
                        "addpd %1, %0\n\t"\
                        "movlpd %0, (%2)\n\t"\
                        \
                        : "+x"(sum1), "+x"(sum2), "+rm"(A)\
                        :\
                        :)
