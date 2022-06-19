# Assignment 2 - Computer Architectures

## Problem statmement
Perform the following matrix operation:
```
C = B x A x A^t + B^t x B
```
Where:
- A and B are square matrices of size NxN
- A is **upper triangular**
- A^t means the transpose of matrix A
- `x` is the matrix multiplication operation
- `+` is the matrix addition operation

We want to implement the above operation in 3 ways:
- `blas` - a variant that only uses functions from [Blas Atlas](http://www.netlib.org/blas/).
- `neopt` - an unoptimised version.
- `opt_m` - an optimised version which has the same algorithmic complexity as `neopt` (You can't use other algorithms for matrix multiplication; eg. Strassen). The optimisations have to target the source code, specific architecture features, CPU instructions, cache usage, register usage, auxiliary variables. pointer arithmetic etc.

## The unoptimised implementation
This implementation used 3 types of matrix multiplication. But all 3 are
multiplication with a tranposed matrix. We first compute A*A^t, knowing that
A is upper triangular. Then, B*(A*A^t). But, since A*A^t is symmetric, we
can also compute B*(A*A^t)^t. Then, we transpose B and compute B*B^t (the
former (B^t*B). Then, we add B*(A*A^t) with B*B^t.

Running time for N=1200 on nehalem: 22.74 s


## The optimised implementation
This implementation uses the same operations as the unoptimised one.

Running time for N=1200 on nehalem: 1.652 s


## The BLAS implementation
This implementation uses dtrmm, dgemm and dsymm. dtrmm performs A*A^t, knowing
that A is upper triangular. dgemm performs B^t*B. dsymm performs 
B*(A*A^t) + (B^t*B), knowing that (A*A^t) is symmetric.

Running time for N=1200 on nehalem: 1.133 s


## Operation 1: B=A*A^t, where A is upper triangular (unoptimised)

The following operation is computed by mul_tr_up_tri() in the source
code. The unoptimised version uses 3 loops. The first loop iterates
through the rows of A, the second iterates through the columns of A^t and
the third performs the dot product between the row and the column. 

First, we don't need to store A^t, because we can just use A. So, by using this
convention, the columns of A^t become the rows of A. So, the second loop, also
iterates through the lines of A, but with another index. This is good, because
it already makes the unoptimised code cache friendly.

Second, knowing that A is upper triangular, the third loop
is shortened, so that we ingore the elements that are 0.

Lastly, we know that A*A^t is symmeric. So, B[i][j] == B[j][i]. This helps
us shorten the second loop, so instead of starting j from 0, it will start
from i.

## Operation 1: B=A*A^t, where A is upper triangular (optimised)

The following operation is computed by mul_tr_uptri() in the source
code.
Optimisations that were made:
    - use of registers
    - incremental computation of &A[i*n + k] and &A[j*n + k]
    - splitting the sum in 2 and making 2 multiplications at each step
    - using SSE instructions
    - writing the SSE instructions in inline assembly

In the second optimisation, the third loop, instead of A[i*n+k] and
A[j*n+k], two pointers were used. At each step, the pointers were incremented,
instead of incrementing k. By doing this, the number of multiplications for
i*n and j*n was decreased from n^3 to n^2. Also, not using k, seemed to 
be a little faster.

The third optimisation was more experimental. Insteading of writing:
```
sum += A[i*n + k] * A[j*n + k];
k++;
```
I did:
```
sum1 += A[i*j + k] * A[j*n + k];
sum2 += A[i*j + k + 1] * A[j*n + k + 1];
k += 2;
```
This is a loop unrolling. Altough, this also seemed to favor SSE instructions,
(or the pipeline?) for some reason. Still, this has a small problem if the data
is not alligned (because j starts at i, not at 0). So, every second step,
A[j*n] will not be alligned to 16 bytes. Unrolling even more gave worse results

The fourth optimisastion halved the number sequential multiplication.

The fifth optimisation let us control the registers and make sure the stack
wasn't used for storing intermediate values from the addition or
multiplication.

## Operation 2: C=A*B^t (unoptimised)
This operation is implemented in mul_tr(). It's a classic matrix
multiplication algorithm. But, for the second matrix, instead of iterating
through columns, we iterate though lines. This is more cache friendly.


## Operation 2: C=A*B^t (optimised)
This operation is implemented in mul_tr() in the source code.

This was the most expensive operation, so the most optimisations were done
here.

Optimisations that were made:
    - use of registers
    - incremental computation of &A[i*n + k] and &B[j*n + k]
    - got rid of k, and used directly a pointer to A for the stop condition
    - used SSE instructions
    - used inline assembly for SSE instructions
    - used block tiling
    - aggressive loop unrolling for arrays of 64 elements

The second optimisation decreased the number of multiplications that were
made at each step. The third one was found by experimentation.

Multiplying by the tranposed is already cache friendly. But with tiling
and proper choose of block dimensions, we can make it faster. The tiling
used in this implementation is a generic one. We consider a block of
size N1xM in matrix A that is multiplied with a block of size N2xM in
matrix B. The result is stored in an N1xN2 block in matrix C.

The idea is to choose N1, N2, M, such that:
 - the 3 blocks from A, B, C can be kept in the L1 cache.
 - M is divisible by the number of elements in one cache line (8 doubles).
 - N1 and N2 should be close, as we can maximize the number of multiplications
    per pair of blocks (which is N1*N2*M) - keep in mind that this also
    increases the number of accesses to matrix C (N1*N2 reads/writes).

Our cache was 32K, so a reasonable choice would be N1 = N2 = M = 32. This
yielded pretty good results. But it was not better than the unoptimised 
implementation that used the tranposed matrix. Decreasing M and increasing
N1 and N2 was even worse. For some reasons, increasing M was finally giving 
better results. The final choice was:
```
N1 = 32
N2 = 32
M = 64
```
With this configuration, all 3 blocks don't fit in the cache. Still, this
was faster.

One disadvantage of tiling, which is usually omitted, is that it increases the
number of reads/writes to the resulting matrix C. Instead of writing once, we
have to write n/M times, because C[i][j] is calculated from more pairs of
blocks. Since A*B^t is already cache friendly, we have to keep a balance
between the cache misses and the number of reads/writes in matrix C. This might
be an explanation for why the configuration with M=64 was better than the one
with M=32.

The fourth optimisation is the most important: It lets us do 2 multiplications
with one instruction.

The other optimisations are explained in the previous operations.

## Operation 3: A=A^t (unoptimised)
This operation is implemented in transpose(). It uses 2 loops and swaps
A[i][j] with A[j][i]. This is not cache friendly, because every time
we increase j, we access a new line in A.

## Operation 3: A=A^t (optimised)
This operation is implemented in transpose(). 

Optimisations that were made:
    - use of registers
    - incremental computation of &A[i*n + k] and &B[j*n + k]
    - block tiling

It performs the transposition on pairs of square blocks. That's because when
swapping A[i][j] with A[j][i], the line starting at A[j][i] gets loaded in the
cache. So it's a good opportunity to also access A[j][i+1], A[j][i+2] etc. The
block size that was used is 40. Also, the blocks on the main diagonal were
treated separately, because they didn't need to swap with other blocks.

The other optimisations are explained in the previous operations.


## Operation 4: B=A*A^t (unoptimised)
This operation is implemented in mul_tr_self(). It uses 3 loops. Since the
result in symmetric, the second loop can start from j = i, so we only compute
the upper triangular part of B. For the lower part, we know that 
B[i][j] = B[j][i].


## Operation 4: B=A*A^t (optimised)
This operation is implemented in mul_tr_self(). It uses the same ideas and
the same tiling from mul_tr (operation 2). But we only consider the blocks in 
B that contain elements above the main diagonal. The other blocks are
obtained by symmetry of B.

Optimisations that were made: same as for Operation 2 (C=A*B^t).


## Operation 5: C=A+B (unoptimsed)
This operation is implemented in add_mat(). It uses 2 loops and performs
addition.


## Operation 5: C=A+B (optimsed)
This operation is implemented in add_mat().

Optimisations that were made:
    - use of registers
    - loop unrolling (second loop with step 2)

The optimisations are explained in the previous operations.


## Operation 1 optimisation (mul_tr_uptri)
Before:
```
NEOPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_neopt: N=1200: Time=2.285007
NEOPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_neopt: N=1600: Time=5.419080
```
After:
```
OPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_opt_m: N=1200: Time=0.185459
OPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_opt_m: N=1600: Time=0.452464

```
Speedup: 12

## Operation 2 optimisation (mul_tr)
Before:
```
NEOPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_neopt: N=1200: Time=13.543858
NEOPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_neopt: N=1600: Time=32.150166
```
After:
```
OPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_opt_m: N=1200: Time=0.914341
OPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_opt_m: N=1600: Time=2.038876

```
Speedup: 14

## Operation 3 optimisation (transpose)
Before:
```
NEOPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_neopt: N=1200: Time=0.011180
NEOPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_neopt: N=1600: Time=0.022798
```
After:
```
OPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_opt_m: N=1200: Time=0.006017
OPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_opt_m: N=1600: Time=0.013180
```
Speedup: 1.5

## Operation 4 optimisation (mul_tr_self)
Before:
```
NEOPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_neopt: N=1200: Time=6.802978
NEOPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_neopt: N=1600: Time=16.232836
```
After:
```
OPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_opt_m: N=1200: Time=0.484430
OPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_opt_m: N=1600: Time=1.118215
```
Speedup: 14

## Operation 5 optimisation (add_mat)
Before:
```
NEOPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_neopt: N=1200: Time=0.013694
NEOPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_neopt: N=1600: Time=0.023239
```
After:
```
OPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_opt_m: N=1200: Time=0.013448
OPT SOLVER
Run=/export/home/acs/stud/m/marius.pricop/asc/tema2/./tema2_opt_m: N=1600: Time=0.022296
```
Speedup: 1.04



## Cachegrind - cache and branch analysis
After running valgrind with `--tool=cachegrind --branch-sim=yes` for every
implementation, we got the following results:

Unoptimised implementation:
```
I   refs:      4,524,145,664
D   refs:      2,476,487,266  (2,366,154,571 rd   + 110,332,695 wr)
D1  misses:       13,828,027  (   13,578,562 rd   +     249,465 wr)
D1  miss rate:           0.6% (          0.6%     +         0.2%  )

Branches:        110,644,960  (  110,404,280 cond +     240,680 ind)
Mispredicts:         344,238  (      344,004 cond +         234 ind)
Mispred rate:            0.3% (          0.3%     +         0.1%   )
```

Optimised implementation:
```
I   refs:      446,680,656
D   refs:      149,713,933  (134,879,074 rd   + 14,834,859 wr)
D1  misses:      2,901,427  (  2,661,042 rd   +    240,385 wr)
D1  miss rate:         1.9% (        2.0%     +        1.6%  )

Branches:       11,880,660  ( 11,639,947 cond +    240,713 ind)
Mispredicts:       162,776  (    162,525 cond +        251 ind)
Mispred rate:          1.4% (        1.4%     +        0.1%   )
```

BLAS implementation:
```
I   refs:      294,534,896
D   refs:      109,651,060  (103,491,252 rd   + 6,159,808 wr)
D1  misses:      1,921,114  (  1,645,320 rd   +   275,794 wr)
D1  miss rate:         1.8% (        1.6%     +       4.5%  )

Branches:        4,623,281  (  4,366,956 cond +   256,325 ind)
Mispredicts:        73,986  (     73,074 cond +       912 ind)
Mispred rate:          1.6% (        1.7%     +       0.4%   )
```

First, we can see that the unoptimised implementation has more data refs
and instructions. That's because of unnecesary read operations or
multiplications. Eg. - for the naive matrix multiplication, we execute
the following line n^3 times:
```
B[i*n + j] += A[i*n + k] * A[j*n + k];
```
This problem is solved in the optimised version because:
    - We use precalculated values of i*n and j*n.
    - We store the sum in registers, so we don't need to do a read from
        B[i*n + j] at each multiplication.
    - We store in registers other variables that are frequently used.

If we apply the above optimisations to the loops that are on the third
level of imbrication, we gain significant improvment:

Unoptimised version with basic optimisations:
```
I   refs:      1,855,475,025
D   refs:        445,974,445  (442,548,354 rd   + 3,426,091 wr)
D1  misses:       13,828,026  ( 13,578,561 rd   +   249,465 wr)
D1  miss rate:           3.1% (        3.1%     +       7.3%  )

Branches:        110,644,952  (110,404,273 cond +   240,679 ind)
Mispredicts:         344,235  (    344,001 cond +       234 ind)
Mispred rate:            0.3% (        0.3%     +       0.1%   )
```

The next aspect to analyse is the number of branches. The mispredict rate is
low even for the unoptimised version because matrix multiplication isn't very
branch dependent. Most branches come from loops.

For the optimised version, the number of branches is smaller with a factor of
20. I didn't analyse them, but these improvements might come from some loop
unrollings in the optimised version. From what I've tested, the decreased
number of branches wasn't impactful. Because, even if the number of branches
was high, the mispredict rate was always low. The actual benefit of loop
unrollings is that they decreased the number of executed instructions (not the
mispredict rate).

The last aspect, and most important in my opinion, is the L1 cache miss
rate. Suprisingly, the unoptimised version has a 0.6% miss rate. That's
because all the multiplications use a transposed matrix as the second operand.
So, both matrices are iterated row by row. Also, since N <= 1600, a row
of the matrix fits in the cache. So, fixing a row in the first matrix,
and iterating through the rows of the second matrix is good because there's
high chance that the row of the first matrix is only loaded once in the cache.

Still, the number of L1 misses is smaller for the optimised version. That's
because it uses tiling for multiplications and overall has less memory
accesses by using registers.

One thing to observe is that after adding just the basic optimisations,
the number of writes decreases to 3 million. On the other hand, for the
optimised version it is 14 million. That's because of matrix tiling. Instead
of computing C[i][j] at once, it is calculated from multiple pairs of
blocks, so the number of reads from C[i][j] is increased.

We can see that the BLAS implementation is better from every aspect. Therefore,
there can stil be done improvements.

Other aspects of the cache analysis, like L3 misses and instruction misses
weren't that impactful. They were always low, even for the unoptimised version.

## Graphs

![Graph for the 3 implementations](running_time_graph "Running times for different implementations")

I plotted the running times for the 3 implementations of the algorithm, using
the following input:
```
13
40 123 out1
120 123 out1
200 766 out1
400 123 out1
500 234 out1
800 456 out2
1000 346 out2
1200 789 out3
1600 957 out3
1720 957 out3
1800 957 out3
1920 957 out3
2000 957 out3
```
The tests were performed on the nehalem architecture.

In the graph, we can observe that for N=500, there is no significant difference
between the versions. From N >= 800, it is clear that the unoptimised version
is slower with a constant, with running time of 106s for N=2000. For N=2000,
there is not such a big difference between the optimised version and the
BLAS version. The optimised version had T=7.84s and the BLAS version had
T=5.17s

