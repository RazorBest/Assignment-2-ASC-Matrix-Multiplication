CC=gcc
DEFINES=
CFLAGS= -Wall -Werror -O0 -g $(DEFINES)
LIBDIRS=-L/usr/lib64/atlas
LIBS=-l :libsatlas.so.3.10

all: compare tema2_neopt tema2_opt_m tema2_basic tema2_op1 tema2_mul tema2_mul_tr tema2_tr tema2_blas 

tema2_mul_tr: solver_mul_tr.c main.c utils.h
	$(CC) $(CFLAGS) $^ -o $@

tema2_mul: solver_mul.c main.c utils.h
	$(CC) $(CFLAGS) $^ -o $@

tema2_op1: solver_op1.c main.c utils.h
	$(CC) $(CFLAGS) $^ -o $@

tema2_basic: solver_basic.c main.c utils.h
	$(CC) $(CFLAGS) $^ -o $@

tema2_blas: solver_blas.c main.c utils.h
	$(CC) $(CFLAGS) $^ $(LIBDIRS) $(LIBS) -o $@

tema2_neopt: solver_neopt.c main.c utils.h
	$(CC) $(CFLAGS) $^ -o $@

tema2_opt_m: solver_opt.c mat_mul.c mat_mul.h main.c utils.h
	$(CC) $(CFLAGS) $^ -o $@

tema2_tr: solver_tr.c main.c utils.h
	$(CC) $(CFLAGS) $^ -o $@

compare: compare.c utils.h
	$(CC) $(OPT_CFLAGS) $^ -lm -o $@

clean:
	rm -rf compare tema2_blas tema2_neopt tema2_opt_m tema2_basic tema2_op1 tema2_tr
