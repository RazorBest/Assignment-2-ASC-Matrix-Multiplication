#!/bin/bash

srun -p nehalem valgrind --tool=cachegrind --branch-sim=yes --cachegrind-out-file=opt_m.cache --log-file=cachegrind.out.opt_m ./tema2_opt_m input_valgrind
srun -p nehalem valgrind --tool=cachegrind --branch-sim=yes --cachegrind-out-file=neopt.cache --log-file=cachegrind.out.neopt ./tema2_neopt input_valgrind

srun -p nehalem valgrind --tool=memcheck --leak-check=full --log-file=opt_m.memory ./tema2_opt_m input_valgrind
srun -p nehalem valgrind --tool=memcheck --leak-check=full --log-file=neopt.memory ./tema2_neopt input_valgrind
