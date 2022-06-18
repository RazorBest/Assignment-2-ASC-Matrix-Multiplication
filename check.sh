#!/bin/bash

if [ $1 == 'neopt' ]; then
    ./tema2_neopt input1 &&
    ./compare out1 out1.good 0.0001 &&
    ./compare out2 out2.good 0.0001 &&
    ./compare out3 out3.good 0.0001
fi

if [ $1 == "opt_m" ]; then
    ./tema2_opt_m input1 &&
    ./compare out1 out1.good 0.0001 &&
    ./compare out2 out2.good 0.0001 &&
    ./compare out3 out3.good 0.0001
fi

if [ $1 == "blas" ]; then
    ./tema2_blas input1 &&
    ./compare out1 out1.good 0.0001 &&
    ./compare out2 out2.good 0.0001 &&
    ./compare out3 out3.good 0.0001
fi
