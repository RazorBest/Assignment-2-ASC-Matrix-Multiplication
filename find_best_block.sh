#!/bin/bash

BL_Y1=16
BL_X=128
BL_Y2=16

for BL_Y1 in {44..52..1}; do
    for BL_X in {12..24..1}; do
        for BL_Y2 in {44..52..1}; do
            rm tema2_mul_tr
            make DEFINES="-DBL_Y1=$BL_Y1 -DBL_X=$BL_X -DBL_Y2=$BL_Y2"
            ./tema2_mul_tr input7 > out
            if [ $? -ne 0 ]; then
                echo "Error"
            fi
            t_elapsed=$(awk 'BEGIN { FS="=" } FNR==2{print $4}' out)

            echo $BL_Y1,$BL_X,$BL_Y2,$t_elapsed >> block_info
        done
    done
done

