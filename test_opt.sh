#!/bin/bash

./tema2_neopt input5 && mv out3 out4 && ./tema2_opt_m input5 && ./compare out3 out4 0.001
