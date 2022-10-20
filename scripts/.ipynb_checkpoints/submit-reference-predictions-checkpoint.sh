#!/usr/bin/bash

-n 16
-t 8:00:00
-A covid-ct

qsub-gpu -A covid-ct -n 16 -t 8:00:00 --env MYVAR=value1 -i inputdata -O Project1_out