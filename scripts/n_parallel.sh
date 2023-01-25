#!/bin/bash

N=10
(
for thing in `seq 1 22` Y; do 
   ((i=i%N)); ((i++==0)) && wait
   echo "$i" & 
done
)