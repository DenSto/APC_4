#!/bin/bash
if [ $# -ne 1 ]; then
    echo "Usage: join_dat.sh nprocs"
    exit
fi

nend=$(expr  $1 - 1)
cat out_0.dat > combined.dat
for i in `seq 1 $nend`; do
    cat out_$i.dat >> combined.dat
done
