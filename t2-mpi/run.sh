#!bin/bash

for i in {5000..80000..5000}
do
    python3 ../genseq.py $i
done

mkdir experimentos
mkdir inputs
mv *.in inputs

mpic++ parallel-mpi.cpp -O3 -o paralelo

STARTALL=$(date +%s.%N)
for (( i = 1; i <= 20; i++ ))
do
    mkdir experimentos/exp"${i}"
    echo "-------------- Starting test number ${i} --------------"
    for j in {5000..80000..5000}
    do
        for k in 1 2 4 6 8 12
        do  
            echo "Executing parallel with size ${j} and with ${k} threads"
            START=$(date +%s.%N)
            mpirun paralelo inputs/"${j}"_A.in inputs/"${j}"_B.in -np "${k}"
            END=$(date +%s.%N)
            DIFF=$(echo "$END - $START" | bc)
            echo $DIFF > experimentos/exp"${i}"/par_${k}_${j}.time
        done
    done
done
ENDALL=$(date +%s.%N)
DIFF=$(echo "$ENDALL - $STARTALL" | bc)
echo "FIN~"
echo $DIFF
