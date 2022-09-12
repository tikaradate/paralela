#!bin/bash

# for i in {5000..80000..5000}
# do
#     python3 ../genseq.py $i
# done

# mkdir experimentos
# mkdir inputs
# mv *.in inputs

mpic++ -O3 parallel-mpi.cpp -o paralelo

STARTALL=$(date +%s.%N)
for (( i = 1; i <= 10; i++ ))
do
    mkdir experimentos/exp"${i}"
    echo "-------------- Starting test number ${i} --------------"
    for j in {10000..80000..10000}
    do
        for k in 1 2 4 6
        do  
            echo "Executing parallel with size ${j} and with ${k} threads"
            START=$(date +%s.%N)
            mpirun -np "${k}" paralelo inputs/"${j}"_A.in inputs/"${j}"_B.in 
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
