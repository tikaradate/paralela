#!bin/bash

g++ -fopenmp parallel.cpp -O3 -o paralelo

STARTALL=$(date +%s.%N)
for (( i = 1; i <= 50; i++))
do
    echo "-------------- Starting test number ${i} --------------"
    for j in {5000..65000..5000}
    do
        echo $DIFF > exp"${i}"/seq_${k}_${j}.time
        for k in 4 5
        do  
            export OMP_NUM_THREADS=${k}
            echo "Executing parallel with size ${j} and with ${k} threads"
            START=$(date +%s.%N)
            ./paralelo "${j}"_A.in "${j}"_B.in
            END=$(date +%s.%N)
            DIFF=$(echo "$END - $START" | bc)
            echo $DIFF > exp"${i}"/par_${k}_${j}.time
        done
    done
done
ENDALL=$(date +%s.%N)
DIFF=$(echo "$ENDALL - $STARTALL" | bc)
echo "FIN"
echo $DIFF
