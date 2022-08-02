#!bin/bash

rm *.in

for i in {5000..65000..5000}
do
    python3 genseq.py $i
done

g++ -fopenmp parallel.cpp -O3 -o paralelo
gcc lcs.c -fopenmp -O3 -o sequencial


STARTALL=$(date +%s.%N)
for (( i = 1; i <= 50; i++))
do
    mkdir exp"${i}"
    echo "-------------- Starting test number ${i} --------------"
    for j in {5000..65000..5000}
    do
        echo "Executing sequential with size ${j}"
        START=$(date +%s.%N)
        ./sequencial "${j}"_A.in "${j}"_B.in
        END=$(date +%s.%N)
        DIFF=$(echo "$END - $START" | bc)
        echo $DIFF > exp"${i}"/seq_${k}_${j}.time
        for k in 1 2 3 6 9 12
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
