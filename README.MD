# Programação Paralela

## Primeiro Trabalho

Paralelização do algoritmo de maior subsequência em comum (Longest Common Subsequence) utilizando o OpenMP.

O algoritmo por Yang et al.(2010) foi utilizado como base. Nesse algoritmo a recorrência foi manipulada de tal forma que é possível calcular uma linha inteira em paralelo, ou de outras palavras, a dependência de colunas foi transformada em um cálculo realizado antes do próprio LCS, a matriz P.

## Segundo Trabalho

O segundo trabalho é uma continuação do primeiro, desta vez é necessário paralelizar o LCS utilizando MPI.