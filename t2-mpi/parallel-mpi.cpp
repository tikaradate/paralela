#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

using std::vector;
using std::string;

// leitura de arquivos para string
string read_seq(string file) {
	std::ifstream t(file);
	std::stringstream buffer;
	while(buffer << t.rdbuf());
	return buffer.str();
}

// string de letras únicas a partir da string 'str'
void buildUniqueLetters(string str, string &unique){
	for(int i = 0; i < str.length(); ++i){
		// checa se a letra já está no vetor de letras únicas
		if(find(unique.begin(), unique.end(), str[i]) == unique.end())
			unique.push_back(str[i]);
	}
}

// acha o índice da primeira letra igual ao caracter 'c' na string 'str'
int firstEqual(string str, char c){
	for(int i = 0; i < str.length(); ++i){
		if(str[i] == c) return i;
	}
	return -1;
}

// calcula o score da maior subsequencia em comum
// processa a matriz P, a matriz que contém quando determinado caracter de 'bString' foi
void calculatePMatrix(string unique, string bString, int *p, int rank, int nrProcs){
	int rows = unique.length();
	int cols = bString.length()+1;

	int size = rows/nrProcs;
	int remaining = (rows % nrProcs);

	char stringBuffer[size];
	int matrixBuffer[size*cols];

	MPI_Scatter(unique.c_str(), size, MPI_CHAR, stringBuffer, size, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Scatter(p, size*cols, MPI_INT, matrixBuffer, size*cols, MPI_INT, 0, MPI_COMM_WORLD);

	for(int i = 0; i < size; ++i){
		for(int j = 1; j < cols; ++j){
			if(bString[j-1] == stringBuffer[i]){
				matrixBuffer[i*(cols) + j] = j;
			} else {
				matrixBuffer[i*(cols) + j] = matrixBuffer[i*(cols) + j-1];
			}
		}
	}
	
    MPI_Gather(matrixBuffer, size*cols, MPI_INT, p, size*cols, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (rank == 0){
        for (int i = rows - remaining ; i < rows; i++){
            for (int j = 1; j < cols; j++){
                if (bString[j - 1] == unique[i]){
                    p[i*cols + j] = j;
                } else {
                    p[i*cols + j] = p[i*cols + j - 1];
                }
            }
        }
    }
}

int LCS(int **scoreMatrix, string a, std:: string b, int *p, string unique, int rank, int nrProcs) {
	int rows = a.length()+1;
	int cols = b.length()+1;
	
	int size = cols/nrProcs;
	int remaining = (cols % nrProcs);

	int start = size*rank;
	int end = size*rank + size;
	
	int scoreBuffer[size];

	MPI_Bcast(p, (unique.length()*cols), MPI_INT, 0, MPI_COMM_WORLD);
	for (int i = 1; i < rows; i++) {
		int c = firstEqual(unique, a[i-1]);

        MPI_Scatter(scoreMatrix[i], size, MPI_INT, scoreBuffer, size, MPI_INT, 0, MPI_COMM_WORLD);

		for (int j = start; j < end; j++) {
			int pValue = p[c*cols + j];
			if(pValue == 0){
				scoreBuffer[j - start] = scoreMatrix[i-1][j];
			} else {
				scoreBuffer[j - start] = max(scoreMatrix[i-1][j], scoreMatrix[i-1][pValue-1] + 1);
			}
		}

		MPI_Allgather(scoreBuffer, size, MPI_INT, scoreMatrix[i], size, MPI_INT, MPI_COMM_WORLD);

		if (rank == 0){
            for (int j = cols - remaining; j < cols; j++){
                int pValue = p[c*cols + j];
                if(pValue == 0){
                    scoreMatrix[i][j] = scoreMatrix[i-1][j];
                } else {
                    scoreMatrix[i][j] = max(scoreMatrix[i-1][j], scoreMatrix[i-1][pValue-1] + 1);
                }
            }
        }
	}
	return scoreMatrix[a.length()][b.length()];
}

void printMatrix(int n, int m, vector<vector<int>> &mat){
	for(int i = 0; i < n; ++i){
	 	for(int j = 0; j < m; ++j){
	 		std::cout << mat[i][j] << ' ';
		}
	 	std::cout << std::endl;
	}
	std::cout << std::endl;
}

void printPMatrix(int n, int m, int *mat){
	for(int i = 0; i < n; ++i){
	 	for(int j = 0; j < m; ++j){
	 		std::cout << mat[i*m + j] << ' ';
		}
	 	std::cout << std::endl;
	}
	std::cout << std::endl;
}

int main(int argc, char ** argv) {
	if(argc < 3){
		fprintf(stderr, "Necessário arquivos de entrada como argumento\nUso: ./paralelo A.in B.in\n");
		exit(1);
	}	

	MPI_Init(&argc, &argv);

	int nrProcs, rank;

	MPI_Comm_size(MPI_COMM_WORLD, &nrProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	string seqA, seqB;
	
	string unique;

	int seqASize, seqBSize, uniqueSize;
	if (rank == 0){
        seqA = read_seq(argv[1]);
		seqB = read_seq(argv[2]);
		
		seqASize = seqA.size(); 
		seqBSize = seqB.size(); 

		buildUniqueLetters(seqA, unique);
		buildUniqueLetters(seqB, unique);
		uniqueSize = unique.size();
	}

	MPI_Bcast(&seqASize, 1, MPI_INT, 0, MPI_COMM_WORLD);	
	MPI_Bcast(&seqBSize, 1, MPI_INT, 0, MPI_COMM_WORLD);	
	MPI_Bcast(&uniqueSize, 1, MPI_INT, 0, MPI_COMM_WORLD);	

	if(rank != 0){
		seqA.resize(seqASize);
		seqB.resize(seqBSize);
		unique.resize(uniqueSize);
	}

	MPI_Bcast(const_cast<char*>(seqA.data()), seqASize, MPI_CHAR, 0, MPI_COMM_WORLD);	
	MPI_Bcast(const_cast<char*>(seqB.data()), seqBSize, MPI_CHAR, 0, MPI_COMM_WORLD);	
	MPI_Bcast(const_cast<char*>(unique.data()), uniqueSize, MPI_CHAR, 0, MPI_COMM_WORLD);	

	int **scoreMatrix = (int **) calloc((seqASize+1), sizeof(int*));
	for(int i = 0; i < (seqASize+1); ++i)
        scoreMatrix[i] = (int *) calloc((seqBSize+1), sizeof(int));
	
	int *currLine = (int *) calloc((seqBSize+1), sizeof(int));
	int *prevLine = (int *) calloc((seqBSize+1), sizeof(int));

	int *p = (int *) calloc((uniqueSize)*(seqBSize+1), sizeof(int));

	calculatePMatrix(unique, seqB, p, rank, nrProcs);

	int score = LCS(scoreMatrix, seqA, seqB, p, unique, rank, nrProcs);
	
	MPI_Finalize();

	if(rank == 0)
		std::cout << "\nScore: " << score << std::endl;

	return EXIT_SUCCESS;
}
