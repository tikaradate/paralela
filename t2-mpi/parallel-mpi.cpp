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

	char stringBuffer[(size+1)];
	int matrixBuffer[(size+1)*cols];

	int strCounts[nrProcs], strDispls[nrProcs];
	int pCounts[nrProcs], pDispls[nrProcs];

	int strInc = 0, pInc = 0;
	for(int i = 0; i < nrProcs; ++i){
		strCounts[i] = (i < remaining)? size + 1: size;
		strDispls[i] = strInc;
		strInc += strCounts[i];

		pCounts[i] = (i < remaining)? size + 1: size;
		pCounts[i] *= cols;
		pDispls[i] = pInc;
		pInc += pCounts[i];
	}

	MPI_Scatterv(unique.c_str(), strCounts, strDispls, MPI_CHAR, stringBuffer, strCounts[rank], MPI_CHAR, 0, MPI_COMM_WORLD);

	MPI_Scatterv(p, pCounts, pDispls, MPI_INT, matrixBuffer, pCounts[rank], MPI_INT, 0, MPI_COMM_WORLD);

	for(int i = 0; i < strCounts[rank]; ++i){
		for(int j = 1; j < cols; ++j){
			if(bString[j-1] == stringBuffer[i]){
				matrixBuffer[i*(cols) + j] = j;
			} else {
				matrixBuffer[i*(cols) + j] = matrixBuffer[i*(cols) + j-1];
			}
		}
	}
	
	MPI_Gatherv(matrixBuffer, pCounts[rank], MPI_INT, p, pCounts, pDispls, MPI_INT, 0, MPI_COMM_WORLD);
}

int LCS(int **scoreMatrix, string a, std:: string b, int *p, string unique, int rank, int nrProcs) {
	int rows = a.length()+1;
	int cols = b.length()+1;
	
	int size = cols/nrProcs;
	int remaining = (cols % nrProcs);
	
	int scoreBuffer[size];
	int displs[nrProcs], counts[nrProcs];

	MPI_Bcast(p, (unique.length()*cols), MPI_INT, 0, MPI_COMM_WORLD);
	for (int i = 1; i < rows; i++) {
		int c = firstEqual(unique, a[i-1]);

		int inc = 0;
		for(int i = 0; i < nrProcs; ++i){
			counts[i] = (i < remaining)? size + 1 : size;
			displs[i] = inc;
			inc += counts[i] ;
		}
		MPI_Scatterv(scoreMatrix[i], counts, displs, MPI_INT, scoreBuffer, counts[rank], MPI_INT, 0, MPI_COMM_WORLD);

		for (int j = displs[rank]; j < displs[rank] + counts[rank]; j++) {
			int pValue = p[c*cols + j];
			if(pValue == 0){
				scoreBuffer[j - displs[rank]] = scoreMatrix[i-1][j];
			} else {
				scoreBuffer[j - displs[rank]] = max(scoreMatrix[i-1][j], scoreMatrix[i-1][pValue-1] + 1);
			}
		}

		MPI_Allgatherv(scoreBuffer, counts[rank], MPI_INT, scoreMatrix[i], counts, displs, MPI_INT, MPI_COMM_WORLD);
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
