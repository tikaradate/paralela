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

typedef unsigned short mtype;

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
void calculatePMatrix(string unique, string bString, mtype *p, int rank, int nrProcs){
	int rows = unique.length();
	int cols = bString.length()+1;

	int size = rows/nrProcs;
	int remaining = rows % nrProcs;

	// char stringBuffer[size];
	mtype matrixBuffer[size*cols];

	// MPI_Scatter(unique.c_str(), size, MPI_CHAR, &stringBuffer, size, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Scatter(p, size*cols, MPI_SHORT, &matrixBuffer, size*cols, MPI_SHORT, 0, MPI_COMM_WORLD);

	for(int i = 0; i < size; ++i){
		for(int j = 1; j < cols; ++j){
			if(bString[j-1] == unique[i]){
				matrixBuffer[i*(cols) + j] = j;
			} else {
				matrixBuffer[i*(cols) + j] = matrixBuffer[i*(cols) + j-1];
			}
		}
	}
	
    MPI_Gather(matrixBuffer, size*cols, MPI_SHORT, p, size*cols, MPI_SHORT, 0, MPI_COMM_WORLD);
	
	if (rank == 0){        
        for (int i = remaining ; i < rows; i++){
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

mtype LCS(vector<vector<mtype>> &scoreMatrix, string a, std:: string b, mtype *p, string unique) {
	for (int i = 1; i < a.length()+1; i++) {
		mtype c = firstEqual(unique, a[i-1]);
		// #pragma omp parallel for
		for (int j = 0; j < b.length()+1; j++) {
			// aplica a programação dinâmica utilizando como base a recursão e 
			// também utilizando a matriz P
			if(i == 0 || j == 0){
				scoreMatrix[i][j] = 0;
			} else if(p[c*(b.length()+1) + j] == 0){
				scoreMatrix[i][j] = max(scoreMatrix[i-1][j], 0);
			} else {
				scoreMatrix[i][j] = max(scoreMatrix[i-1][j], scoreMatrix[i-1][p[c*(b.length()+1) + j]-1] + 1);
			}
		}
	}
	return scoreMatrix[a.length()][b.length()];
}

void printMatrix(int n, int m, vector<vector<mtype>> &mat){
	for(mtype i = 0; i < n; ++i){
	 	for(mtype j = 0; j < m; ++j){
	 		std::cout << mat[i][j] << ' ';
		}
	 	std::cout << std::endl;
	}
	std::cout << std::endl;
}

void printPMatrix(int n, int m, mtype *mat){
	for(mtype i = 0; i < n; ++i){
	 	for(mtype j = 0; j < m; ++j){
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

	// std::cout << rank << std::endl;
	string seqA, seqB;
	
	string unique;

	int seqSize, uniqueSize;
	if (rank == 0){
        seqA = read_seq(argv[1]);
		seqB = read_seq(argv[2]);
		// seqA = "stringA";
		// seqB = "stringB";

		seqSize = seqA.size(); 

		buildUniqueLetters(seqA, unique);
		buildUniqueLetters(seqB, unique);
		uniqueSize = unique.size();
	}

	MPI_Bcast(&seqSize, 1, MPI_INT, 0, MPI_COMM_WORLD);	
	MPI_Bcast(&uniqueSize, 1, MPI_INT, 0, MPI_COMM_WORLD);	

	if(rank != 0){
		seqA.resize(seqSize);
		seqB.resize(seqSize);
		unique.resize(uniqueSize);
	}

	MPI_Bcast(const_cast<char*>(seqA.data()), seqSize, MPI_CHAR, 0, MPI_COMM_WORLD);	
	MPI_Bcast(const_cast<char*>(seqB.data()), seqSize, MPI_CHAR, 0, MPI_COMM_WORLD);	
	MPI_Bcast(const_cast<char*>(unique.data()), uniqueSize, MPI_CHAR, 0, MPI_COMM_WORLD);	

	vector<vector<mtype>> scoreMatrix (seqA.length()+1, vector<mtype>(seqB.length()+1));
	// mtype *scoreMatrix = (mtype *) malloc((seqSize+1)*(seqSize+1)*sizeof(mtype));

	// vector<vector<mtype>> p (unique.length(), vector<mtype>(seqB.length()+1));
	mtype *p = (mtype *) malloc((uniqueSize)*(seqSize+1)*sizeof(mtype));

	calculatePMatrix(unique, seqB, p, rank, nrProcs);

	// mtype score = LCS(scoreMatrix, seqA, seqB, p, unique);
	mtype score = LCS(scoreMatrix, seqA, seqB, p, unique);
	
	MPI_Finalize();

	//printMatrix(seqA.length(), seqB.length(), scoreMatrix);

	std::cout << "\nScore: " << score << std::endl;

	return EXIT_SUCCESS;
}
