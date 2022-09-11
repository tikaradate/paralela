#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#define MAX_CHAR 256

typedef unsigned short mtype;

// leitura de arquivos para string
char* read_seq(char *fname) {
	//file pointer
	FILE *fseq = NULL;
	//sequence size
	long size = 0;
	//sequence pointer
	char *seq = NULL;
	//sequence index
	int i = 0;

	//open file
	fseq = fopen(fname, "rt");
	if (fseq == NULL ) {
		printf("Error reading file %s\n", fname);
		exit(1);
	}

	//find out sequence size to allocate memory afterwards
	fseek(fseq, 0, SEEK_END);
	size = ftell(fseq);
	rewind(fseq);

	//allocate memory (sequence)
	seq = (char *) calloc(size + 1, sizeof(char));
	if (seq == NULL ) {
		printf("Erro allocating memory for sequence %s.\n", fname);
		exit(1);
	}

	//read sequence from file
	while (!feof(fseq)) {
		seq[i] = fgetc(fseq);
		if ((seq[i] != '\n') && (seq[i] != EOF))
			i++;
	}
	//insert string terminator
	seq[i] = '\0';

	//close file
	fclose(fseq);

	//return sequence pointer
	return seq;
}


// string de letras únicas a partir da string 'str'
char *buildUniqueLetters(char *seqA, char *seqB){
	char *seq = malloc(sizeof("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")+1);
	seq = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	return seq;
}

// acha o índice da primeira letra igual ao caracter 'c' na string 'str'
int firstEqual(char *str, char c){
	for(int i = 0; i < strlen(str); ++i){
		if(str[i] == c) return i;
	}
	return -1;
}

// calcula o score da maior subsequencia em comum
// processa a matriz P, a matriz que contém quando determinado caracter de 'bString' foi
void calculatePMatrix(char *unique, char *bString, mtype *p, int rank, int nrProcs){
	int rows = strlen(unique);
	int cols = strlen(bString)+1;
	int size = rows/nrProcs;
	int remaining = rows % nrProcs;

	char stringBuffer[size];
	mtype matrixBuffer[size*cols];

	// MPI_Scatter(unique.c_str(), size, MPI_CHAR, &stringBuffer, size, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Scatter(p, size*cols, MPI_SHORT, matrixBuffer, size*cols, MPI_SHORT, 0, MPI_COMM_WORLD);

	for(int i = 0; i < size; ++i){
		for(int j = 1; j < cols; ++j){
			if(bString[j-1] == unique[i]){
				matrixBuffer[i*(cols) + j] = j;
			} else {
				matrixBuffer[i*(cols) + j] = matrixBuffer[i*(cols) + j-1];
			}
		}
	}
	// MPI_Barrier(MPI_COMM_WORLD);
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

mtype LCS(mtype **scoreMatrix, char *a, char *b, mtype *p, char *unique, int rank, int nrProcs) {
	int rows = strlen(a)+1;
	int cols = strlen(b)+1;
	
	int size = cols/nrProcs;
	int remaining = cols % nrProcs;
	int start = size*rank;
	if(start == 0) start = 1;
	int end = size*rank + size;
	int scoreBuffer[size];

	MPI_Bcast(p, (strlen(unique)*cols), MPI_SHORT, 0, MPI_COMM_WORLD);
	for (int i = 1; i < rows; i++) {
		mtype c = firstEqual(unique, a[i-1]);

        MPI_Scatter(scoreMatrix[i], size, MPI_SHORT, scoreBuffer, size, MPI_SHORT, 0, MPI_COMM_WORLD);

		for (int j = start; j < end; j++) {
			int pValue = p[c*cols + j];
			if(pValue == 0){
				// scoreMatrix[i][j] = max(scoreMatrix[i-1][j], 0);
				scoreBuffer[j - start] = max(scoreMatrix[i-1][j], 0);
			} else {
				// scoreMatrix[i][j] = max(scoreMatrix[i-1][j], scoreMatrix[i-1][pValue-1] + 1);
				scoreBuffer[j - start] = max(scoreMatrix[i-1][j], scoreMatrix[i-1][pValue-1] + 1);
			}
		}
		
		MPI_Allgather(scoreBuffer, size, MPI_SHORT, scoreMatrix[i], size, MPI_SHORT, MPI_COMM_WORLD);
	
		if (rank == 0){
            for (int j = remaining; j < cols; j++){
                int pValue = p[c*cols + j];
                if(pValue){
                    scoreMatrix[i][j] = max(scoreMatrix[i-1][j], scoreMatrix[i-1][pValue - 1] + 1);
                } else {
                    scoreMatrix[i][j] = scoreMatrix[i-1][j];
                }
            }
        }
	}
	return scoreMatrix[strlen(a)][strlen(b)];
}

void printMatrix(int n, int m, mtype **mat){
	for(mtype i = 0; i < n; ++i){
		for(mtype j = 0; j < m; ++j){
			printf("%d ", mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}
void printPMatrix(int n, int m, mtype *mat){
	for(mtype i = 0; i < n; ++i){
		for(mtype j = 0; j < m; ++j){
			printf("%d ", mat[i*m + j]);
		}
		printf("\n");
	}
	printf("\n");
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

	char *seqA, *seqB, *unique;

	int seqASize, seqBSize, uniqueSize;
	
	if(rank == 0){
		seqA = read_seq(argv[1]);
		seqB = read_seq(argv[2]);

		unique = buildUniqueLetters(seqA, seqB);

		seqASize = strlen(seqA)+1;
		seqASize = strlen(seqB)+1;
		uniqueSize = strlen(unique)+1;
	}

	MPI_Bcast(&seqASize, 1, MPI_INT, 0, MPI_COMM_WORLD);	
	MPI_Bcast(&seqBSize, 1, MPI_INT, 0, MPI_COMM_WORLD);	
	MPI_Bcast(&uniqueSize, 1, MPI_INT, 0, MPI_COMM_WORLD);	

	if(rank != 0){
		seqA = malloc(seqASize);
		seqB = malloc(seqBSize);
		unique = malloc(uniqueSize);
	}

	MPI_Bcast(seqA, seqASize, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(seqB, seqBSize, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(unique, uniqueSize, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	mtype **scoreMatrix = (mtype **) calloc(seqASize, sizeof(mtype *));
	for(int i = 0; i < seqASize; ++i)
		scoreMatrix[i] = (mtype *) calloc(seqBSize, sizeof(mtype));

	mtype *p = (mtype *) calloc((uniqueSize)*(seqBSize), sizeof(mtype));

	calculatePMatrix(unique, seqB, p, rank, nrProcs);

	mtype score = LCS(scoreMatrix, seqA, seqB, p, unique, rank, nrProcs);
	
	if(rank == 0)
		printf("\nScore: %d\n", score);

	MPI_Finalize();
	return EXIT_SUCCESS;
}