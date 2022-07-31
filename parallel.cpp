#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

typedef unsigned short mtype;

/* Read sequence from a file to a char vector.
 Filename is passed as parameter */

std::string read_seq(std::string file) {
	std::ifstream t(file);
	std::stringstream buffer;
	while(buffer << t.rdbuf());
	//return sequence pointer
	return buffer.str();
}

void buildUniqueLetters(std::string string, std::string &unique){
	for(int i = 0; i < string.length(); ++i){
		// checa se a letra já está no vetor de letras únicas
		if(find(unique.begin(), unique.end(), string[i]) == unique.end())
			unique.push_back(string[i]);
	}
}

int firstEqual(std::string str, char c){
	for(int i = 0; i < str.length(); ++i){
		if(str[i] == c) return i;
	}
}

void pMatrix(std::string unique, std::string bString, std::vector<std::vector<mtype>> &p){
	#pragma omp parallel for
	for(int i = 0; i < unique.length(); ++i){
		for(int j = 1; j < bString.length()+1; ++j){
			if(bString[j-1] == unique[i]){
				p[i][j] = j;
			} else {
				p[i][j] = p[i][j-1];
			}
		}
	}
}

mtype LCS(std::vector<std::vector<mtype>> &scoreMatrix, std::string a, std:: string b, std::vector<std::vector<mtype>> p, std::string unique) {
	for (int i = 1; i < a.length(); i++) {
		mtype c = firstEqual(unique, a[i-1]);
		#pragma omp parallel for schedule(static)
		for (int j = 0; j < b.length(); j++) {
			mtype m;
			if(p[c][j] == 0){
				m = max(scoreMatrix[i-1][j], 0);
			} else {
				m = max(scoreMatrix[i-1][j], scoreMatrix[i-1][p[c][j]-1] + 1);
			}
			scoreMatrix[i][j] = m;
		}
	}
	return scoreMatrix[a.length()-1][b.length()-1];
}

int main(mtype argc, char ** argv) {
	// sequence pointers for both sequences
	std::string seqA, seqB;
	
	std::string unique = "AGTC";

	//read both sequences
	seqA = read_seq("fileA.in");
	seqB = read_seq("fileB.in");

	// std::cout << seqA << ' ' << seqB << std::endl;

	// buildUniqueLetters(seqA, unique);
	// buildUniqueLetters(seqB, unique);

	// for(mtype i = 0; i < unique.length(); ++i){
	// 	std::cout << unique[i] << ' ';
	// }
	std::cout << unique << std::endl;

	// allocate LCS score matrix
	std::vector<std::vector<mtype>> scoreMatrix(seqA.length(), std::vector<mtype>(seqB.length()));
	// std::cout << "allocated score matrix" << std::endl;

	std::vector<std::vector<mtype>> p (unique.length(), std::vector<mtype>(seqB.length()+1));
	// std::cout << "allocated p matrix" << std::endl;

	pMatrix(unique, seqB, p);
	// std::cout << "finished processing p matrix " << std::endl;
	//fill up the rest of the matrix and return final score (element locate at the last line and collumn)
	mtype score = LCS(scoreMatrix, seqA, seqB, p, unique);

	/* if you wish to see the entire score matrix,
	 for debug purposes, define DEBUGMATRIX. */
	
	// for(mtype i = 0; i < seqA.length(); ++i){
	// 	for(mtype j = 0; j < seqB.length(); ++j){
	// 		std::cout << scoreMatrix[i][j] << ' ';
	// 	}
	// 	std::cout << std::endl;
	// }
	// std::cout << std::endl;


	// for(mtype i = 0; i < unique.length(); ++i){
	// 	for(mtype j = 0; j < seqB.length(); ++j){
	// 		std::cout << p[i][j] << ' ';
	// 	}
	// 	std::cout << std::endl;
	// }

	//prmtype score
	printf("\nScore: %d\n", score);

	//free score matrix
	// freeScoreMatrix(scoreMatrix, sizeB);

	return EXIT_SUCCESS;
}
