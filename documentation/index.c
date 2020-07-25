#include <stdio.h>
#include <stdlib.h>

int k_up(int i) {
	return i+M;
};

int k_dn(int i) {
	return ;
};

int main(void) {

	// request M and N
	M = 5;
	N = 5;

	Z = 

	// method 1
	float z = 0.;
	for (int i=-M; i<K; i++) {
		for (int j=k_dn(i); j<=k_up(i); j++) {
			z += Z[i,j];
			printf("%d,%d",i,j);
		}
	}
	
	// method 2

	return 0;
};
