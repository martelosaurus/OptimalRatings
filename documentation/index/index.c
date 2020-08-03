#include <stdio.h>
#define M 5
#define N 5

int k_up(int i) {
	return 0;
};

int k_dn(int i) {
	return 0;
};

int print_matrix(float **A,int m,int n) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			printf("%f ",A[i,j]);
		printf("\n");
		return 0;
};

int main(void) {

	// request M and N
	//int M, N;
	//printf("M:");
	//scanf("%d",&M);
	//printf("N:");
	//scanf("%d",&N);
	
	float A[M][N];
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			A[i,j] = rand(); 

	// number of interior breakpoints
	int K = M*N;
	printf("%d",K);

	print_matrix(A,M,N);

	// allocate space for matrices
	//float **A = malloc((M+K)*sizeof(**char));
	//float **B = malloc(sizeof(**char));
	//float **C = malloc(sizeof(**char));

	/*
	// method 1
	float z = 0.;
	for (int i=-M; i<K; i++) {
		for (int j=k_dn(i); j<=k_up(i); j++) {
			z += Z[i,j];
			printf("%d,%d",i,j);
		}
	}
	
	*/
	// method 2
	
	// print matrices

	return 0;
};
