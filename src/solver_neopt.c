/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

double* my_solver(int N, double *A, double* B) {
	printf("NEOPT SOLVER\n");

	double *At = calloc(N * N, sizeof(double)); 
	double *Bt = calloc(N * N, sizeof(double));

	double *AxB = calloc(N * N, sizeof(double));
	double *AxBxAt = calloc(N * N, sizeof(double));
	double *BtxBt = calloc(N * N, sizeof(double));

	double *C = calloc(N * N, sizeof(double));

	int i, j, k;

	// A transpose
	for (i = 0; i < N; i++) {
		for (j = i; j < N; j++) {
			At[j * N + i] = A[i * N + j];
		}
	}

	// B transpose
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			Bt[j * N + i] = B[i * N + j];
		}
	}

	double sum = 0;

	// A * B
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			sum = 0;
			for (k = i; k < N; k++) {
				sum += A[i * N + k] * B[k * N + j];
			}
			AxB[i * N + j] = sum;
		}
	}

	// A * B * At
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			sum = 0;
			for (k = j; k < N; k++) {
				sum += AxB[i * N + k] * At[k * N + j];
			}
			AxBxAt[i * N + j] = sum;
		}
	}

	// Bt * Bt
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			sum = 0;
			for (k = 0; k < N; k++) {
				sum += Bt[i * N + k] * Bt[k * N + j];
			}
			BtxBt[i * N + j] = sum;
		}
	}


	// C = A * B * At + Bt * Bt
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			C[i * N + j] = AxBxAt[i * N + j] + BtxBt[i * N + j];
		}
	}

	free(At);
	free(Bt);
	free(AxB);
	free(AxBxAt);
	free(BtxBt);

	return C;
}
