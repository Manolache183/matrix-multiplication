/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

double* my_solver(int N, double *A, double* B) {
	printf("OPT SOLVER\n");

	double *At = calloc(N * N, sizeof(double)); 
	double *Bt = calloc(N * N, sizeof(double));

	double *AxB = calloc(N * N, sizeof(double));
	double *AxBxAt = calloc(N * N, sizeof(double));
	double *BtxBt = calloc(N * N, sizeof(double));
	
	double *C = calloc(N * N, sizeof(double));

	int i, j, k;
	int blockSize = 40;

	// A transpose
	register double *p = A;
	for (i = 0; i < N; i++) {
		p += i; // since j starts from i (A is upper triangular)
		register double *pt = At + i * N + i; // since j starts from i (A is upper triangular)
		
		for (j = i; j < N; j++) {
			*pt = *p;
			
			p++;
			pt += N;
		}
	}

	// B transpose
	p = B;
	for (i = 0; i < N; i++) {
		register double *pt = Bt + i;
		
		for (j = 0; j < N; j++) {
			*pt = *p;
			
			p++;
			pt += N;
		}
	}

	//A * B, i k j
	int bi, bk, bj;
	for(bi = 0; bi < N; bi += blockSize) {
        for(bk = bi; bk < N; bk += blockSize) { // A upper triang
            for(bj = 0; bj < N; bj += blockSize) {
                for(i = 0; i < blockSize; i++) {
					register int is = (i + bi) * N;
					
					register double *pt = A + is + bk;

                    for(k = 0; k < blockSize; k++) {
						register int ks = (k + bk) * N;

						register double *p = AxB + is + bj;
						register double *ptt = B + ks + bj;

                        for(j = 0; j < blockSize; j++) {
                           *p += (*pt) * (*ptt);
						   p++;
						   ptt++;
						}

						pt++;
					}

					pt += blockSize;
					p += blockSize;
				}
	}}}

	// A * B * At
	p = AxBxAt;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			register double *pt = AxB + i * N + j; // since k starts from j (A is upper triangular)
			register double *ptt = At + j * N + j; // since k starts from j (A is upper triangular)
			for (k = j; k < N; k++) {
				*p += (*pt) * (*ptt);
				pt++;
				ptt += N;
			}

			p++;
		}
	}

	// Bt * Bt
	for(bi = 0; bi < N; bi += blockSize) {
        for(bk = 0; bk < N; bk += blockSize) {
            for(bj = 0; bj < N; bj += blockSize) {
                for(i = 0; i < blockSize; i++) {
					register int is = (i + bi) * N;
					
					register double *pt = Bt + is + bk;

                    for(k = 0; k < blockSize; k++) {
						register int ks = (k + bk) * N;

						register double *p = BtxBt + is + bj;
						register double *ptt = Bt + ks + bj;

                        for(j = 0; j < blockSize; j++) {
                           *p += (*pt) * (*ptt);
						   p++;
						   ptt++;
						}

						pt++;
					}

					pt += blockSize;
					p += blockSize;
				}
	}}}

	// C = A * B * At + Bt * Bt
	p = C;
	register double *pt = AxBxAt;
	register double *ptt = BtxBt;
	for(bi = 0; bi < N; bi += blockSize) {
        for(bj = 0; bj < N; bj += blockSize) {
            for(i = 0; i < blockSize; i++) {
				for(j = 0; j < blockSize; j++) {
					*p = *pt + *ptt;
					
					p++;
					pt++;
					ptt++;
				}
			}
	}}

	free(At);
	free(Bt);
	free(AxB);
	free(AxBxAt);
	free(BtxBt);

	return C;	
}
