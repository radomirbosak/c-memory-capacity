#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h> // for memcpy
//#include <stdbool.h>


#include "foo.h"
#include "randomvar.h"
#include "matrix.h"

#define ITERATIONS 2000
#define ITERATIONS_SKIPPED 1000
#define ITERATIONS_COEF_MEASURE 1000
#define RUNS 1
#define INPUT_LOW -1.0
#define INPUT_HIGH 1.0


int main(int argc, char** argv) {
	printf("-- Doing one computation of MC! --\n");
	//printf("ste hej som");

	clock_t start, end;
	double cpu_time_used;


	int q = 100;
	double W[q * q];
	double WI[q];
	
	double tau = 0.1;
	double sigma = 0.095;

	initrand();
	int i;
	for (i = 0; i < q*q; i++)
		W[i] = normal(0.0, sigma);

	for (i = 0; i < q; i++) 
		WI[i] = unif(-tau, tau);

	double mc;
	
	start = clock();
	for (i = 0; i < RUNS; i++)
		mc = MC(W, WI, sigma, tau, q);

	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC / ((double)RUNS);
	printf("time used: %f\n", cpu_time_used);
	printf("MC = %f\n", mc);

	return 0;
}

void print_vec(double* vec, int size) {
	int i;
	printf("[");
	for (i = 0; i < size; i++)
		printf("%f, ", vec[i]);

	printf("]\n");
}

void print_mat(double* mat, int rows, int cols) {
	int i,j;
	printf("%d x %d matrix\n", rows, cols);
	printf("[");
	for (j = 0; j < rows; j++) {
		print_vec(mat + j * cols, cols);
	}
	printf("]\n");
}

clock_t start, end;
int time_id;
void time_start() {
	time_id = 0;
	start = clock();
}

void time_check() {

	double cpu_time_used;
	cpu_time_used = ((double) (clock() - start)) / CLOCKS_PER_SEC;
	printf("time #%d: %f\n", time_id, cpu_time_used);

	time_id++;
	start = clock();
}

double MC(double* W, double* WI, double sigma, double tau, int memory_max) {

	int ITERATIONS_MEASURED = ITERATIONS - ITERATIONS_SKIPPED;
	 

	int q = 100, p = 1;
	int i,h; // for iteration

	assert(ITERATIONS >= memory_max + ITERATIONS_COEF_MEASURE);

	double* u = newVector(ITERATIONS);
	for (i = 0; i < ITERATIONS; i++)
		u[i] = unif01();

	double* X = newVectorZ(q);
	double* X1 = newVectorZ(q);
	double* S = newMatrixZ(ITERATIONS_MEASURED, q); //  vymenene: vela riadkov, 100 stlpcov
	double* Sbeh = S;

	time_start();

	printf("First run\n");
	for (i = 0; i < ITERATIONS; i++) {
		dotMV(W, X, q, q, X1);
		dotVS(WI, u[i], q, X);
		plusVV(X, X1, q, X);
		tanhV(X, q, X); 

		if (i >= ITERATIONS_SKIPPED) {
			memcpy(Sbeh, X, q * sizeof(double)); // priestor na optimalizaciu
			Sbeh = Sbeh + q;
		}
	}

	time_check();

	assert(memory_max < ITERATIONS_SKIPPED);
	double* D = newMatrix(memory_max, ITERATIONS_MEASURED);

	// tu sa da tiez asi optimalizovat (D sa bude nasobit)
	for (i = 0; i < memory_max; i++) 
		memcpy(D + i * ITERATIONS_MEASURED, u + (ITERATIONS_SKIPPED - i), ITERATIONS_MEASURED * sizeof(double));
	
	double* S_T = newMatrix(q, ITERATIONS_MEASURED);
	transpose(S, ITERATIONS_MEASURED, q, S_T);
	

	time_check();

	printf("Computing pinverse...\n");
	double* S_PINV_T = newMatrix(ITERATIONS_MEASURED, q);
	pinv(S_T, q, ITERATIONS_MEASURED, S_PINV_T);
	
	time_check();


	double* WO = newMatrix(memory_max, q);
	dotMM(D, S_PINV_T, memory_max, ITERATIONS_MEASURED, q, WO);
	
	time_check();

	//print_mat(D, memory_max, ITERATIONS_MEASURED);
	// measure memory capacity
	// let u be
	//for (i = 0; i < ITERATIONS; i++)
	//	u[i] = unif01();
	printf("Output run & correlation coefficient\n");
	memset(X, 0, q * sizeof(double));

	double* o = newMatrixZ(ITERATIONS_COEF_MEASURE, memory_max); // vymenene
	double* obeh = o;
	// new run, but we use the old inputs
	for (i = 0; i < ITERATIONS_COEF_MEASURE + memory_max; i++) {
		dotMV(W, X, q, q, X1);
		dotVS(WI, u[i], q, X);
		plusVV(X, X1, q, X);
		tanhV(X, q, X);

		if (i >= memory_max) {
			dotMV(WO, X, memory_max, q, obeh);
			obeh = obeh + memory_max;


		}
	}

	/*
	for (i = 0; i < q; i++) 
		printf("o[%d] = %f\n", i, o[i + (ITERATIONS_COEF_MEASURE - 1) * q]);/**/

	// medzi u[minule] a o[terajsie?]
	double mcsum = 0, mymc;
	for (h = 0; h < memory_max; h++) {
		// vypocitat corrcoef z
		// u[memory_max - h : memory_max + iterations_coef_measure - h]
		// a
		// o[h, : ]))
		// corrcoef - space for optimization
		mymc = corrcoefVTV(o + h, u + (memory_max - h), ITERATIONS_COEF_MEASURE, memory_max);
		//printf("%f\n", mymc);
		mcsum +=mymc;

	}

	time_check();
	printf("\n");

	free(o);      // (ITERATIONS_COEF_MEASURED, memory_max) [op]
	free(WO);     // (memory_max, q)
	free(S_PINV_T);
	free(S_T);
	free(D);
	free(u);
	free(X);      // (q)
	free(S);
	return mcsum;
}
