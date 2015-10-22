
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

double dotVV(double* vec1, double* vec2, int size) {
	double sum = 0;
	int i;
	for (i = 0; i < size; i++)
		sum = sum + vec1[i] * vec2[i];
	return sum;
}

void dotMV(double* mat, double* vec, int rows, int cols, double* result) {
	double sum;
	int i,j; 
	for (j = 0; j < rows; j++) {
		sum = 0;
		for (i = 0; i < cols; i++)
			sum = sum + mat[j * cols + i] * vec[i];
		result[j] = sum;
	}
}

void dotMM(double* mat1, double* mat2, int rows1, int cols1, int cols2, double* result) {
	int i,j,k;
	double sum;
	for (j = 0; j < rows1; j++) {
		for (i = 0; i < cols2; i++) {
			sum = 0;
			for (k = 0; k < cols1; k++) {
				sum = sum + mat1[j * cols1 + k] * mat2[i + k * cols2];
			}
			result[j * cols2 + i] = sum;
		}
	}
}

void dotVS(double* vec, double factor, int size, double* result) {
	int i;
	for (i = 0; i < size; i++)
		result[i] = vec[i] * factor;
}

void plusVV(double* vec1, double* vec2, int size, double* result) {
	int i;
	for (i = 0; i < size; i++)
		result[i] = vec1[i] + vec2[i];
}

void plusMM(double* mat1, double* mat2, int rows, int cols, double* result) {
	int i;
	for (i = 0; i < rows*cols; i++)
		result[i] = mat1[i] + mat2[i];
}

void timesVV(double* vec1, double* vec2, int size, double* result) {
	int i;
	for (i = 0; i < size; i++)
		result[i] = vec1[i] * vec2[i];
}

void timesMM(double* mat1, double* mat2, int rows, int cols, double* result) {
	int i;
	for (i = 0; i < rows*cols; i++)
		result[i] = mat1[i] * mat2[i];
}

bool equalsVV(double* vec1, double* vec2, int size) {
	int i;
	for (i = 0; i < size; i++)
		if (vec1[i] != vec2[i])
			return false;
	return true;
}

bool equalsMM(double* mat1, double* mat2, int rows, int cols) {
	return equalsVV(mat1, mat2, rows*cols);
}



/*
 vectorized functions
*/

void tanhV(double* vec, int size, double* result) {
	int i;
	for (i = 0; i < size; i++)
		result[i] = tanh(vec[i]);
}

/*
 Allocation functions
*/

double* newMatrix(int rows, int cols) {
	return (double*)malloc(rows * cols * sizeof(double));
}

double* newVector(int size) {
	return (double*)malloc(size * sizeof(double));
}

double* newMatrixZ(int rows, int cols) {
	return (double*)calloc(rows * cols, sizeof(double));
}

double* newVectorZ(int size) {
	return (double*)calloc(size, sizeof(double));
}

/*
 Other
*/


void transpose(double* mat, int rows, int cols, double* result) {
	int i,j;
	for (j=0; j < rows; j++) {
		for (i = 0; i < cols; i++) {
			result[i * rows + j] = mat[j * cols + i];
		}
	}
}

double corrcoef(double* vec1, double* vec2, int size) {
	double n = (double)size;
	double xsum = 0, ysum = 0, xysum = 0, xsum2 = 0, ysum2 = 0;
	double x,y;
	int i;
	for (i = 0; i < size; i++) {
		x = vec1[i]; y = vec2[i];
		xsum += x;
		xsum2 += x*x;
		ysum += y;
		ysum2 += y*y;
		xysum += x*y;
	}
	return (n * xysum - xsum * ysum) / 
		(sqrt(n*xsum2 - xsum*xsum) * sqrt(n*ysum2 - ysum*ysum));
}

double corrcoefVTV(double* vec1, double* vec2, int size, int rowskip) {
	double n = (double)size;
	double xsum = 0, ysum = 0, xysum = 0, xsum2 = 0, ysum2 = 0;
	double x,y;
	double* vec1beh = vec1;
	int i;
	for (i = 0; i < size; i++) {
		x = *vec1beh; y = vec2[i];
		xsum += x;
		xsum2 += x*x;
		ysum += y;
		ysum2 += y*y;
		xysum += x*y;

		vec1beh += rowskip;
	}
	return (n * xysum - xsum * ysum) / 
		(sqrt(n*xsum2 - xsum*xsum) * sqrt(n*ysum2 - ysum*ysum));
}

void pinv(double* mat, int rows, int cols, double* result) {
	// if A = U S V*
	// => A^ = V S^ U* 
	// where S^ has reciprocals on diagonal, except zeros- they stay

	// int n_fil = 3;
	// int n_col = 8;
	// double mat[3][8] =  {{ 50,4.5, -23,  12,  1,  0, -1,  0},		// Example Rank-deficient matrix
	// 	             	 {  1,  2,   3,   4,  5, 1, 0, 0},
	//     	             {  2,  4,   6,   8, 10, 2, 0, 0}};

	int n_fil = rows;
	int n_col = cols;

	unsigned i = 0;
	unsigned j = 0;
	gsl_matrix * gA = gsl_matrix_alloc (n_fil, n_col);
	for (i = 0; i < n_fil; i++)
		for (j = 0; j < n_col; j++)
	   		gsl_matrix_set (gA, i, j, (double)mat[i * cols + j]);


	gsl_matrix * gA_t = gsl_matrix_alloc (n_col, n_fil);
	gsl_matrix_transpose_memcpy (gA_t, gA);					// Computing the transpose of gA
		
	gsl_matrix * U = gsl_matrix_alloc (n_col, n_fil);
	gsl_matrix * V= gsl_matrix_alloc (n_fil, n_fil);
	gsl_vector * S = gsl_vector_alloc (n_fil);


	// Computing the SVD of the transpose of A
	// The matrix 'gA_t' will contain 'U' after the function is called
	gsl_vector * work = gsl_vector_alloc (n_fil);
	gsl_linalg_SV_decomp (gA_t, V, S, work);
	gsl_vector_free(work);
		
	gsl_matrix_memcpy (U, gA_t);
	
				
	//Inverting S//
	//----------------------------------------------------------
	// Matrix 'S' is diagonal, so it is contained in a vector.
	// We operate to convert the vector 'S' into the matrix 'Sp'.
	//Then we invert 'Sp' to 'Spu'
	//----------------------------------------------------------
	gsl_matrix * Sp = gsl_matrix_alloc (n_fil, n_fil);
	gsl_matrix_set_zero (Sp);
	for (i = 0; i < n_fil; i++)
		gsl_matrix_set (Sp, i, i, gsl_vector_get(S, i));	// Vector 'S' to matrix 'Sp'
	
	gsl_permutation * p = gsl_permutation_alloc (n_fil);
	int signum;
	gsl_linalg_LU_decomp (Sp, p, &signum);				// Computing the LU decomposition
	
	gsl_matrix * SI = gsl_matrix_calloc (n_fil, n_fil);

	for (i = 0; i < n_fil; i++) {
	  //std::cout << "S [" << i << "] = " << gsl_vector_get (S, i) << std::endl;

	  if (gsl_vector_get (S, i) > 0.0000000001)
	    gsl_matrix_set (SI, i, i, 1.0 / gsl_vector_get (S, i));
	}
	//end Inverting S//
		
		
		
	gsl_matrix * VT = gsl_matrix_alloc (n_fil, n_fil);
	gsl_matrix_transpose_memcpy (VT, V);					// Tranpose of V
		
		
	//THE PSEUDOINVERSE//
	//----------------------------------------------------------
	//Computation of the pseudoinverse of trans(A) as pinv(A) = U·inv(S).trans(V)   with trans(A) = U.S.trans(V)
	//----------------------------------------------------------
	gsl_matrix * SIpVT = gsl_matrix_alloc (n_fil, n_fil);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,				// Calculating  inv(S).trans(V)
                	1.0, SI, VT,
                	0.0, SIpVT);

			
	gsl_matrix * pinvm = gsl_matrix_alloc (n_col, n_fil);	// Calculating  U·inv(S).trans(V)
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                	1.0, U, SIpVT,
                	0.0, pinvm);

	//end THE PSEUDOINVERSE//
	
   	for (i = 0; i < n_col; i++)  
    	for (j = 0; j < n_fil; j++)
    		result[i * n_fil + j] = (double) gsl_matrix_get (pinvm, i, j);

	gsl_matrix_free(VT);
	gsl_matrix_free(SI);
	gsl_matrix_free(SIpVT);
	gsl_matrix_free(gA_t);
	gsl_matrix_free(U);
	gsl_matrix_free(gA);
	gsl_matrix_free(V);
	gsl_vector_free(S);

	return;
}