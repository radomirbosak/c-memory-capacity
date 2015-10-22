#ifndef __MATRIX_H__
#define __MATRIX_H__
#endif

#include <stdbool.h>

double dotVV(double* vec1, double* vec2, int size);
void dotMV(double* mat, double* vec, int rows, int cols, double* result);
void dotMM(double* mat1, double* mat2, int rows1, int cols1, int cols2, double* result);
void dotVS(double* vec, double factor, int size, double* result);

void plusVV(double* vec1, double* vec2, int size, double* result);
void plusMM(double* mat1, double* mat2, int rows, int cols, double* result);
void timesVV(double* vec1, double* vec2, int size, double* result);
void timesMM(double* mat1, double* mat2, int rows, int cols, double* result);

bool equalsVV(double* vec1, double* vec2, int size);
bool equalsMM(double* mat1, double* mat2, int rows, int cols);

double* newMatrix(int rows, int cols);
double* newVector(int size);
double* newMatrixZ(int rows, int cols);
double* newVectorZ(int size);

void tanhV(double* vec, int size, double* result);

void pinv(double* mat, int rows, int cols, double* result);
void transpose(double* mat, int rows, int cols, double* result);

double corrcoef(double* vec1, double* vec2, int size);

double corrcoefVTV(double* vec1, double* vec2, int size, int rowskip);