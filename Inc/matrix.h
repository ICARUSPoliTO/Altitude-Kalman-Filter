#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>
#include <stdio.h>

#define MAX_SIZE 2 // Maximum size for matrices

/* Matrix Operation Functions */
void createIdentityMatrix(float I[MAX_SIZE][MAX_SIZE], int size);
void createDiagonalMatrix(float diag[MAX_SIZE][MAX_SIZE], float *vector, int size);
void zerosMatrix(float Z[MAX_SIZE][MAX_SIZE], int rows, int cols);
void addFloatMatrices(float A[MAX_SIZE][MAX_SIZE], float B[MAX_SIZE][MAX_SIZE],
                      float result[MAX_SIZE][MAX_SIZE], int rows, int cols);
void subFloatMatrices(float A[MAX_SIZE][MAX_SIZE], float B[MAX_SIZE][MAX_SIZE],
                      float result[MAX_SIZE][MAX_SIZE], int rows, int cols);
void multiplyFloatMatrices(float A[MAX_SIZE][MAX_SIZE], int rowsA, int colsA,
                           float B[MAX_SIZE][MAX_SIZE], int rowsB, int colsB,
                           float result[MAX_SIZE][MAX_SIZE]);
void transposeFloatMatrix(float src[MAX_SIZE][MAX_SIZE], float dest[MAX_SIZE][MAX_SIZE],
                          int rows, int cols);
float determinant(float matrix[MAX_SIZE][MAX_SIZE], int n);
void adjointMatrix(float matrix[MAX_SIZE][MAX_SIZE], float adj[MAX_SIZE][MAX_SIZE], int n);
int inverseMatrix(float matrix[MAX_SIZE][MAX_SIZE], float inverse[MAX_SIZE][MAX_SIZE], int n);

#endif