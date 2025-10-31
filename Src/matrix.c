#ifndef MATRIX_C
#define MATRIX_C

#include "matrix.h"


void createIdentityMatrix(float I[MAX_SIZE][MAX_SIZE], int size)
{
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            I[i][j] = (i == j) ? 1.0f : 0.0f;
}

void createDiagonalMatrix(float diag[MAX_SIZE][MAX_SIZE], float *vector, int size)
{
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            diag[i][j] = (i == j) ? vector[i] : 0.0f;
}

void zerosMatrix(float Z[MAX_SIZE][MAX_SIZE], int rows, int cols)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            Z[i][j] = 0.0f;
}

void addFloatMatrices(float A[MAX_SIZE][MAX_SIZE], float B[MAX_SIZE][MAX_SIZE],
                      float result[MAX_SIZE][MAX_SIZE], int rows, int cols)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result[i][j] = A[i][j] + B[i][j];
}

void subFloatMatrices(float A[MAX_SIZE][MAX_SIZE], float B[MAX_SIZE][MAX_SIZE],
                      float result[MAX_SIZE][MAX_SIZE], int rows, int cols)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result[i][j] = A[i][j] - B[i][j];
}

void multiplyFloatMatrices(float A[MAX_SIZE][MAX_SIZE], int rowsA, int colsA,
                           float B[MAX_SIZE][MAX_SIZE], int rowsB, int colsB,
                           float result[MAX_SIZE][MAX_SIZE])
{
    if (colsA != rowsB)
        return; // dimension mismatch

    for (int i = 0; i < rowsA; i++)
        for (int j = 0; j < colsB; j++)
        {
            result[i][j] = 0.0f;
            for (int k = 0; k < colsA; k++)
                result[i][j] += A[i][k] * B[k][j];
        }
}

void transposeFloatMatrix(float src[MAX_SIZE][MAX_SIZE], float dest[MAX_SIZE][MAX_SIZE],
                          int rows, int cols)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            dest[j][i] = src[i][j];
}

float determinant(float matrix[MAX_SIZE][MAX_SIZE], int n)
{
    if (n == 1)
        return matrix[0][0];

    if (n == 2)
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    float sub[MAX_SIZE][MAX_SIZE];
    float det = 0.0f;

    for (int p = 0; p < n; p++)
    {
        int subi = 0;
        for (int i = 1; i < n; i++)
        {
            int subj = 0;
            for (int j = 0; j < n; j++)
            {
                if (j == p)
                    continue;
                sub[subi][subj] = matrix[i][j];
                subj++;
            }
            subi++;
        }
        det += matrix[0][p] * powf(-1, p) * determinant(sub, n - 1);
    }
    return det;
}

void adjointMatrix(float matrix[MAX_SIZE][MAX_SIZE], float adj[MAX_SIZE][MAX_SIZE], int n)
{
    if (n == 1)
    {
        adj[0][0] = 1.0f;
        return;
    }

    float temp[MAX_SIZE][MAX_SIZE];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int subi = 0;
            for (int row = 0; row < n; row++)
            {
                if (row == i)
                    continue;
                int subj = 0;
                for (int col = 0; col < n; col++)
                {
                    if (col == j)
                        continue;
                    temp[subi][subj] = matrix[row][col];
                    subj++;
                }
                subi++;
            }
            adj[j][i] = powf(-1, i + j) * determinant(temp, n - 1);
        }
    }
}

int inverseMatrix(float matrix[MAX_SIZE][MAX_SIZE], float inverse[MAX_SIZE][MAX_SIZE], int n)
{
    float det = determinant(matrix, n);
    if (fabs(det) < 1e-6f)
        return 0; // singular matrix, no inverse

    float adj[MAX_SIZE][MAX_SIZE];
    adjointMatrix(matrix, adj, n);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inverse[i][j] = adj[i][j] / det;

    return 1;
}

#endif