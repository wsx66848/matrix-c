#ifndef MY_MATRIX_H
#define MY_MATRIX_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

struct _Matrix {
    double* data;
    int row;
    int col;
};

typedef struct _Matrix Matrix;

Matrix* newMatrix(int row, int col,double* data);
void    freeMatrix(Matrix* matrix);
int     getRow(Matrix* matrix);
int     getCol(Matrix* matrix);
double  getData(Matrix* matrix, int row, int col);
void    setData(Matrix* matrix, double value, int row, int col);
Matrix* invert(Matrix* matrix);
Matrix* multiplyDouble(Matrix* matrix, double value);
Matrix* multiplyMatrix(Matrix* matrix1, Matrix* matrix2);
double  getMod(Matrix* matrix);
Matrix* addMatrix(Matrix* matrix1, Matrix* matrix2, int operator);
double  getMax(Matrix* matrix);
double  getMin(Matrix* matrix);
void    show(Matrix* matrix);
double* convertToArray(Matrix* matrix);

#endif