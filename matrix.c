#include "matrix.h"

Matrix* newMatrix(int row, int col, double* data) {
    Matrix* matrix = (Matrix*)malloc(sizeof(Matrix));
    if(data != NULL) {
      matrix->data = data;
      matrix->row = row;
      matrix->col = col;
    } else {
      matrix->data = (double*)malloc(sizeof(double) * row * col);
      matrix->row = row;
      matrix->col = col;
    }
    return matrix;
}

void freeMatrix(Matrix* matrix) {
    free(matrix->data);
    free(matrix);
}

int getRow(Matrix* matrix) {
    return matrix->row;
}

int getCol(Matrix* matrix) {
    return matrix->col;
}

Matrix* invert(Matrix* matrix) {
    int row = getRow(matrix);
    int col = getCol(matrix);
    Matrix* res = newMatrix(col, row, NULL);
    int i = 1, j = 1;
    for(i=1;i<=row;i++) {
        for(j=1;j<=col;j++) {
            setData(res, getData(matrix, i, j), j, i);
        }
    }
    return res;
}

double getData(Matrix* matrix, int row, int col) {
    double* data = matrix->data;
    int m_col = getCol(matrix);
    return data[(row-1)*m_col + (col-1)];
}

void setData(Matrix* matrix, double value, int row, int col) {
    double* data = matrix->data;
    int m_col = getCol(matrix);
    data[(row-1)*m_col + (col-1)] = value;
}

Matrix* multiplyMatrix(Matrix* matrix1, Matrix* matrix2) {
    int row = getRow(matrix1);
    int col = getCol(matrix2);
    int m_count = getCol(matrix1);
    Matrix* res = newMatrix(row, col, NULL);
    int i = 1, j = 1, k=1;
    for(i=1;i<=row;i++) {
        for(j=1;j<=col;j++) {
            double temp = 0;
            for(k=1;k<=m_count;k++) {
                temp += getData(matrix1,i,k) * getData(matrix2, k, j);
            }
            setData(res, temp, i, j);
        }
    }
    return res;
}

Matrix* multiplyDouble(Matrix* matrix, double value) {
    int row = getRow(matrix);
    int col = getCol(matrix);
    Matrix* res = newMatrix(row,col,NULL);
    int i = 1, j = 1;
    for(i=1;i<=row;i++) {
        for(j=1;j<=col;j++) {
            setData(res, getData(matrix, i, j) * value, i, j);
        }
    }
    return res;
}

double getMod(Matrix* matrix) {
    int row = getRow(matrix);
    int col = getCol(matrix);
    double res = 0;
    int i = 1, j = 1;
    for(i=1;i<=row;i++) {
        for(j=1;j<=col;j++) {
            res += pow(getData(matrix, i, j), 2);
        }
    }
    return sqrt(res);
}

Matrix* addMatrix(Matrix* matrix1, Matrix* matrix2, int operator) {
    int row = getRow(matrix1);
    int col = getCol(matrix2);
    int i = 1,j = 1;
    Matrix* res = newMatrix(row, col, NULL);
    for(i=1;i<=row;i++) {
        for(j=1;j<=col;j++) {
            if(operator == 0) {
                setData(res,getData(matrix1,i,j) + getData(matrix2, i, j), i, j);
            } else {
                setData(res,getData(matrix1,i,j) - getData(matrix2, i, j), i, j);
            }
        }
    }
    return res;
}

double getMax(Matrix* matrix) {
    int max = getData(matrix,1,1);
    int row = getRow(matrix);
    int col = getCol(matrix);
    int i = 1, j = 1;
    for(i=1;i<=row;i++) {
        for(j=1;j<=col;j++) {
            double value = getData(matrix, i, j);
            if(value >= max) {
                max = value; 
            }
        }
    }
    return max;
}

double getMin(Matrix* matrix) {
    int min = getData(matrix,1,1);
    int row = getRow(matrix);
    int col = getCol(matrix);
    int i = 1, j = 1;
    for(i=1;i<=row;i++) {
        for(j=1;j<=col;j++) {
            double value = getData(matrix, i, j);
            if(value < min) {
                min = value; 
            }
        }
    }
    return min;
}

double* convertToArray(Matrix* matrix) {
    return matrix->data;
}

void show(Matrix* matrix) {
    int row = getRow(matrix);
    int col = getCol(matrix);
    int i = 1, j = 1;
    for(i=1;i<=row;i++) {
        for(j=1;j<=col;j++) {
            printf("%.2f ", getData(matrix, i, j));
        }
        printf("\n");
    }
}
