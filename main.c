#include "matrix.h"

typedef int bool;

#define true 1
#define false 0
#define delta 0.00001
#define T     1000000
#define PI    3.1415926

bool power_eng(double* pld, double* env, double* a, int n);
bool lu(double* a, int* pivot, int n);
bool jacobi_eng(double* ev, double* a, int n);
bool inv_power_eng(double* pld, double* env, double* a, int n);

int main(int argc, char** argv) {
    //jacobi
    /* 
    double a[9]={1,2,5,0.5,1,2,0.2,0.5,1};
    double env[3] = {0};
    jacobi_eng(env, a, 3);
    show(newMatrix(3,1,env));
    */
     //反幂法
    double a[9]={1,1,0.5,1,1,0.25,0.5,0.25,2};
    double pld;
    double env[3] = {1,1,1};
    inv_power_eng(&pld,env,a,3);
    printf("pld: %.2f\n", pld);
    show(newMatrix(3,1,env));
     //幂法
     /* 
    double a[9]={1,2,5,0.5,1,2,0.2,0.5,1};
    double pld;
    double env[3] = {1,1,1};
    power_eng(&pld, env, a, 3);
    printf("pld: %.2f\n", pld);
    show(newMatrix(3,1,env));
    */
    //LU分解 
    /* 
    double l[9] = {0};
    double u[9] = {0};
    double a[9] = {1,2,3,2,5,7,3,5,3};
    lu(a,l,u, 3);
    Matrix* lm = newMatrix(3,3,l);
    Matrix* um = newMatrix(3,3,u);
    show(lm);
    show(um);
    show(multiplyMatrix(lm,um));
    */
   //矩阵测试
    /* 
     double* a = (double*)malloc(sizeof(double) * 8);
     int i;
     for(i=1;i<=8;i++) {
         a[i] = i;
     }
     Matrix* ma = newMatrix(2,4,a);
     show(ma);
     Matrix* in = invert(ma);
     show(in);
     double max = getMax(ma);
     show(multiplyDouble(ma, max));
     show(multiplyDouble(ma, getMin(ma)));
     show(multiplyMatrix(ma, in));
     show(addMatrix(ma, ma, 0));
     show(addMatrix(ma, ma, 1));
     printf("%.2f\n", getMod(ma));
     freeMatrix(ma);
     */
}

bool lu(double* a, int* pivot, int n) {
	int i, j, k;
	double max = 0, temp = 0.0;      
	for (i = 0; i < n - 1; i++){       
		max = fabs(a[n*i + i]);         
		pivot[i]=i; 
		for (j = i + 1; j < n; j++){             
			if (fabs(a[n*j  + i]) > max) {            
				max = fabs(a[n*j + i]);
				pivot[i] = j;
			}
		}
		if (fabs(a[n*n - 1]) < 1e-14) {
			printf("Matrix LU decomposition failed!\n");
			return true;
		}
		if(pivot[i]!=i) {         
			for (j = i; j < n; j++){       
				temp = a[n*i  + j];
				a[n*i  + j] = a[n*pivot[i] + j];                 
				a[n*pivot[i] + j] = temp;
			}
		}
		for (j = i + 1; j < n; j++)             
			a[n*j + i]=a[n*j + i]/a[n*i + i]; 
		for (j = i + 1; j < n; j++)             
			for(k=i+1; k<n; k++)  
				a[n*j  + k] = a[n*j + k] - a[n*j + i] * a[n*i + k];     
	}
	for (i = 0; i < n - 2; i++)          
		for(k=i+1; k<n-1;k++) {          
			temp = a[n*pivot[k] + i];
			a[n*pivot[k] + i] = a[k*n  + i];             
			a[k*n  + i] = temp;         
		}
	if (fabs(a[n*n - 1]) < 1e-14) {
		printf("Matrix LU decomposition failed!\n");
		return true;
	}else
		return false;
}

bool guass(double const* lu, int const* p, double* b, int n) {
	int i, j;
	double temp;
	for (i = 0; i < n; i++){
		if (fabs(lu[i*n + i]) < 1e-14) {
			printf("Matrix singularity!\n");
			return true;
		}
	}

	for (i = 0; i < n - 1; i++) {      
		temp  = b[p[i]];         
		b[p[i]] = b[i];         
		b[i] = temp;
	}
	 
	for(i=0; i<n; i++)          
		for(j=0; j<i; j++) 
			b[i] = b[i] - lu[n*i + j] * b[j];       
	
	for (i = n - 1; i >= 0; i--)  {    
		for (j = n - 1; j > i; j--)
			b[i] = b[i] - lu[n*i + j] * b[j];          
		b[i] = b[i] / lu[n*i + i];
	}
	return false;
}

bool inv_power_eng(double* pld, double* env, double* a, int n) {
    Matrix* eng = newMatrix(n,1,env);
    Matrix* enk;
    Matrix* temp;
    int* pivot = (int*)malloc(sizeof(int) * n);
    lu(a,pivot,n);
    int i,j;
    int iter_num = 0;
    double eps = 1;
    while(eps > delta) {
        if(iter_num > T) {
            printf("inv_power_eng limitation exceed\n");
            return true;
        }
        temp = eng;
        enk = multiplyDouble(eng, 1.0 / getAbsMax(eng));
        guass(a,pivot,convertToArray(enk),n);
        eng = enk;
        eps = fabs(getAbsMax(eng) - getAbsMax(temp));
        iter_num++;
    }
    double* data = convertToArray(multiplyDouble(temp, 1.0 / getAbsMax(temp)));
    for(i=0;i<n;i++) {
        env[i] = data[i];
    }
    *pld = 1.0 / getMax(eng);
    return false;

}

bool power_eng(double* pld, double* env, double* a, int n) {
    Matrix* mat = newMatrix(n,n,a);
    Matrix* eng = newMatrix(n,1,env);
    Matrix* enk;
    int iter_num = 0;
    double eps = 1;
    while(eps >= delta) {
        if(iter_num > T) {
            return true;
        }
        Matrix* temp = eng;
        enk = multiplyDouble(eng, 1.0 / getAbsMax(eng));
        eng = multiplyMatrix(mat,enk);
        eps = fabs(getAbsMax(eng) - getAbsMax(temp));
        iter_num++;
    }
    double* data = convertToArray(enk);
    int i;
    for(i=0;i<n;i++) {
        env[i] = data[i];
    }
    *pld = getMax(eng);
    return false;
}

bool jacobi_eng(double* ev, double* a, int n) {
     Matrix* mat = newMatrix(n,n,a);
     Matrix* eng = newMatrix(n,n,NULL);
     int i,j;
     int iter_num = 0;
     double max = getData(mat, 1, 2);
     int row = 1;
     int col = 2;
     for(i=1;i<=n;i++) {
         for(j=1;j<=n;j++) {
             if(i != j) {
                 setData(eng,0,i,j);
                 double data = getData(mat, i, j);
                 if(data > max) {
                     max = data;
                     row = i;
                     col = j;
                 }
             } else {
                 setData(eng,1,i,j);
             } 
         }
     }
     while(max > delta) {
         if(iter_num > T) {
             return true;
         }
         double theta;
         double data = getData(mat,row,col);
         double rowV = getData(mat, row, row);
         double colV = getData(mat, col, col);
         if(rowV == colV) {
             theta = PI / 4;
         } else {
             theta = 0.5 * atan(2.0 * data / (rowV - colV));
         }
         double sin_theta = sin(theta);
         double cos_theta = cos(theta);
         double sin_2theta = sin(2*theta);
         double cos_2theta = cos(2*theta);
         setData(mat, rowV * pow(cos_theta,2) + colV * pow(sin_theta, 2) + data * sin_2theta, row, row);
         setData(mat, rowV * pow(cos_theta,2) + colV * pow(sin_theta, 2) - data * sin_2theta, col, col);
         setData(mat, 0.5 * (colV - rowV) * sin_2theta + data * cos_2theta, row, col);
         setData(mat, getData(mat, row, col), col, row);
         int k;
         for(k=1;k<=n;k++) {
             if(k != row && k!= col) {
                 double rowk = getData(mat, row, k);
                 double colk = getData(mat, col, k);
                 setData(mat, rowk * cos_theta + colk * sin_theta, row, k);
                 setData(mat, getData(mat, row, k), k, row);
                 setData(mat, colk * cos_theta - rowk * sin_theta, col, k);
                 setData(mat, getData(mat, col, k), k, col);
             }
         }
         if(iter_num == 0) {
             setData(eng, cos_theta, row, row);
             setData(eng, -sin_theta, row, col);
             setData(eng, sin_theta, col, row);
             setData(eng, cos_theta, col, col);
         } else {
             for(k=1;k<=n;k++) {
                 double prow = getData(eng,k,row);
                 double pcol = getData(eng,k,col);
                 setData(eng,prow * cos_theta + pcol * sin_theta, k, row);
                 setData(eng,prow * cos_theta - pcol * sin_theta, k, col);
             }
         }
         for(i=1;i<=n;i++) {
             for(j=1;j<=n;j++) {
                 if(i != j) {
                     double data = getData(mat, i, j);
                     if(data > max) {
                         max = data;
                         row = i;
                         col = j;
                     }
                 } 
             }
         }
         iter_num++;
     }
     for(i=0;i<n;i++) {
         ev[i] = getData(mat, i+1, i+1); 
     }
}