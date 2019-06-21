#include "matrix.h"
#include <time.h>

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
    int i,j;
    int pos = 0;
    //A
   /*  
    int n = 8;
    double a[64] = {611,196,-192,407,-8,-52,-49,29,
                    196,899,113,-192,-71,-43,-8,-44,
                    -192,133,899,196,61,49,8,52,
                    407,-192,196,611,8,44,59,-23,
                    -8,-71,61,8,411,-599,208,208,
                    -52,-43,49,44,-599, 411,208,208,
                    -49,-8,8,59,208,208,99,-911,
                    29,-44,52,-23,208,208,-911, 99};
                    */
    //B
    /* 
    int n = 10;
    double a[100] = {0};
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            a[pos] = 1.0/(i+j+1);
            pos++;
        }
    }
    */
    //C
    /* 
    int n = 12;
    double a[144] = {0};
    double count;
    for(i=0;i<n;i++) {
        count = 0;
        for(j=0;j<n;j++) {
            double init = n-i;
            if(j<=i) {
                a[pos] = init;
            } else {
                a[pos] = init-(++count);
            }
            pos++;
        }
    }
    */

    //D
    int n = 20;
    double a[400] = {0};
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            a[pos] = sqrt(2 * 1.0 / 21) * sin((i+1)*(j+1)*PI/21);
            pos++;
        }
    }

    //show(newMatrix(n,n,a));
    struct timeval time_first;  
    mingw_gettimeofday(&time_first,NULL);
    double pld;
    double* env = (double*)malloc(sizeof(double) * n);
    for(i=0;i<n;i++) {
        env[i] = 1;
    }
    //power_eng(&pld, env, a, n);  //幂法
    //inv_power_eng(&pld, env, a, n); //反幂法
    jacobi_eng(env, a, n);  //jacobi
    printf("pld: %10.7f\n", pld); 
    show(newMatrix(n,1,env));
    struct timeval time_second;  
    mingw_gettimeofday(&time_second,NULL);
    printf("run time: %ldus\n", time_second.tv_usec - time_first.tv_usec);
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
    double eps = 1;
    int i,j,k;
    int iter_num = 0;
    while(eps > delta) {
        if(iter_num > T) {
            printf("limitation exceed\n");
            return true;
        }
        int p = 1, q = 1; 
        double mmax = 0;
        for(i=1;i<=n;i++) {
            for(j=1;j<=n;j++) {
                if(i!=j) {
                    double value = fabs(getData(mat,i,j));
                    if(value > mmax) {
                        mmax = value;
                        p = i;
                        q = j;
                    }
                }
            }
        }
        eps = mmax;
        double va, vt, vc, vs;
        va = (getData(mat,q,q) - getData(mat,p,p)) / (2 * getData(mat, p, q));
        if(va >= 0) {
            vt = 1.0 / (fabs(va) + sqrt(1+pow(va,2))); 
        } else {
            vt = -1.0 / (fabs(va) + sqrt(1+pow(va,2)));
        }
        vc = 1.0 / sqrt(1+pow(vt,2));
        vs = vt * vc;
        Matrix* ms = newMatrix(n,n,NULL);
        for(i=1;i<=n;i++) {
            for(j=1;j<=n;j++) {
                if(i==j) {
                    setData(ms,1,i,j);
                } else {
                    setData(ms,0,i,j);
                }
            }
        }
        setData(ms,vc,p,p);
        setData(ms,vs,p,q);
        setData(ms,-vs,q,p);
        setData(ms,vc,q,q);
        Matrix* ms_invert = invert(ms);
        Matrix* mt = multiplyMatrix(ms_invert, mat);
        mat = multiplyMatrix(mt, ms);
        //printf("%d time\n", iter_num+1);
        //show(mat);
        iter_num++;
    }
    for(i=0;i<n;i++) {
        ev[i] = getData(mat,i+1,i+1);
    }
    return false;
}