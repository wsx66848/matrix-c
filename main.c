#include "matrix.h"

int main(int argc, char** argv) {
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
}