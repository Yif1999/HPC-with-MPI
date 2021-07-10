#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void Block_Divide(int n,int *a){
    double ref=sqrt(n);
    int refINT=int(ceil(ref));
    for (int i=refINT;i>0;i--){
        if (1.0*n/i==int(1.0*n/i)){
            if (i<n/i){
                a[0]=i;
                a[1]=n/i;
            }else{
                a[1]=i;
                a[0]=n/i;
            }
            return;
        }
    }
}