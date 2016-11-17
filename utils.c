#include <stdlib.h>
#include "assert.h"

double **new_contiguous_2dArray(int nx, int ny){
    double *newarr;
    double **ret;
    int i;

    newarr = malloc(nx*ny*sizeof(double));
    ret= malloc(nx*sizeof(double));
    assert(newarr != NULL && ret != NULL);

    for(i = 0; i < nx; i++){
        ret[i] = &(newarr[i*ny]);
    }
    return ret;
}

void free_contiguous_2dArray(double **array){
    free(array[0]);
    free(array);
}

void shape_1d_to_2d(double *in, double **out, int nx, int ny){

}
