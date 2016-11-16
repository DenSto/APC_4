#ifndef UTILS_H
#define UTILS_H


double **new_contiguous_2dArray(int nx, int ny);

void free_contiguous_2dArray(double **array);

void shape_1d_to_2d(double *in, double **out, int nx, int ny);


#endif
