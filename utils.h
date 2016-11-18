/*
 * Utils : Helpful utilities.
 *
 * 			Specifically: constructor of 2D array that's guaranteed to be contiguous
 * 			              in memory (unnecessary?).
 */
#ifndef UTILS_H
#define UTILS_H


double **new_contiguous_2dArray(int nx, int ny);

void free_contiguous_2dArray(double **array);

#endif
