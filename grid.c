#include <stdlib.h>
#include "assert.h"
#include "grid.h"
#include "utils.h"


Grid *new_grid(int nx, int ny, double lX, double lY, int ghost_cells, double offsets[]){
    int i,j;
    Grid* newGrid = malloc(sizeof(Grid));
    assert(newGrid != NULL);
    newGrid->nx = nx;
    newGrid->ny = ny;
    newGrid->length_x = lX;
    newGrid->length_y = lY;
    newGrid->dx = lX/((double) nx);
    newGrid->dy = lY/((double) ny);
    newGrid->nghost = ghost_cells;
    newGrid->x = new_contiguous_2dArray(nx + 2*ghost_cells,ny + 2*ghost_cells);
    newGrid->y = new_contiguous_2dArray(nx + 2*ghost_cells,ny + 2*ghost_cells);

    for (i = 0; i < nx; i++){
        for (j = 0; j < ny; j++){
            newGrid->x[i + ghost_cells][j + ghost_cells] = offsets[0] + newGrid->dx*i;
            newGrid->y[i + ghost_cells][j + ghost_cells] = offsets[1] + newGrid->dy*j;
        }
    }

    // handle ghost zones specific to our problem. Not the best way to do it.
    for (i = 0; i < nx; i++){
        newGrid->x[i + ghost_cells][0] = offsets[0] + newGrid->dx*i;
        newGrid->y[i + ghost_cells][0] = offsets[1];

        newGrid->x[i + ghost_cells][ny + ghost_cells] = offsets[0] + newGrid->dx*i;
        newGrid->y[i + ghost_cells][ny + ghost_cells] = offsets[1] + lY;
    }

    for (j = 0; j <= ny; j++){
        newGrid->x[0][j + ghost_cells] = offsets[0];
        newGrid->y[0][j + ghost_cells] = offsets[1] + newGrid->dx*j;

        newGrid->x[nx + ghost_cells][j + ghost_cells] = offsets[0] + lY;
        newGrid->y[nx + ghost_cells][j + ghost_cells] = offsets[1] + newGrid->dy*j;
    }

    return newGrid;
}

void free_grid(Grid* grid){
    free_contiguous_2dArray(grid->x);
    free_contiguous_2dArray(grid->y);
    free(grid);
}

