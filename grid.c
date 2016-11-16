#include <stdlib.h>
#include "grid.h"


Grid *new_grid(int nx, int ghost_cells, double length, double offsets[]){
    Grid* newGrid = malloc(sizeof(Grid));
    newGrid->nx = nx;
    newGrid->length = length;
    newGrid->dx = length/((double) nx);
    newGrid->nghost = ghost_cells;
    newGrid->data = malloc((nx + 2*ghost_cells)*sizeof(double));
    for(int i = 0; i < nx + 2*ghost_cells; i++){
        newGrid->data[i] = malloc((nx+2*ghost_cells)*sizeof(double));
    }
    return newGrid;
}

void free_grid(Grid* grid){
    for(int i = 0; i < grid->nx + 2*grid->nghost; i++){
        free(grid->data[i]);
    }
    free(grid->data);
    free(grid);
}

