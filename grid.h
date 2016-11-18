/*
 * Grid : Structure that represents the grid on which fields live.
 * 	      Contains grid spacing, lenghts and information on the ghost cells. 
 */
#ifndef GRID_H
#define GRID_H


typedef struct grid_t Grid;

struct grid_t{
        double **x;
        double **y;
        int nx;
        int ny;
        double length_x;
        double length_y;
        double dx;
        double dy;
        int nghost;
};

Grid *new_grid(int nx, int ny, double length_x, double length_y, int ghost_cells, double offsets[]);

void free_grid(Grid * grid);

static inline double grid_X(Grid* grid, int i, int j){
    return grid->x[i+grid->nghost][j+grid->nghost];
}

static inline double grid_Y(Grid* grid, int i, int j){
    return grid->y[i+grid->nghost][j+grid->nghost];
}



#endif
