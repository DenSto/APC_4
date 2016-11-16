#ifndef GRID_H
#define GRID_H


typedef struct grid_t Grid;

struct grid_t{
        double **data;
        int nx;
        double length;
        double dx;
        int nghost;
};

Grid *new_grid(int nx, int ghost_cells, double length, double offsets[]);

void free_grid(Grid * grid);


#endif
