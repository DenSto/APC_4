#include <stdlib.h>
#include "assert.h"
#include "field.h"
#include "utils.h"


Field *new_field(int nx, int ny, int ghost_cells){
    int i,j;
    Field* newField = malloc(sizeof(Field));
    assert(newField!=NULL);
    newField->nx = nx;
    newField->ny = ny;
    newField->nghost = ghost_cells;
    newField->data = new_contiguous_2dArray(nx + 2*ghost_cells,ny + 2*ghost_cells);
    
    for(i = 0; i < nx + 2*ghost_cells; i++){
        for(j = 0; i < ny + 2*ghost_cells; i++){
            newField->data[i][j] = 0.0;
        }
    }
    return newField;
}

void free_field(Field* field){
    free_contiguous_2dArray(field->data);
    free(field);
}

void set_value(Field* field, int i, int j, double value){
    field->data[i+field->nghost][j+field->nghost] = value;
}

double get_value(Field* field, int i, int j){
    return field->data[i+field->nghost][j+field->nghost];
}
