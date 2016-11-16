#include <stdlib.h>
#include "field.h"


Field *new_field(int nx, int ghost_cells){
    Field* newField = malloc(sizeof(Field));
    newField->nx = nx;
    newField->nghost = ghost_cells;
    newField->data = malloc((nx + 2*ghost_cells)*sizeof(double));
    for(int i = 0; i < nx + 2*ghost_cells; i++){
        newField->data[i] = malloc((nx+2*ghost_cells)*sizeof(double));
    }
    
    for(int i = 0; i < nx + 2*ghost_cells; i++){
        for(int j = 0; i < nx + 2*ghost_cells; i++){
            newField->data[i][j] = 0.0;
        }
    }
    return newField;
}

void free_field(Field* field){
    for(int i = 0; i < field->nx + 2*field->nghost; i++){
        free(field->data[i]);
    }
    free(field->data);
    free(field);
}

void set_value(Field* field, int i, int j, double value){
    field->data[i+field->nghost][j+field->nghost] = value;
}

double get_value(Field* field, int i, int j){
    return field->data[i+field->nghost][j+field->nghost];
}
