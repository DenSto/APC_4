#include <stdlib.h>
#include "assert.h"
#include "field.h"
#include "utils.h"
#include "grid.h"


Field *new_field(Grid* grid){
    int i,j;
    Field* newField = malloc(sizeof(Field));
    assert(grid != NULL && newField!=NULL);
    newField->grid = grid;
    newField->nghost = grid->nghost;
	newField->data = new_contiguous_2dArray(grid->nx + 2*grid->nghost,grid->ny + 2*grid->nghost);
    
    for(i = 0; i <grid->nx + 2*grid->nghost; i++){
        for(j = 0; j < grid->ny + 2*grid->nghost; j++){
            newField->data[i][j] = 0.0;
        }
    }
    return newField;
}

void free_field(Field* field){
	field->grid=NULL;
    free_contiguous_2dArray(field->data);
    free(field);
}
