/*
 * Field : Structure of field quantities (Temperature, density, etc...)
 *         Requires a grid object to construct.
 */
#ifndef FIELD_H
#define FIELD_H

#include "grid.h"

typedef struct field_t Field;

struct field_t{
        double **data;
		int nghost; // kept for speed!
		Grid* grid;
};

Field *new_field(Grid* grid);

void free_field(Field* field);
void set_value(Field* field, int i, int j, double value);
double get_value(Field* field, int i, int j);

static inline void set_field_value(Field* field, int i, int j, double value){
    field->data[i+field->nghost][j+field->nghost] = value;
}

static inline double get_field_value(Field* field, int i, int j){
    return field->data[i+field->nghost][j+field->nghost];
}

#endif
