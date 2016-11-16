#ifndef FIELD_H
#define FIELD_H


typedef struct field_t Field;

struct field_t{
        double **data;
        int nx;
        int nghost;
};

Field *new_field(int nx, int ghost_cells);

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
