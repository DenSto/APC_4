#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "assert.h"
#include "grid.h"
#include "field.h"

#define GHOST_CELLS 1
#define KAPPA 1.0

double dx=0;


int initialize_problem(int nx, double length, Grid** grid, Field** field, Field** rhs);
int compute_rhs(Grid* grid, Field* field, Field* rhs);
int compute_boundaries(Grid* grid, Field* field);
int time_step(double dt, Grid* grid, Field* field, Field* rhs);

double calculate_temperature(Grid* grid, Field* field);
int output_field(Grid* grid, Field* field);


int main(int argc, char* argv[]){
    
// Check to make sure usage is correct
    int nx;
    long step=0;

    if(argc != 2){
	printf("Usage: %s nx\n",argv[0]);
	return 1;
    }


    nx = atoi(argv[1]);
    
    dx = M_PI/nx;

    double dt = 0.5 * dx * dx / (4*KAPPA);
    double t = 0;
    double tF = 0.5 * M_PI * M_PI / KAPPA;
    printf("Total Time: %f     Time Step:%f Steps to take:%d\n",tF,dt, (int) (tF/dt));

    Grid* grid;
    Field *temp, *rhs;

    initialize_problem(nx,M_PI,&grid, &temp, &rhs);

    while(t <= tF){
        if(step % 1000 == 0) printf("Step: %ld\n",step);
        compute_rhs(grid,temp,rhs);

        time_step(dt,grid,temp,rhs);

        compute_boundaries(grid,temp);
        
        t += dt;
        step++;
    } 
    double ave_temp = calculate_temperature(grid,temp);
    printf("Final temperature:  %e\n", ave_temp);

    output_field(grid,temp);

    free_grid(grid);
    free_field(temp);
    free_field(rhs);


}

int initialize_problem(const int nx, const double length, Grid** grid, Field** temp, Field** rhs){
    double offset[2] = {0.0,0.0};

    *temp = new_field(nx,nx,GHOST_CELLS);
    *rhs = new_field(nx,nx,GHOST_CELLS);
    *grid = new_grid(nx,nx,length,length,GHOST_CELLS,offset);

    compute_boundaries(*grid, *temp);

    return 1;
}

int compute_boundaries(Grid* grid, Field* field){
    assert(field != NULL); 
    int nx = field->nx; 
    int ny = field->ny; 
    for(int i = 0; i < nx; i++){
        // Bottom cos(x)**2
        double x = grid_X(grid,i,-1);
        set_field_value(field,i,-1,cos(x)*cos(x));

        // Top sin(x)**2
        x = grid_X(grid,i,ny);
        set_field_value(field,i,ny,sin(x)*sin(x));

        // Left - Periodic

        set_field_value(field,-1,i,get_field_value(field,nx-1,i));
    
        // Right - Periodic
        set_field_value(field,nx,i,get_field_value(field,0,i));
    }
    return 1;
}

int compute_rhs(Grid* grid, Field* field, Field* rhs){
    assert(field != NULL && rhs != NULL); 
    for(int i = 0; i < field->nx; i++){
        for(int j = 0; j < field->ny; j++){
           double u,d,l,r,n,c,val;
           c= get_field_value(field,i,j);
           u= get_field_value(field,i,j+1);
           d= get_field_value(field,i,j-1);
           l= get_field_value(field,i-1,j);
           r= get_field_value(field,i+1,j);
           val = KAPPA * (u+d+l+r-4.0*c)/(dx*dx);
           set_field_value(rhs,i,j,val);
        } 
    }
    return 1;
}

int time_step(double dt, Grid* grid, Field* field, Field* rhs){
    assert(field != NULL && rhs != NULL); 
    for(int i = 0; i < field->nx; i++){
        for(int j = 0; j < field->ny; j++){
            double old = get_field_value(field,i,j);
            double new = dt*get_field_value(rhs,i,j);
            set_field_value(field,i,j,old+new);
        }
    } 
    return 1;
}

double calculate_temperature(Grid* grid, Field* field){
    double result=0;
    for(int i =0; i < field->nx; i++){
        for(int j =0; j < field->ny; j++){
            result += get_field_value(field,i,j);
        }
    }
    return result / (field->nx*field->ny);
}


int output_field(Grid* grid, Field* field){
    FILE *fp = fopen("out.dat","w");
    for(int i =0; i <= field->nx; i++){
        for(int j =0; j <= field->ny; j++){
            fprintf(fp,"%e %e %e\n",grid_X(grid,i,j),grid_Y(grid,i,j),get_field_value(field, i,j));
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    return 1;



}
