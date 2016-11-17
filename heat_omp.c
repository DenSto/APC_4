#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "grid.h"
#include "field.h"
#include <omp.h>

#define GHOST_CELLS 1
#define KAPPA 1.0
#define OUTPUT_CADENCE 10000


int nthreads=0;

int initialize_problem(int nx, double length, Grid** grid, Field** field, Field** rhs);
int compute_rhs(Grid* grid, Field* field, Field* rhs);
int compute_boundaries(Grid* grid, Field* field);
int time_step(double dt, Grid* grid, Field* field, Field* rhs);

double calculate_temperature(Grid* grid, Field* field);
int output_field(Grid* grid, Field* field);


int main(int argc, char* argv[]){
    
// Check to make sure usage is correct
    if(argc != 3){
	printf("Usage: %s nx nthreads\n",argv[0]);
	return 1;
    }

    int nx = atoi(argv[1]);
    nthreads=atoi(argv[2]);

    long step=0;

    double step_time=0.0;
    struct timeval tv_curr, tv_prev, tv_st, tv_fin;
    gettimeofday(&tv_prev,NULL);
    tv_st = tv_prev;

    omp_set_num_threads(nthreads);

    Grid* grid;
    Field *temp, *rhs;

    initialize_problem(nx,M_PI,&grid, &temp, &rhs);


    double dt = 0.5 * grid->dx * grid->dy / (4*KAPPA);
    double t = 0;
    double tF = 0.5 * M_PI * M_PI / KAPPA;
    printf("Total Time: %f     Time Step:%.2e Steps to take:%d\n",tF,dt, (int) (tF/dt));

    #pragma omp parallel shared(grid,temp,rhs,tF,tv_curr,tv_prev) firstprivate(step,step_time,t,dt)
    {
        while(t <= tF){
            #pragma omp master
            {
                if(step % OUTPUT_CADENCE == 0){
                    gettimeofday(&tv_curr,NULL);
                    step_time = (double)(tv_curr.tv_sec - tv_prev.tv_sec)+
                         1.0e-6*(double)(tv_curr.tv_usec - tv_prev.tv_usec);
                    tv_prev = tv_curr;
                    printf("Step: %ld   Step time: %.2e s \n",step,step_time);
                } 
            }

            compute_rhs(grid,temp,rhs);
            time_step(dt,grid,temp,rhs);

            #pragma omp barrier 

            compute_boundaries(grid,temp);
        
            t += dt;
            step++;
        }
    } 

    double ave_temp = calculate_temperature(grid,temp);
    output_field(grid,temp);


    gettimeofday(&tv_fin,NULL);
    double total_time = (double)(tv_fin.tv_sec - tv_st.tv_sec)+
                    1.0e-6*(double)(tv_fin.tv_usec - tv_st.tv_usec);

    printf("Final temperature:  %e   Total simulation time:%.2e s\n", ave_temp, total_time);

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
    int i;
    int nx = field->nx; 
    int ny = field->ny; 
   
    {
    #pragma omp for
        for(i = 0; i < nx; i++){
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
    }
    return 1;
}

int compute_rhs(Grid* grid, Field* field, Field* rhs){
    int nx = field->nx;
    int ny = field->ny;
    int i,j;
    double val;
    {
    #pragma omp for
    for(i = 0; i < nx; i++){
        for(j = 0; j < ny; j++){
           val= -4*get_field_value(field,i,j);
           val+= get_field_value(field,i,j+1);
           val+= get_field_value(field,i,j-1);
           val+= get_field_value(field,i-1,j);
           val+= get_field_value(field,i+1,j);
           val = KAPPA * val/(grid->dx*grid->dx);
           set_field_value(rhs,i,j,val);
        }
      }
    }
    return 1;
}

int time_step(double dt, Grid* grid, Field* field, Field* rhs){
    int i,j;
     {  
    #pragma omp for
    for(i = 0; i < field->nx; i++){
        for(j = 0; j < field->ny; j++){
            double old = get_field_value(field,i,j);
            double new = dt*get_field_value(rhs,i,j);
            set_field_value(field,i,j,old+new);
        }
      }
    } 
    return 1;
}

double calculate_temperature(Grid* grid, Field* field){
    double result=0;
    int i,j;
    #pragma omp parallel for shared(field,grid) private (i,j) reduction(+:result)
    for(i =0; i < field->nx; i++){
        for(j =0; j < field->ny; j++){
            result += get_field_value(field,i,j);
        }
    }
    return result / (field->nx*field->ny);
}


int output_field(Grid* grid, Field* field){
    int i,j;
    FILE *fp = fopen("out.dat","w");
    for(i =0; i <= field->nx; i++){
        for(j =0; j <= field->ny; j++){
            fprintf(fp,"%e %e %e\n",grid_X(grid,i,j),grid_Y(grid,i,j),get_field_value(field, i,j));
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    return 1;

