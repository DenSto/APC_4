#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "grid.h"
#include "field.h"

#ifdef OPENMP
#include <omp.h>
int nthreads=0;
#endif


#ifdef MPI
#include <mpi.h>
int  numtasks,rc;

#endif

#define GHOST_CELLS 1
#define KAPPA 1.0
#define OUTPUT_CADENCE 10000

int rank =0; 

int initialize_problem(Grid* grid, Field** field, Field** rhs);
int compute_rhs(Field* field, Field* rhs);
int compute_boundaries(Field* field);
int time_step(double dt, Field* field, Field* rhs);

double calculate_temperature(Grid* grid, Field* field);
int output_field(Field* field);


int main(int argc, char* argv[]){
    
// Check to make sure usage is correct
#ifdef OPENMP
    if(argc != 3){
	printf("Usage: %s nx nthreads\n",argv[0]);
	return 1;
    }
    nthreads=atoi(argv[2]);
    omp_set_num_threads(nthreads);
#else
    if(argc != 2){
	printf("Usage: %s nx\n",argv[0]);
	return 1;
    }
#endif

    int nx = atoi(argv[1]);
    int ny = nx;
    double length_x=M_PI, length_y=M_PI;
    double offset[2] = {0.0,0.0};

#ifdef MPI
    rc = MPI_Init(&argc,&argv);
    if(rc != MPI_SUCCESS){
        printf("Error Starting MPI Program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD,rc);
    } 

    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// Contruct information on the subgrid that constitutes our domain decomposition.
	// Decomposition is on the x direction. Can handle cases where nx is not divisible by nproc!
    int nx_global = nx;
    int chunk, mod;

    chunk = nx/numtasks;
    mod = nx % numtasks;

    nx=chunk;
    double off_n= (double)(rank*nx + mod);
    if(rank < mod){
        nx++;
        off_n= (double)(rank*(nx+1));
    }
    
    length_x  = M_PI*((double)nx/(double)nx_global);
    offset[0] = M_PI*(off_n/nx_global);

    Grid* grid = new_grid(nx_global,nx_global,M_PI,length_y,GHOST_CELLS,offset);
    Grid* subgrid = new_grid(nx,ny,length_x,length_y,GHOST_CELLS,offset);
    Field *temp, *rhs;
    initialize_problem(subgrid, &temp, &rhs);

#else

    Grid* grid = new_grid(nx,ny,length_x,length_y,GHOST_CELLS,offset);
    Field *temp, *rhs;
    initialize_problem(grid, &temp, &rhs);

#endif
    

	// Timestep information. CFL number is 0.5 here.
    double dt = 0.5 * grid->dx * grid->dy / (4*KAPPA);
    double t = 0;
    double tF = 0.5 * M_PI * M_PI / KAPPA;
    long step=0;

	// Timer for runtime and step cadences
    double step_time=0.0;
    struct timeval tv_curr, tv_prev, tv_st, tv_fin;
    gettimeofday(&tv_prev,NULL);
    tv_st = tv_prev;

    if(rank==0) {
        printf("Total Time: %f     Time Step:%.2e Steps to take:%d\n",tF,dt, (int) (tF/dt));
    }

	// Main time loop!
    #pragma omp parallel shared(grid,temp,rhs,tF,tv_curr,tv_prev) firstprivate(step,step_time,t,dt)
    {
        while(t <= tF){
            #pragma omp master
            {
                if(step % OUTPUT_CADENCE == 0 && rank == 0){
                    gettimeofday(&tv_curr,NULL);
                    step_time = (double)(tv_curr.tv_sec - tv_prev.tv_sec)+
                         1.0e-6*(double)(tv_curr.tv_usec - tv_prev.tv_usec);
                    tv_prev = tv_curr;
                    printf("Step: %ld   Step time: %.2e s \n",step,step_time);
                } 
            }

            compute_rhs(temp,rhs);
            time_step(dt,temp,rhs);

            #pragma omp barrier 

            compute_boundaries(temp);
        
            t += dt;
            step++;
        }
    } 

    double ave_temp = calculate_temperature(grid,temp);
   
	// Print-out temperature field to file (1 file per processor in MPI)
    output_field(temp);

	// Calculate runtime.
    gettimeofday(&tv_fin,NULL);
    double total_time = (double)(tv_fin.tv_sec - tv_st.tv_sec)+
                    1.0e-6*(double)(tv_fin.tv_usec - tv_st.tv_usec);

    if(rank == 0){
        printf("Final temperature:  %e   Total simulation time:%.2e s\n", ave_temp, total_time);
    }


	// Clean-up
    free_grid(grid);
    free_field(temp);
    free_field(rhs);

#ifdef MPI
    free_grid(subgrid);
    MPI_Finalize();
#endif 

}

// Construct fields
int initialize_problem(Grid* grid, Field** temp, Field** rhs){
    *temp = new_field(grid);
    *rhs  = new_field(grid);
    compute_boundaries(*temp);

    return 1;
}

// Compute boundary conditions (Periodic in x)
int compute_boundaries(Field* field){
	Grid* grid=field->grid;
    int i,j;
    int nx = grid->nx; 
    int ny = grid->ny; 
    int nghost=grid->nghost;
   
    {
#ifdef MPI

    int left= rank-1;
    int right=rank+1;
	int ierr1, ierr2, ierr3, ierr4;
    MPI_Status stat1,stat2;
    MPI_Request req1,req2;

    if(left<0) left=numtasks-1;
    if(right==numtasks) right=0;

    ierr1 = MPI_Isend(field->data[nghost],ny+2*nghost,MPI_DOUBLE,left,1,MPI_COMM_WORLD,&req1);
    ierr2 = MPI_Isend(field->data[nx],ny+2*nghost,MPI_DOUBLE,right,2,MPI_COMM_WORLD,&req2);
    ierr3 = MPI_Recv(field->data[0],ny+2*nghost,MPI_DOUBLE,left,2,MPI_COMM_WORLD,&stat1);
    ierr4 = MPI_Recv(field->data[nx+nghost],ny+2*nghost,MPI_DOUBLE,right,1,MPI_COMM_WORLD,&stat2);

#else
    #pragma omp for
        for(j = 0; j < ny; j++){
            // Left - Periodic
            set_field_value(field,-1,j,get_field_value(field,nx-1,j));
    
            // Right - Periodic
            set_field_value(field,nx,j,get_field_value(field,0,j));
        }
#endif
    #pragma omp for
        for(i = -nghost; i < nx+nghost; i++){
            // Bottom cos(x)**2
            double x = grid_X(field->grid,i,-1);
            set_field_value(field,i,-1,cos(x)*cos(x));

            // Top sin(x)**2
            x = grid_X(grid,i,ny);
            set_field_value(field,i,ny,sin(x)*sin(x));
        }
    }
    return 1;
}

// Compute right-hand side of field(temperature)-evolution equation
int compute_rhs(Field* field, Field* rhs){
    int nx = field->grid->nx;
    int ny = field->grid->ny;
	double dx = field->grid->dx;
	double dy = field->grid->dy;
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
           val = KAPPA * val/(dx*dy);
           set_field_value(rhs,i,j,val);
        }
      }
    }
    return 1;
}

// Time stepper. 1st order forward-Euler.
int time_step(double dt, Field* field, Field* rhs){
    int i,j;
	int nx = field->grid->nx;
	int ny = field->grid->ny;
     {  
    #pragma omp for
    for(i = 0; i < nx; i++){
        for(j = 0; j < ny; j++){
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
    int nx = field->grid->nx;
    int ny = field->grid->ny;
    #pragma omp parallel shared(field,grid) private (i,j) reduction(+:result)
	{
		#pragma omp for
    	for(i =0; i < nx; i++){
        	for(j =0; j < ny; j++){
            	result += get_field_value(field,i,j);
        	}
    	}
	}
#ifdef MPI
    double red_temp;
    MPI_Allreduce(&result,&red_temp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    result=red_temp;
#endif
    return result / (grid->nx*grid->ny);
}

int output_field(Field* field){
    int i,j;
	Grid* grid = field->grid;
    char fname[50] ="out.dat";
#ifdef MPI
    sprintf(fname,"out_%d.dat",rank);    
#endif 

    FILE *fp = fopen(fname,"w");
    for(i =0; i <= field->grid->nx; i++){
        for(j =0; j <= field->grid->ny; j++){
            fprintf(fp,"%e %e %e\n",grid_X(grid,i,j),grid_Y(grid,i,j),get_field_value(field, i,j));
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    return 1;
}
