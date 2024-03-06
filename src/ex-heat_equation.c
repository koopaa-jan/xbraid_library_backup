/**
 * Example:       ex-heat_equation.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-heat_equation
 * 
 * Sample run:    mpirun -np 2 --host n01,n02 ex-heat_equation 1.0 2
 * 
 **/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>


#include "braid_dyn.h"
#include "mpi.h"

/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

#define SIZE 50000 // number of temperature points, number of spatical points in the grid, attention > 100000 too big
#define NTIME 2000 // number of time steps to simulate attention, > 5000 crashed computer
#define ALPHA 1.0 // thermal diffusivity of the material, indicated how quickly the material conducts heat
#define DT 0.01 // time step size
#define DX 1.0 // distance between spacial points 
#define TSTOP 20.0 // end of time period, NTIME * DT

// CFL: alpha * DT / DX^2 <= 1/2 then the simulation is considered stable

/* App structure can contain anything, and be named anything as well */
typedef struct _braid_App_struct
{
    int       rank;
} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
    double value[SIZE]; // temperatures at one time step
} my_Vector;

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
    double tstart;             /* current time */
    double tstop;              /* evolve to this time*/
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

    double val_old_prior = u->value[0];

    for (int i = 1; i < SIZE - 1; i++) {
        double val = val_old_prior;
        val_old_prior = u->value[i];
        
        u->value[i] += (tstop - tstart) * ALPHA * (u->value[i+1] - 2 * u->value[i] + val) / (DX * DX);
    }

    return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
    my_Vector *u;

    u = (my_Vector *) malloc(sizeof(my_Vector));

    if (t == 0.0) {
        // example when the heat is only comming from the middle of the grid
        // for (int i = 0; i < SIZE; i++) {
        //     u->value[i] = 1.0; // Initial temperatures
        // }
        // Set initial condition
        // u->value[0] = 1.0; // Heat source on the left
        // u->value[SIZE / 2] = 1.0; // Heat source in the middle

        // sin wave between 0 and 1
        for (int i = 0; i < SIZE; i++) {
            u->value[i] = 0.5 + 0.5 * sin(2 * M_PI * i / SIZE); // Sine wave between 0 and 1
        }
    } else {
        for (int i = 0; i < SIZE; i++) {
            u->value[i] = 0.5; // random
        }
    }

    // printf("for t = %f values are: %f  %f  %f  %f  %f\n", t, u->value[0], u->value[1], u->value[2], u->value[3], u->value[4]);

    *u_ptr = u;

    return 0;
}

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
    my_Vector *v;

    v = (my_Vector *) malloc(sizeof(my_Vector));

    for (int i = 0; i < SIZE; i++) {
        (v->value[i]) = (u->value[i]);
    }
    
    *v_ptr = v;

    return 0;
}

int
my_Free(braid_App    app,
        braid_Vector u)
{
    free(u);
    return 0;
}

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{
    for (int i = 0; i < SIZE; i++) {
        (y->value[i]) = alpha*(x->value[i]) + beta*(y->value[i]);
    }
    return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
    double sum = 0.0;
    
    for (int i = 0; i < SIZE; i++) {
        sum += u->value[i] * u->value[i];
    }

    *norm_ptr = sqrt(sum);
    return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
    int        index;
    char       filename[255];
    FILE      *file;
    int        iteration = 0;

    braid_AccessStatusGetTIndex(astatus, &index);
    sprintf(filename, "%s.%04d.%04d.%03d", "ex-h_e.out", iteration, index, app->rank);

    file = fopen(filename, "r");

    while (file != NULL) {
        iteration++;
        fclose(file);

        sprintf(filename, "%s.%04d.%04d.%03d", "ex-h_e.out", iteration, index, app->rank);

        file = fopen(filename, "r");
    }
    if (!(iteration != 0 && index == 0)) {
        file = fopen(filename, "w");
        for (int i = 0; i < SIZE; i++) {
            fprintf(file, "%f ", u->value[i]);
        }
        fprintf(file, "\n");
        fflush(file);
        fclose(file);
    }

    return 0;
}

double
my_GetValue(braid_Vector u) {
    for (int i = 0; i < SIZE; i++) {
        printf("%f  ", u->value[i]);
    }
    printf("\n");
    return 0.;
}

int
my_BufSize(braid_App          app,
           int                *size_ptr,
           braid_BufferStatus bstatus)
{
    *size_ptr = sizeof(double) * SIZE;
    return 0;
}

int
my_BufPack(braid_App          app,
           braid_Vector       u,
           void               *buffer,
           braid_BufferStatus bstatus)
{
    double *dbuffer = buffer;

    memcpy(dbuffer, u->value, sizeof(double) * SIZE);
    braid_BufferStatusSetSize( bstatus, sizeof(double) * SIZE );

    return 0;
}

int
my_BufUnpack(braid_App          app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus bstatus)
{
    double    *dbuffer = buffer;
    my_Vector *u;

    u = (my_Vector *) malloc(sizeof(my_Vector));
    memcpy(u->value, dbuffer, sizeof(double) * SIZE);
    *u_ptr = u;

    return 0;
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{  
    braid_Core_dyn core_dyn;
    my_App       *app;
    double        tstart, tstop;
    int           ntime;

    /* Define time domain: ntime intervals */
    ntime  = NTIME;
    tstart = 0.0;
    tstop  = TSTOP;

    /* set up app structure */
    app = (my_App *) malloc(sizeof(my_App));
    // only used for my_access
    (app->rank)   = 0;

    //newDyn look for new parameter in command
    double interval_len = 0.0;
    char *endptr;
    // and define how many procs is the maximum
    int max_procs = 0;
    if (argc > 1) {
        interval_len = strtod(argv[1], &endptr);

        if (argc > 2) {
            max_procs = atoi(argv[2]);
        }
    }
    if (argc <= 1 || interval_len <= tstart || interval_len > tstop) {
        // make sure interval_len is valid
        interval_len = tstop;
    }

    //newDyn
    braid_Init_Dyn(tstart, tstop, ntime, interval_len, max_procs, app,
               my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, my_GetValue,
               my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core_dyn);

    /* Set some typical Braid parameters */
    braid_SetPrintLevel( braid_dyn_get_original_core(core_dyn), 2);
    braid_SetMaxLevels(braid_dyn_get_original_core(core_dyn), 2);
    braid_SetAbsTol(braid_dyn_get_original_core(core_dyn), 1.0e-06);

    //newDyn
    braid_Drive_Dyn(core_dyn);

    braid_Destroy_Dyn(core_dyn);
    free(app);

    for (int i = 0; i < argc; ++i) {
        printf("argv[%d]: %s\n", i, argv[i]);
    }

    return (0);
}