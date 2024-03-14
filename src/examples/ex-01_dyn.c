/**
 * Example:       ex-01.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-01
 *
 * Help with:     this is the simplest example available, read the source
 *
 * Sample run:    mpirun -np 2 ex-01
 *
 * Description:   solve the scalar ODE 
 *                   u' = lambda u, 
 *                   with lambda=-1 and y(0) = 1
 *                in a very simplified XBraid setting.
 *                
 *                When run with the default 10 time steps, the solution is:
 *                $ ./ex-01
 *                $ cat ex-01.out.00*
 *                  1.00000000000000e+00
 *                  6.66666666666667e-01
 *                  4.44444444444444e-01
 *                  2.96296296296296e-01
 *                  1.97530864197531e-01
 *                  1.31687242798354e-01
 *                  8.77914951989026e-02
 *                  5.85276634659351e-02
 *                  3.90184423106234e-02
 *                  2.60122948737489e-02
 *                  1.73415299158326e-02
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

/* App structure can contain anything, and be named anything as well */
typedef struct _braid_App_struct
{
    int       rank;
} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
    double value;
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

    /* Use backward Euler to propagate solution */
    (u->value) = 1./(1. + tstop-tstart)*(u->value);

    return 0;
}

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
    my_Vector *u;

    u = (my_Vector *) malloc(sizeof(my_Vector));
    if (t == 0.0) /* Initial condition */
    {
        (u->value) = 1.0;
    }
    else /* All other time points set to arbitrary value */
    {
        (u->value) = 0.456;
    }
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
    (v->value) = (u->value);
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
    (y->value) = alpha*(x->value) + beta*(y->value);
    return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
    double dot;

    dot = (u->value)*(u->value);
    *norm_ptr = sqrt(dot);
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
    sprintf(filename, "%s.%d.%04d.%03d", "ex-01_dyn.out", iteration, index, app->rank);

    file = fopen(filename, "r");

    while (file != NULL) {
        iteration++;
        fclose(file);

        sprintf(filename, "%s.%d.%04d.%03d", "ex-01_dyn.out", iteration, index, app->rank);

        file = fopen(filename, "r");
    }
    if (!(iteration != 0 && index == 0)) {
        file = fopen(filename, "w");
        fprintf(file, "%.14e\n", (u->value));
        fflush(file);
        fclose(file);
    }

    return 0;
}

double
my_GetValue(braid_Vector u) {
    printf("value of u: %f\n", u->value);
    return 0.;
}

int
my_BufSize(braid_App          app,
           int                *size_ptr,
           braid_BufferStatus bstatus)
{
    *size_ptr = sizeof(double);
    return 0;
}

int
my_BufPack(braid_App          app,
           braid_Vector       u,
           void               *buffer,
           braid_BufferStatus bstatus)
{
    double *dbuffer = buffer;

    dbuffer[0] = (u->value);
    braid_BufferStatusSetSize( bstatus, sizeof(double) );

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
    (u->value) = dbuffer[0];
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
    //newDyn
    int           ntime;//, rank, size;

    /* Define time domain: ntime intervals */
    ntime  = 10;
    tstart = 0.0;
    tstop  = tstart + ntime/2.;

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
    braid_SetPrintLevel(braid_dyn_get_original_core(core_dyn), 2);
    braid_SetMaxLevels(braid_dyn_get_original_core(core_dyn), 2);
    braid_SetAbsTol(braid_dyn_get_original_core(core_dyn), 1.0e-06);
    braid_SetCFactor(braid_dyn_get_original_core(core_dyn), -1, 2);

    /* Set the info for reconfigurations */
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "mpi_num_procs_add", "2");
    braid_Set_Info(info);
    MPI_Info_free(&info);

    //newDyn
    braid_Drive_Dyn(core_dyn);

    braid_Destroy_Dyn(core_dyn);
    free(app);

    for (int i = 0; i < argc; ++i) {
        printf("argv[%d]: %s\n", i, argv[i]);
    }

    return (0);
}