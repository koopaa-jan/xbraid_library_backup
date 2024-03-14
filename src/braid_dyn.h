#ifndef braid_dyn_HEADER
#define braid_dyn_HEADER

#include "braid.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef braid_Real
(*braid_PtFcnGetValue)(braid_Vector u);


struct _braid_Core_struct_dyn;
/**
 * points to the core structure defined in _braid_dyn.h 
 **/
typedef struct _braid_Core_struct_dyn *braid_Core_dyn;

/**
 * Create a core object with the required initial data.
 *
 * This core is used by XBraid for dynamic internal data structures. 
 * The output is *core_ptr* which points to the newly created 
 * braid_Core structure. 
 **/
braid_Int
braid_Init_Dyn(braid_Real             tstart,      /**< start time */
               braid_Real             tstop,       /**< End time*/
               braid_Int              ntime,       /**< Initial number of temporal grid values*/
               braid_Real             interval_len,/**< length of lookup interval for change in processes*/
               braid_Int              max_procs,   /**< maximum number of processes to run in this simulation*/
               braid_App              app,         /**< User-defined _braid_App structure */
               braid_PtFcnStep        step,        /**< User time stepping routine to advance a braid_Vector forward one step */
               braid_PtFcnInit        init,        /**< Initialize a braid_Vector on the finest temporal grid*/
               braid_PtFcnClone       clone,       /**< Clone a braid_Vector*/
               braid_PtFcnFree        free,        /**< Free a braid_Vector*/
               braid_PtFcnSum         sum,         /**< Compute vector sum of two braid_Vectors*/
               braid_PtFcnSpatialNorm spatialnorm, /**< Compute norm of a braid_Vector, this is a norm only over space */
               braid_PtFcnGetValue    getValue,
               braid_PtFcnAccess      access,      /**< Allows access to XBraid and current braid_Vector */
               braid_PtFcnBufSize     bufsize,     /**< Computes size for MPI buffer for one braid_Vector */
               braid_PtFcnBufPack     bufpack,     /**< Packs MPI buffer to contain one braid_Vector*/
               braid_PtFcnBufUnpack   bufunpack,   /**< Unpacks MPI buffer into a braid_Vector */
               braid_Core_dyn         *core_ptr    /**< Pointer to braid_Core (_braid_Core) struct*/   
           );

braid_Int
braid_Destroy_Dyn(braid_Core_dyn core_dyn);

braid_Int
braid_Update_Dyn_Procs(braid_Int *iteration, MPI_Comm comm_world);

/*
 * Does the same as braid_Drive but is adjusted for dynamic processes
*/
braid_Int
braid_Drive_Dyn_Iterate(braid_Core  core,     /**< braid_Core (_braid_Core) struct*/
                        braid_Vector solVector
            );

/**
 * Carry out a simulation with XBraid, but dynamically newDyn. Integrate in time.
 **/
braid_Int
braid_Drive_Dyn(braid_Core_dyn  core_dyn            /**< braid_Core (_braid_Core) struct*/
            );

braid_Core braid_dyn_get_original_core(braid_Core_dyn core_dyn);

braid_Int
braid_Set_Info(MPI_Info info);

#ifdef __cplusplus
}
#endif

#endif

