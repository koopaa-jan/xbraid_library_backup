#include "_braid_dyn.h"
#include "util.h"
#include "dmr.h"
#include <stdbool.h>

#ifndef DEBUG_DYN
#define DEBUG_DYN 0
#endif

// counts the iteration
static braid_Int iteration = 0;

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_Update_Dyn_Procs(braid_Int *iter, MPI_Comm comm_world)
{
   // sending counter of iteration from first process to all other including the dynamically added process
   // so they get the knowledge of how far the algorithm is so they can calculate their parameters correctly
   MPI_Bcast(iter, 1, MPI_INT, 0, comm_world);
   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_Drive_Dyn_Iterate(braid_Core  core, braid_Vector transfer_vector)
{
   MPI_Comm             comm_world      = _braid_CoreElt(core, comm_world);
   braid_Int            myid            = _braid_CoreElt(core, myid_world);
   braid_Real           tstart          = _braid_CoreElt(core, tstart);
   braid_Real           tstop           = _braid_CoreElt(core, tstop);
   braid_Int            ntime           = _braid_CoreElt(core, ntime);
   braid_Int            warm_restart    = _braid_CoreElt(core, warm_restart);
   braid_Int            print_level     = _braid_CoreElt(core, print_level);
   braid_App            app             = _braid_CoreElt(core, app);
   braid_Int            obj_only        = _braid_CoreElt(core, obj_only);
   braid_Int            adjoint         = _braid_CoreElt(core, adjoint);
   braid_Int            delta           = _braid_CoreElt(core, delta_correct);
   braid_SyncStatus     sstatus         = (braid_SyncStatus)core;

   braid_Int      ilower, iupper, i;
   braid_Real    *ta;
   _braid_Grid   *grid;
   braid_Real     localtime, globaltime;
   braid_Real     timer_drive_init;

   timer_drive_init = _braid_MPI_Wtime(core, 2);
   /* Check for non-supported adjoint features */
   if (adjoint)
   {
      _braid_AdjointFeatureCheck(core);
   }
   /* Check for non-supported Delta correction features */
   if (delta)
   {
      _braid_DeltaFeatureCheck(core);
   }
   if (myid == 0 )
   {
      if (!warm_restart && print_level > 0)
      {
         _braid_printf("\n  Braid: Begin simulation, %d time steps\n",
                       _braid_CoreElt(core, gupper));
      }
      if ( adjoint && print_level > 0 )
      {
         if (_braid_CoreElt(core, max_levels) > 1)
         {
            _braid_printf("\n");
            _braid_printf("  Braid:      || r ||      || r_adj ||     Objective\n");
            _braid_printf("  Braid:---------------------------------------------\n");
         }
         else
         {
            _braid_printf("  Braid: Serial time-stepping. \n\n");
         }
      }
   }

   /* Start timer */
   localtime = _braid_MPI_Wtime(core, 1);

   /* Allocate and initialize grids */
   if ( !warm_restart )
   {
      /* Create fine grid */
      _braid_GetDistribution(core, &ilower, &iupper);

      _braid_GridInit(core, 0, ilower, iupper, &grid);

      /* Set t values */
      ta = _braid_GridElt(grid, ta);
      if ( _braid_CoreElt(core, tgrid) != NULL )
      {
         /* Call the user's time grid routine */
         _braid_BaseTimeGrid(core, app, ta, &ilower, &iupper);
      }
      else
      {
         for (i = ilower; i <= iupper; i++)
         {
            ta[i-ilower] = tstart + (((braid_Real)i)/ntime)*(tstop-tstart);
            // printf("myid: %d and ta[%d]: %f\n", myid, i-ilower, ta[i-ilower]);
         }
      }

      /* Create a grid hierarchy */
      _braid_InitHierarchy(core, grid, 0);

      /* Set initial values */
      _braid_InitGuess(core, 0);

      //setting solution vector for next iteration
      if (transfer_vector != NULL && myid == 0) {
         _braid_USetVector_Dyn(core, 0, 0, transfer_vector, 0);
      }

      /* Let the users sync after initialization */
      _braid_SyncStatusInit(0, // Iteration
                            0, // Level
                            _braid_CoreElt(core, nrefine),
                            _braid_CoreElt(core, gupper),
                            0, // done
                            braid_ASCaller_Drive_AfterInit,
                            sstatus);
      _braid_Sync(core, sstatus);
   }

   /* Initialize sensitivity computation */
   if ( adjoint )
   {
      if (!warm_restart)
      {
         /* Initialize and allocate the adjoint variables */
         _braid_InitAdjointVars(core, grid);
      }
      else
      {
         /* Prepare for next adjoint iteration in case of warm_restart */
         _braid_CoreElt(core, optim)->sum_user_obj  = 0.0;
         _braid_CoreElt(core, optim)->f_bar         = 0.0;
         if (!obj_only)
         {
            _braid_CoreFcn(core, reset_gradient)(_braid_CoreElt(core, app));
         }
      }

      if ( obj_only )
      {
         _braid_CoreElt(core, record) = 0;
      }
      else
      {
         _braid_CoreElt(core, record) = 1;
      }
   }
   _braid_CoreElt(core, timer_drive_init) += _braid_MPI_Wtime(core, 2) - timer_drive_init;

   /* Reset from previous calls to braid_drive() */
   _braid_CoreElt(core, done) = 0;

   /* Solve with MGRIT */
   _braid_Drive(core, localtime);

   /* Stop timer */
   localtime = _braid_MPI_Wtime(core, 1) - localtime;
   MPI_Allreduce(&localtime, &globaltime, 1, braid_MPI_REAL, MPI_MAX, comm_world);
   _braid_CoreElt(core, localtime)  = localtime;
   _braid_CoreElt(core, globaltime) = globaltime;

   /* Print statistics for this run */
   if ( (print_level > 1) && (myid == 0) )
   {
      braid_PrintStats(core);
   }

   /* Print more intrusive timing information, only if level is >= 2 */
   if (_braid_CoreElt(core, timings) >= 2)
   {
      braid_PrintTimers(core);
   }

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_Drive_Dyn(braid_Core_dyn  core_dyn)
{
   braid_Core core = _braid_CoreDynElt(core_dyn, original_core);

   MPI_Comm             comm_world      = _braid_CoreElt(core, comm_world);
   braid_Int            myid            = _braid_CoreElt(core, myid_world);
   braid_Real           tstart          = _braid_CoreElt(core, tstart);
   braid_Real           tstop           = _braid_CoreElt(core, tstop);
   braid_Int            ntime           = _braid_CoreElt(core, ntime);
   //newDyn
   braid_Real           interval_len    = _braid_CoreDynElt(core_dyn, interval_len);
   // braid_Int            max_procs       = _braid_CoreDynElt(core_dyn, max_procs);

   braid_App            app             = _braid_CoreElt(core, app);

   braid_Real globaltstart = tstart;
   braid_Real globaltstop = tstop;
   braid_Real globalinterval = globaltstop - globaltstart;

   // size of current communicator
   braid_Int size;


   //making sure interval_len is a valid timestemp
   braid_Real trange_per_ts = globalinterval / ntime;
   
   if (interval_len <= globaltstart - tstart || interval_len > globalinterval) {
      interval_len = globalinterval;
   } else {
      braid_Int rem = ceil(interval_len / trange_per_ts);

      interval_len = trange_per_ts * rem + tstart;

      if (interval_len == tstart) {
         interval_len = trange_per_ts + tstart;
      }
   }

   // printf("-------------------------------------------------interval_len: %f trange_per_ts: %f\n", interval_len, trange_per_ts);

   //helps to get the solution vector and store for next iteration
   braid_BaseVector transfer_vector = _braid_TAlloc(_braid_BaseVector, 1);
   transfer_vector->userVector = NULL;
   transfer_vector->bar        = NULL;
   transfer_vector->basis      = NULL;

   // counts at which time step the program currently is
   braid_Real current_ts = globaltstart;

   // gupper is global size of fine grid
   _braid_CoreElt(core, gupper) = (interval_len / trange_per_ts);
   current_ts = globaltstart + (interval_len * iteration);


   // repeat as long as the next iteration(which has the length of interval_len) wouldnt exceed tstop
   for (; current_ts < globaltstop - interval_len; current_ts += interval_len) {
      // sleep(1);
      ++iteration;
      //set parameters
      _braid_CoreElt(core, tstart) = current_ts;
      _braid_CoreElt(core, tstop) = current_ts + interval_len;
      _braid_CoreElt(core, ntime) = interval_len / trange_per_ts;

      // printf("1 ++++++++++++++ tstart: %f tstop: %f ntime: %f gupper: %f\n",
      //  current_ts, current_ts + interval_len, interval_len / trange_per_ts,
      //   (interval_len / trange_per_ts));


      braid_Drive_Dyn_Iterate(core, transfer_vector->userVector);

      MPI_Barrier(comm_world);

      //get solution Vector of previous run, store in transfer_vector and set as start vector in braid_Drive_Dyn_Iterate
      //get last vector of the process which is responsible for it, which is the process with the highest id
      MPI_Comm_size(comm_world, &size);


      // starting time measurement for duration of the solution vector distribution
      braid_Real gtime_sol_vec;
      braid_Real ltime_sol_vec = _braid_MPI_Wtime(core, 1);



      // adjusting sol_vec distribution so not the last one but the last one with a valid sol_vec is distributing
      // if number of procs available is higher than number of ts to be computed, then the vector with the sol vector is
      // not the last process
      braid_Int num_procs_needed = (interval_len / trange_per_ts) + 1;
      braid_Int sol_vec_id = size - 1;
      if (num_procs_needed < size) {
         sol_vec_id = (interval_len / trange_per_ts);
      }
      
      // allocating parameters for transfering solution vector
      // and getting size of solution vector
      braid_BufferStatus bstatus = (braid_BufferStatus)core;
      braid_Int sol_vec_size;

      char *buffer;

      if (((myid == sol_vec_id) || (myid == 0)) && size > 1) {
         _braid_CoreFcn(core, bufsize)(app, &sol_vec_size, bstatus);

         buffer = (char *)malloc(sol_vec_size);
      }

      if (myid == sol_vec_id) {
         if (myid == size - 1) {
            _braid_UGetLast(core, &transfer_vector);
         } else {
            // get sol vector from processes with the last time step that is not the last process, so UGetLast cant be used
            _braid_Grid       **grids    = _braid_CoreElt(core, grids);
            braid_Int           cfactor  = _braid_GridElt(grids[0], cfactor);

            braid_BaseVector bv;
            if ( (_braid_CoreElt(core, storage) < 0) && !(_braid_IsCPoint(gupper, cfactor)) ) { 
               bv = _braid_GridElt(grids[0], ulast);
               transfer_vector->userVector = bv->userVector;
            }
            else {
               _braid_UGetVectorRef(core, 0, gupper, &bv);
               transfer_vector->userVector = bv->userVector;
            }
         }

         // only send when there are more than one processes
         if (myid != 0) {
            // serializing data of solution vector
            _braid_CoreFcn(core, bufpack)(app, transfer_vector->userVector, buffer, bstatus);
            MPI_Send(buffer, sol_vec_size, MPI_BYTE, 0, 9, comm_world);
         }
      }
      
      if (size > 1 && myid == 0) {
         MPI_Recv(buffer, sol_vec_size, MPI_BYTE, sol_vec_id, 9, comm_world, MPI_STATUS_IGNORE);
         // deserializing data of solution vector
         _braid_CoreFcn(core, bufunpack)(app, buffer, &transfer_vector->userVector, bstatus);
      }

      if ((myid == sol_vec_id || myid == 0) && size > 1) {
         free(buffer);
      }

      // ending time measure for distribution of solution vector
      ltime_sol_vec = _braid_MPI_Wtime(core, 1) - ltime_sol_vec;
      MPI_Allreduce(&ltime_sol_vec, &gtime_sol_vec, 1, braid_MPI_REAL, MPI_MAX, comm_world);

      if (myid == 0) {
         printf("time for distribution of solution vector was: %f with %d processes\n", gtime_sol_vec, size);
      }
      

      // starting time measurement for duration of processes change and update
      braid_Real ltime_procs = _braid_MPI_Wtime(core, 1);
      
      // changing amount of processes
      // braid_Int num_procs_sub = 0;
      // braid_Int num_procs_add = 2;

      // if adding or removing isnt possible, continue without adjusting
      //if (((num_procs_add == 0) && (size - num_procs_sub > 0)) || ((size + num_procs_add <= max_procs) && (num_procs_sub == 0))) {
         //printf("-+-+-+--+--+-+-++--+---++- myid is: %d and old size: %d +-+-+-+-+-+-+-+-+\n", myid, size);

      DMR_RECONFIGURATION(
            braid_Update_Dyn_Procs(&iteration, DMR_INTERCOMM),
            NULL, 
            NULL, 
            NULL);

      //update new processes
      //braid_Update_Dyn_Procs(&iteration, DMR_INTERCOMM);

      // updating parameters
      comm_world = DMR_INTERCOMM;
      myid = DMR_comm_rank;

      _braid_CoreElt(core, comm_world)      = comm_world;
      _braid_CoreElt(core, comm)            = comm_world;
      _braid_CoreElt(core, myid_world)      = myid;
      _braid_CoreElt(core, myid)            = myid;
      //} else {
      //   printf("reconfiguration was skipped as removing or adding processes wasnt possible!\n");
      //}


      // ending time measure for duration of processes change and update
      ltime_procs = _braid_MPI_Wtime(core, 1) - ltime_procs;
      
      if (myid == 0) {
         printf("time for duration of processes change and update was: %f\n", ltime_procs);
      }
      
      //MPI_Comm_size(comm_world, &size);
      //printf("-+-+-+--+--+-+-++--+---++- my new id is: %d new size after rearrange is: %d ++-+-+-+-+--+-+-+-+\n", myid, size);
   }

   if (DMR_FINALIZE_FLAG) {
      return _braid_error_flag;
   }

   if (current_ts < globaltstop) {
      // compute the remainding time until tstop

      _braid_CoreElt(core, tstart) = current_ts;
      _braid_CoreElt(core, tstop) = globaltstop;

      _braid_CoreElt(core, ntime) = (globaltstop - current_ts) / trange_per_ts;

      _braid_CoreElt(core, gupper) = ((globaltstop - current_ts) / trange_per_ts);


      // printf("2 ++++++++++++++ tstart: %f tstop: %f ntime: %f gupper: %f myid: %d\n",
      //  current_ts, globaltstop, (globaltstop - current_ts) / trange_per_ts,
      //   ((globaltstop - current_ts) / trange_per_ts), myid);

      braid_Drive_Dyn_Iterate(core, transfer_vector->userVector);

   }

   free(transfer_vector);

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_Destroy_Dyn(braid_Core_dyn core_dyn) {

   braid_Core core = _braid_CoreDynElt(core_dyn, original_core);

   if (core)
   {
      braid_App               app             = _braid_CoreElt(core, app);
      braid_Int               nlevels         = _braid_CoreElt(core, nlevels);
      _braid_Grid           **grids           = _braid_CoreElt(core, grids);
      braid_Int               cfactor         = _braid_GridElt(grids[0], cfactor);
      braid_Int               gupper          = _braid_CoreElt(core, gupper);
      braid_Int               richardson      = _braid_CoreElt(core, richardson);
      braid_Int               est_error       = _braid_CoreElt(core, est_error);
      char                   *timer_file_stem = _braid_CoreElt(core, timer_file_stem);
      braid_Int               level;

      _braid_TFree(_braid_CoreElt(core, nrels));
      _braid_TFree(_braid_CoreElt(core, CWts));
      _braid_TFree(_braid_CoreElt(core, rnorms));
      _braid_TFree(_braid_CoreElt(core, full_rnorms));
      _braid_TFree(_braid_CoreElt(core, cfactors));
      _braid_TFree(_braid_CoreElt(core, rfactors));
      _braid_TFree(_braid_CoreElt(core, tnorm_a));
      _braid_TFree(_braid_CoreElt(core, rdtvalues));

      /* Destroy stored Lyapunov exponents */
      if (_braid_CoreElt(core, lyap_exp))
      {

         braid_Int npoints;
         if (nlevels == 1)
         {
            npoints = _braid_GridElt(grids[0], iupper) - _braid_GridElt(grids[0], ilower);
         }
         else
         {
            npoints = _braid_GridElt(grids[0], ncpoints);
         }
         for (braid_Int i = 0; i < npoints; i++)
         {
            _braid_TFree(_braid_CoreElt(core, local_exponents)[i]);
         }
         _braid_TFree(_braid_CoreElt(core, local_exponents));
      }

      /* Destroy the optimization structure */
      _braid_CoreElt(core, record) = 0;
      if (_braid_CoreElt(core, adjoint))
      {
         _braid_OptimDestroy( core );
         _braid_TFree(_braid_CoreElt(core, optim));
      }

      /* Free last time step, if set */
      if ( (_braid_CoreElt(core, storage) < 0) && !(_braid_IsCPoint(gupper, cfactor)) )
      {
         if (_braid_GridElt(grids[0], ulast) != NULL)
         {
           _braid_BaseFree(core, app, _braid_GridElt(grids[0], ulast));
         }
      }

      /* Destroy Richardson estimate structures */
      if ( richardson || est_error )
      {
         _braid_TFree(_braid_CoreElt(core, dtk));
         if (est_error)
         {
            _braid_TFree(_braid_CoreElt(core, estimate));
         }
      }

      for (level = 0; level < nlevels; level++)
      {
         _braid_GridDestroy(core, grids[level]);
      }

      _braid_TFree(grids);

      _braid_TFree(core);
      _braid_TFree(core_dyn);

      if (timer_file_stem != NULL)
      {
         _braid_TFree(timer_file_stem);
      }

   }

   if (_braid_printfile != NULL)
   {
      fclose(_braid_printfile);
   }

   //newDyn
   DMR_FINALIZE0();

   return _braid_error_flag;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

braid_Int
braid_Init_Dyn(
           braid_Real             tstart,
           braid_Real             tstop,
           braid_Int              ntime,
           //newDyn
           braid_Real             interval_len,
           braid_Int              max_procs,
           braid_App              app,
           braid_PtFcnStep        step,
           braid_PtFcnInit        init,
           braid_PtFcnClone       clone,
           braid_PtFcnFree        free,
           braid_PtFcnSum         sum,
           braid_PtFcnSpatialNorm spatialnorm,
           //newDyn
           braid_PtFcnGetValue    getValue,
           braid_PtFcnAccess      access,
           braid_PtFcnBufSize     bufsize,
           braid_PtFcnBufPack     bufpack,
           braid_PtFcnBufUnpack   bufunpack,
           braid_Core_dyn         *core_ptr)
{
   _braid_Core_dyn           *core_dyn;
   _braid_Core               *original_core;

   /* Braid default values */
   braid_Int              cfdefault        = 2;              /* Default coarsening factor */
   braid_Int              nrdefault        = 1;              /* Default number of FC sweeps on each level */
   braid_Real             CWt_default      = 1.0;            /* Default C-relaxtion weight  */
   braid_Int              fmg              = 0;              /* Default fmg (0 is off) */
   braid_Int              nfmg             = -1;             /* Default fmg cycles is -1, indicating all fmg-cycles (if fmg=1) */
   braid_Int              nfmg_Vcyc        = 1;              /* Default num V-cycles at each fmg level is 1 */
   braid_Int              max_iter         = 100;            /* Default max_iter */
   braid_Int              max_levels       = 30;             /* Default max_levels */
   braid_Int              incr_max_levels  = 0;              /* Default increment max levels is false */
   braid_Int              min_coarse       = 2;              /* Default min_coarse */
   braid_Int              relax_only_cg    = 0;              /* Default to no relaxation on coarsest grid (alternative to serial solve) */
   braid_Int              seq_soln         = 0;              /* Default initial guess is from user's Init() function */
   braid_Int              print_level      = 2;              /* Default print level */
   braid_Int              io_level         = 1;              /* Default output-to-file level */
   braid_Int              access_level     = 1;              /* Default access level */
   braid_Int              finalFCrelax     = 0;              /* Default final FCrelax */
   braid_Int              tnorm            = 2;              /* Default temporal norm */
   braid_Real             tol              = 1.0e-09;        /* Default absolute tolerance */
   braid_Int              warm_restart     = 0;              /* Default is no warm restart */
   braid_Int              rtol             = 1;              /* Use relative tolerance */
   braid_Int              skip             = 1;              /* Default skip value, skips all work on first down-cycle */
   braid_Int              max_refinements  = 200;            /* Maximum number of F-refinements */
   braid_Int              tpoints_cutoff   = braid_Int_Max;  /* Maximum number of time steps, controls FRefine()*/
   braid_Int              delta_correct    = 0;              /* Default Delta correction: Turned off */
   braid_Int              delta_rank       = 0;              /* Default Delta correction rank (always set by setDelta())*/
   braid_Int              estimate_lyap    = 0;              /* Default estimation of Lyapunov vectors is turned off*/
   braid_Int              relax_lyap       = 0;              /* Default propagation of Lyapunov vectors during FCRelax is turned off*/
   braid_Int              lyap_exp         = 0;              /* Default Lyapunov exponents are not saved */
   braid_Int              delta_defer_lvl  = 0;              /* Default Delta correction defer level (set by setDeltaDefer())*/
   braid_Int              delta_defer_iter = 0;              /* Default Delta correction defer iteration (set by setDeltaDefer())*/
   braid_Int              adjoint          = 0;              /* Default adjoint run: Turned off */
   braid_Int              record           = 0;              /* Default action recording: Turned off */
   braid_Int              obj_only         = 0;              /* Default objective only: Turned off */
   braid_Int              reverted_ranks   = 0;              /* Default objective only: Turned off */
   braid_Int              verbose_adj      = 0;              /* Default adjoint verbosity Turned off */

   braid_Int              myid_world,  myid, world_size;

   MPI_Comm comm_world;
   MPI_Comm comm;

   // newDyn
   DMR_INIT(-1, NULL, braid_Update_Dyn_Procs(&iteration, DMR_INTERCOMM));

   comm_world = DMR_INTERCOMM;
   comm = DMR_INTERCOMM;

   myid_world = DMR_comm_rank;
   myid = DMR_comm_rank;
   world_size = DMR_comm_size;

   // printf("Size:%d and rank:%d\n", world_size, myid_world);

   core_dyn = _braid_CTAlloc(_braid_Core_dyn, 1);
   original_core = _braid_CTAlloc(_braid_Core, 1);

   //newDyn
   _braid_CoreDynElt(core_dyn, interval_len)    = interval_len;
   _braid_CoreDynElt(core_dyn, getValue)        = getValue;
   _braid_CoreDynElt(core_dyn, max_procs)       = max_procs;
   _braid_CoreDynElt(core_dyn, original_core)   = original_core;

   _braid_Core *core = _braid_CoreDynElt(core_dyn, original_core);

   _braid_CoreElt(core, comm_world)      = comm_world;

   _braid_CoreElt(core, comm)            = comm;
   _braid_CoreElt(core, myid_world)      = myid_world;
   _braid_CoreElt(core, myid)            = myid;
   _braid_CoreElt(core, tstart)          = tstart;
   _braid_CoreElt(core, tstop)           = tstop;
   _braid_CoreElt(core, ntime)           = ntime;
   _braid_CoreElt(core, done)            = 0;
   _braid_CoreElt(core, app)             = app;

   _braid_CoreElt(core, step)            = step;
   _braid_CoreElt(core, init)            = init;
   _braid_CoreElt(core, sinit)           = NULL;
   _braid_CoreElt(core, init_basis)      = NULL;
   _braid_CoreElt(core, clone)           = clone;
   _braid_CoreElt(core, sclone)          = NULL;
   _braid_CoreElt(core, free)            = free;
   _braid_CoreElt(core, sfree)           = NULL;
   _braid_CoreElt(core, sum)             = sum;
   _braid_CoreElt(core, spatialnorm)     = spatialnorm;
   _braid_CoreElt(core, inner_prod)      = NULL;
   _braid_CoreElt(core, access)          = access;
   _braid_CoreElt(core, bufsize)         = bufsize;
   _braid_CoreElt(core, bufpack)         = bufpack;
   _braid_CoreElt(core, bufunpack)       = bufunpack;
   _braid_CoreElt(core, residual)        = NULL;
   _braid_CoreElt(core, scoarsen)        = NULL;
   _braid_CoreElt(core, srefine)         = NULL;
   _braid_CoreElt(core, tgrid)           = NULL;
   _braid_CoreElt(core, sync)            = NULL;
   _braid_CoreElt(core, bufalloc)        = NULL;
   _braid_CoreElt(core, buffree)         = NULL;

   _braid_CoreElt(core, access_level)    = access_level;
   _braid_CoreElt(core, finalFCrelax)    = finalFCrelax;
   _braid_CoreElt(core, tnorm)           = tnorm;
   _braid_CoreElt(core, print_level)     = print_level;
   _braid_CoreElt(core, io_level)        = io_level;
   _braid_CoreElt(core, max_levels)      = 0; /* Set with SetMaxLevels() below */
   _braid_CoreElt(core, incr_max_levels) = incr_max_levels;
   _braid_CoreElt(core, min_coarse)      = min_coarse;
   _braid_CoreElt(core, relax_only_cg)   = relax_only_cg;
   _braid_CoreElt(core, seq_soln)        = seq_soln;
   _braid_CoreElt(core, tol)             = tol;
   _braid_CoreElt(core, rtol)            = rtol;
   _braid_CoreElt(core, warm_restart)    = warm_restart;

   _braid_CoreElt(core, nrels)           = NULL; /* Set with SetMaxLevels() below */
   _braid_CoreElt(core, nrdefault)       = nrdefault;
   _braid_CoreElt(core, CWts)            = NULL; /* Set with SetMaxLevels() below */
   _braid_CoreElt(core, CWt_default)     = CWt_default;

   _braid_CoreElt(core, cfactors)        = NULL; /* Set with SetMaxLevels() below */
   _braid_CoreElt(core, cfdefault)       = cfdefault;

   _braid_CoreElt(core, max_iter)        = 0; /* Set with SetMaxIter() below */
   _braid_CoreElt(core, niter)           = 0;
   _braid_CoreElt(core, fmg)             = fmg;
   _braid_CoreElt(core, nfmg)            = nfmg;
   _braid_CoreElt(core, nfmg_Vcyc)       = nfmg_Vcyc;

   _braid_CoreElt(core, storage)         = -1;            /* only store C-points */
   _braid_CoreElt(core, useshell)         = 0;

   _braid_CoreElt(core, gupper)          = 0; /* Set with SetPeriodic() below */

   _braid_CoreElt(core, refine)          = 0;  /* Time refinement off by default */
   _braid_CoreElt(core, rfactors)        = NULL;
   _braid_CoreElt(core, rdtvalues)       = NULL;
   _braid_CoreElt(core, r_space)         = 0;
   _braid_CoreElt(core, rstopped)        = -1;
   _braid_CoreElt(core, nrefine)         = 0;
   _braid_CoreElt(core, max_refinements) = max_refinements;
   _braid_CoreElt(core, tpoints_cutoff)  = tpoints_cutoff;

   _braid_CoreElt(core, nlevels)         = 0;
   _braid_CoreElt(core, grids)           = NULL; /* Set with SetMaxLevels() below */

   _braid_CoreElt(core, skip)            = skip;

   /* Delta correction */
   _braid_CoreElt(core, delta_correct)    = delta_correct;
   _braid_CoreElt(core, delta_rank)       = delta_rank;
   _braid_CoreElt(core, estimate_lyap)    = estimate_lyap;
   _braid_CoreElt(core, relax_lyap)       = relax_lyap;
   _braid_CoreElt(core, lyap_exp)         = lyap_exp;
   _braid_CoreElt(core, delta_defer_iter) = delta_defer_iter;
   _braid_CoreElt(core, delta_defer_lvl)  = delta_defer_lvl;


   /* Adjoint */
   _braid_CoreElt(core, adjoint)               = adjoint;
   _braid_CoreElt(core, record)                = record;
   _braid_CoreElt(core, obj_only)              = obj_only;
   _braid_CoreElt(core, reverted_ranks)        = reverted_ranks;
   _braid_CoreElt(core, verbose_adj)           = verbose_adj;
   _braid_CoreElt(core, actionTape)            = NULL;
   _braid_CoreElt(core, userVectorTape)        = NULL;
   _braid_CoreElt(core, barTape)               = NULL;
   _braid_CoreElt(core, optim)                 = NULL;
   _braid_CoreElt(core, objectiveT)            = NULL;
   _braid_CoreElt(core, objT_diff)             = NULL;
   _braid_CoreElt(core, step_diff)             = NULL;
   _braid_CoreElt(core, reset_gradient)        = NULL;
   _braid_CoreElt(core, postprocess_obj)       = NULL;
   _braid_CoreElt(core, postprocess_obj_diff)  = NULL;

   /* Residual history and accuracy tracking for StepStatus*/
   _braid_CoreElt(core, rnorm0)              = braid_INVALID_RNORM;
   _braid_CoreElt(core, rnorms)              = NULL; /* Set with SetMaxIter() below */
   _braid_CoreElt(core, full_rnorm_res)      = NULL;
   _braid_CoreElt(core, full_rnorm0)         = braid_INVALID_RNORM;
   _braid_CoreElt(core, full_rnorms)         = NULL; /* Set with SetMaxIter() below */
   _braid_CoreElt(core, old_fine_tolx)       = -1.0;
   _braid_CoreElt(core, tight_fine_tolx)     = 1;

   /* Richardson Error Estimation */
   _braid_CoreElt(core, est_error)       = 0;     /* Error Estimation off by default */
   _braid_CoreElt(core, richardson)      = 0;     /* Richardson Extrapolation off by default */
   _braid_CoreElt(core, order)           = 2;     /* 2nd order local time stepper by defualt */
   _braid_CoreElt(core, dtk)             = NULL;  /* Set in _braid_InitHierarchy */
   _braid_CoreElt(core, estimate)        = NULL;  /* Set in _braid_InitHierarchy */

   /* Timers for key parts of code */
   _braid_CoreElt(core, timings) = 1;
   _braid_CoreElt(core, timer_coarse_solve) = 0.0;
   _braid_CoreElt(core, timer_drive_init)  = 0.0;
   _braid_CoreElt(core, timer_user_step)  = 0.0;
   _braid_CoreElt(core, timer_user_init)  = 0.0;
   _braid_CoreElt(core, timer_user_clone)  = 0.0;
   _braid_CoreElt(core, timer_user_free)  = 0.0;
   _braid_CoreElt(core, timer_user_sum)  = 0.0;
   _braid_CoreElt(core, timer_user_spatialnorm)  = 0.0;
   _braid_CoreElt(core, timer_user_access)  = 0.0;
   _braid_CoreElt(core, timer_user_sync)  = 0.0;
   _braid_CoreElt(core, timer_user_bufsize)  = 0.0;
   _braid_CoreElt(core, timer_user_bufpack)  = 0.0;
   _braid_CoreElt(core, timer_user_bufunpack)  = 0.0;
   _braid_CoreElt(core, timer_user_residual)  = 0.0;
   _braid_CoreElt(core, timer_user_scoarsen)  = 0.0;
   _braid_CoreElt(core, timer_user_srefine)  = 0.0;
   _braid_CoreElt(core, timer_user_inner_prod)  = 0.0;
   _braid_CoreElt(core, timer_MPI_recv)  = 0.0;
   _braid_CoreElt(core, timer_MPI_wait)  = 0.0;
   _braid_CoreElt(core, timer_MPI_wait_coarse)  = 0.0;
   _braid_CoreElt(core, timer_MPI_send)  = 0.0;
   _braid_CoreElt(core, timer_file_stem_len) = 13; /* length of file stem, not including the null terminator */
   _braid_CoreElt(core, timer_file_stem)  = malloc(14 * sizeof(char));
   sprintf(_braid_CoreElt(core, timer_file_stem), "braid_timings");

   braid_SetMaxLevels(core, max_levels);
   braid_SetMaxIter(core, max_iter);
   braid_SetPeriodic(core, 0);

   *core_ptr = core_dyn;

   return _braid_error_flag;
}

braid_Int
braid_Set_Info(MPI_Info info){
   DMR_SET_INFO(info);
   return _braid_error_flag;
}

braid_Core braid_dyn_get_original_core(braid_Core_dyn core_dyn) {
    if (core_dyn != NULL) {
         return _braid_CoreDynElt(core_dyn, original_core);
    } else {
        return NULL;
    }
}
