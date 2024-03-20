#include "_braid_dyn.h"
#include "util.h"

/*----------------------------------------------------------------------------
 * Stores the u-vector on grid 'level' at point 'index'.  If 'index' is my "send
 * index", a send is initiated to a neighbor processor.  If 'move' is true, the
 * u-vector is moved into core storage instead of copied.  If the u-vector is
 * not stored, nothing is done or only the shell is copied/moved when the shellvector
 * feature is used.
 *----------------------------------------------------------------------------*/

braid_Int
_braid_USetVector_Dyn(braid_Core        core,
                  braid_Int         level,
                  braid_Int         index,
                  braid_Vector  u_,
                  braid_Int         move)
{
   //newDyn
   braid_BaseVector u = _braid_TAlloc(_braid_BaseVector, 1);
   u->userVector = u_;
   u->bar        = _braid_TAlloc(_braid_VectorBar, 1);
   u->basis      = NULL;



   braid_App            app         = _braid_CoreElt(core, app);
   _braid_Grid        **grids       = _braid_CoreElt(core, grids);
   braid_BaseVector    *ua          = _braid_GridElt(grids[level], ua);
   braid_Int            send_index  = _braid_GridElt(grids[level], send_index);
   _braid_CommHandle   *send_handle = _braid_GridElt(grids[level], send_handle);
   braid_Int            done        = _braid_CoreElt(core, done);
   braid_Int            gupper      = _braid_CoreElt(core, gupper);
   braid_Int            storage     = _braid_CoreElt(core, storage);
   braid_Int            cfactor     = _braid_GridElt(grids[level], cfactor);
   braid_Int            iu, sflag;

   u->bar->useCount = 1;
   // _braid_CoreFcn(core, init)(app, 1500.0, &(u->bar->userVector));
   u->bar->userVector = u_;

   if (index == send_index)
   {
      /* Post send to neighbor processor */
      _braid_CommSendInit(core, level, index, u, &send_handle);
      _braid_GridElt(grids[level], send_index)  = _braid_SendIndexNull;
      _braid_GridElt(grids[level], send_handle) = send_handle;
   }

   _braid_UGetIndex(core, level, index, &iu, &sflag);
   if (sflag == 0) // We have a full point
   {
      if (ua[iu] != NULL)
      {
         _braid_BaseFree(core, app,  ua[iu]);
      }
      if (move)
      {
         ua[iu] = u;                                   /* move the vector */
      }
      else
      {
         _braid_BaseClone(core, app,  u, &ua[iu]); /* copy the vector */
      }
   }
   else if (sflag == -1) // We have a shell
   {
      if (ua[iu] != NULL)
      {
         _braid_BaseFree(core, app,  ua[iu]);
      }
      if (move)
      {
         // We are on an F-point, with shellvector option. We only keep the shell.
         _braid_BaseSFree(core,  app, u);
         ua[iu] = u;                                   /* move the vector */
      }
      else
      {
         _braid_BaseSClone(core,  app, u, &ua[iu]); /* copy the vector */
      }
   }
   else if (move) // We store nothing
   {
      _braid_BaseFree(core, app,  u);              /* free the vector */
   }

   /* If braid is finished, make sure the last time point is stored, i.e.,
    * store the last time point in ulast if storage is not enabled for F-points
    * OR the last time point is not a C-point */

   if (done && (index == gupper) && (level == 0) && (storage < 0) && !(_braid_IsCPoint(gupper, cfactor)))
   {
      if (_braid_GridElt(grids[level], ulast) != NULL)
      {
         _braid_BaseFree(core, app, _braid_GridElt(grids[level], ulast));
         _braid_GridElt(grids[level], ulast) = NULL;
      }
      _braid_BaseClone(core, app,  u, &(_braid_GridElt(grids[level], ulast)));
  }

// newDyn
   // free(u->bar);
   // free(u);

   return _braid_error_flag;
}
