#ifndef _braid_dyn_HEADER
#define _braid_dyn_HEADER

#include "braid_dyn.h"
#include "_braid.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _braid_Core_struct_dyn
{
   braid_Real             interval_len;
   braid_PtFcnGetValue    getValue;
   braid_Int              max_procs;
   _braid_Core            *original_core;
}_braid_Core_dyn;

/** 
 * Accessor for _braid_Core_dyn attributes 
 **/
#define _braid_CoreDynElt(core_dyn, elt)     (  (core_dyn)  -> elt )

/** 
 * Accessor for _braid_Core_dyn functions 
 **/
#define _braid_CoreDynFcn(core_dyn, fcn)     (*((core_dyn)  -> fcn))


braid_Int
_braid_USetVector_Dyn(braid_Core        core,
                  braid_Int         level,
                  braid_Int         index,
                  braid_Vector  u,
                  braid_Int         move);

#ifdef __cplusplus
}
#endif

#include "status.h"
#include "base.h"
#include "adjoint.h"
#include "delta.h"

#endif

