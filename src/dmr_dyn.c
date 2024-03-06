#include "dmr_dyn.h"

char _keys[][MPI_MAX_INFO_KEY] = {"mpi_dyn", "mpi_primary", "inter_pset", "mpi_included", "dmr://finalize"};
char *keys[5] = {_keys[0], _keys[1], _keys[2], _keys[3], _keys[4]};
MPI_Info info = MPI_INFO_NULL;
char main_pset[MPI_MAX_PSET_NAME_LEN];
char delta_pset[MPI_MAX_PSET_NAME_LEN];
char old_main_pset[MPI_MAX_PSET_NAME_LEN];
char final_pset[MPI_MAX_PSET_NAME_LEN];
char mpi_world_pset[] = "mpi://WORLD";
char boolean_string[6];
int flag;
bool dynamic_proc = false;
bool primary_proc = false;
char **input_psets = NULL, **output_psets = NULL, **q_output_psets = NULL;
char str[MPI_MAX_PSET_NAME_LEN];
MPI_Group wgroup = MPI_GROUP_NULL;
int noutput = 0, noutput2 = 0, op_req = MPI_PSETOP_NULL, op_query = MPI_PSETOP_NULL;

MPI_Comm DMR_INTERCOMM = MPI_COMM_NULL;
MPI_Comm DMR_COMM_OLD = MPI_COMM_NULL;
MPI_Comm DMR_COMM_NEW = MPI_COMM_NULL;
int DMR_comm_rank, DMR_comm_size, DMR_comm_prev_size;
MPI_Session DMR_session = MPI_SESSION_NULL;
MPI_Request psetop_req = MPI_REQUEST_NULL;

void DMR_Set_parameters(MPI_Info mpi_info){
    info = mpi_info;
}

void free_string_array(char **array, int size)
{
    int i;
    if (0 == size)
    {
        return;
    }
    for (i = 0; i < size; i++)
    {
        free(array[i]);
    }
    free(array);
}