#ifndef dmr_HEADER
#define dmr_HEADER

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <signal.h>
#include <string.h>
#include <stdbool.h>

#define STR_SIZE 2048 /**< @brief Standard size for strings. **/

extern char _keys[][MPI_MAX_INFO_KEY];
extern char *keys[5];
extern MPI_Info info;
extern char main_pset[MPI_MAX_PSET_NAME_LEN];
extern char delta_pset[MPI_MAX_PSET_NAME_LEN];
extern char old_main_pset[MPI_MAX_PSET_NAME_LEN];
extern char final_pset[MPI_MAX_PSET_NAME_LEN];
extern char mpi_world_pset[];
extern char boolean_string[6];
extern int flag;
extern bool dynamic_proc;
extern bool primary_proc;
extern char **input_psets;
extern char **output_psets;
extern char **q_output_psets;
extern char str[MPI_MAX_PSET_NAME_LEN];
extern MPI_Group wgroup;
extern int noutput;
extern int noutput2;
extern int op_req;
extern int op_query;

extern MPI_Comm DMR_INTERCOMM;
extern MPI_Comm DMR_COMM_OLD;
extern MPI_Comm DMR_COMM_NEW;
extern int DMR_comm_rank;
extern int DMR_comm_size;
extern int DMR_comm_prev_size;
extern MPI_Session DMR_session;
extern MPI_Request psetop_req;

#define GET_MACRO(_0, _1, NAME, ...) NAME
#define DMR_FINALIZE(...) GET_MACRO(_0, ##__VA_ARGS__, DMR_FINALIZE1, DMR_FINALIZE0)(__VA_ARGS__)

#define DMR_INIT()                                                                                                                 \
    {                                                                                                                              \
        MPI_Session_init(MPI_INFO_NULL, MPI_ERRORS_ARE_FATAL, &DMR_session);                                                       \
        strcpy(main_pset, "mpi://WORLD");                                                                                          \
        strcpy(delta_pset, "");                                                                                                    \
        MPI_Session_get_pset_info(DMR_session, main_pset, &info);                                                                  \
        MPI_Info_get(info, "mpi_dyn", 6, boolean_string, &flag);                                                                   \
        if (flag && 0 == strcmp(boolean_string, "True"))                                                                           \
        {                                                                                                                          \
            MPI_Info_get(info, "mpi_primary", 6, boolean_string, &flag);                                                           \
            if (flag && 0 == strcmp(boolean_string, "True"))                                                                       \
            {                                                                                                                      \
                primary_proc = true;                                                                                               \
            }                                                                                                                      \
            MPI_Info_free(&info);                                                                                                  \
            dynamic_proc = true;                                                                                                   \
            MPI_Session_get_pset_data(DMR_session, main_pset, main_pset, (char **)&keys[2], 1, true, &info);                       \
            MPI_Info_get(info, "inter_pset", MPI_MAX_PSET_NAME_LEN, main_pset, &flag);                                             \
            MPI_Info_free(&info);                                                                                                  \
            if (0 == strcmp(main_pset, "dmr://finalize"))                                                                          \
            {                                                                                                                      \
                if (primary_proc)                                                                                                  \
                {                                                                                                                  \
                    MPI_Info_create(&info);                                                                                        \
                    MPI_Info_set(info, "dmr://finalize", "ack");                                                                   \
                    MPI_Session_set_pset_data(DMR_session, "mpi://WORLD", info);                                                   \
                    MPI_Info_free(&info);                                                                                          \
                }                                                                                                                  \
                MPI_Session_finalize(&DMR_session);                                                                                \
                return 0;                                                                                                          \
            }                                                                                                                      \
            MPI_Session_get_pset_info(DMR_session, main_pset, &info);                                                              \
        }                                                                                                                          \
        MPI_Info_get(info, "mpi_primary", 6, boolean_string, &flag);                                                               \
        if (flag && 0 == strcmp(boolean_string, "True"))                                                                           \
        {                                                                                                                          \
            primary_proc = true;                                                                                                   \
        }                                                                                                                          \
        else                                                                                                                       \
        {                                                                                                                          \
            primary_proc = false;                                                                                                  \
        }                                                                                                                          \
        MPI_Info_free(&info);                                                                                                      \
        MPI_Group wgroup = MPI_GROUP_NULL;                                                                                         \
        MPI_Group_from_session_pset(DMR_session, main_pset, &wgroup);                                                              \
        MPI_Comm_create_from_group(wgroup, "mpi.forum.example", MPI_INFO_NULL, MPI_ERRORS_RETURN, &DMR_INTERCOMM);                 \
        DMR_COMM_NEW = DMR_INTERCOMM;                                                                                              \
        MPI_Comm_rank(DMR_INTERCOMM, &DMR_comm_rank);                                                                              \
        MPI_Comm_size(DMR_INTERCOMM, &DMR_comm_size);                                                                              \
        MPI_Group_free(&wgroup);                                                                                                   \
        strcpy(final_pset, main_pset);                                                                                             \
    }

#define DMR_RECONFIGURATION(finalize_flag)                                                                                                               \
    {                                                                                                                                                    \
            noutput = 0;                                                                                                                                 \
            if (primary_proc && psetop_req == MPI_REQUEST_NULL)                                                                                          \
            {                                                                                                                                            \
                input_psets = (char **)malloc(1 * sizeof(char *));                                                                                       \
                input_psets[0] = strdup(main_pset);                                                                                                      \
                op_req = MPI_PSETOP_REPLACE;                                                                                                             \
                MPI_Session_dyn_v2a_psetop_nb(DMR_session, &op_req, input_psets, 1, &output_psets, &noutput, info, &psetop_req);                         \
                MPI_Info_free(&info);                                                                                                                    \
            }                                                                                                                                            \
            /* Query if there is a resource change */                                                                                                    \
            noutput2 = 0;                                                                                                                                \
            MPI_Session_dyn_v2a_query_psetop(DMR_session, main_pset, main_pset, &op_query, &q_output_psets, &noutput2);                                  \
            if (MPI_PSETOP_NULL != op_query)                                                                                                             \
            {                                                                                                                                            \
                MPI_Comm_dup(DMR_INTERCOMM, &DMR_COMM_OLD);                                                                                              \
                DMR_COMM_NEW = MPI_COMM_NULL;                                                                                                            \
                /* Publish name of inter-pset */                                                                                                         \
                if (primary_proc)                                                                                                                        \
                {                                                                                                                                        \
                    MPI_Request_free(&psetop_req);                                                                                                       \
                    MPI_Info_create(&info);                                                                                                              \
                    MPI_Info_set(info, "inter_pset", output_psets[2]);                                                                                   \
                    MPI_Session_set_pset_data(DMR_session, output_psets[1], info);                                                                       \
                    MPI_Info_free(&info);                                                                                                                \
                    free_string_array(output_psets, noutput);                                                                                            \
                    free_string_array(input_psets, 1);                                                                                                   \
                }                                                                                                                                        \
                strcpy(final_pset, q_output_psets[0]);                                                                                                   \
                MPI_Session_get_pset_data(DMR_session, main_pset, q_output_psets[1], (char **)&keys[2], 1, true, &info);                                 \
                strcpy(old_main_pset, main_pset);                                                                                                        \
                MPI_Info_get(info, "inter_pset", MPI_MAX_PSET_NAME_LEN, main_pset, &flag);                                                               \
                MPI_Info_free(&info);                                                                                                                    \
                /* Check if this process is the primary process of the new pset */                                                                       \
                MPI_Session_get_pset_info(DMR_session, main_pset, &info);                                                                                \
                MPI_Info_get(info, "mpi_primary", 6, boolean_string, &flag);                                                                             \
                primary_proc = (0 == strcmp(boolean_string, "True"));                                                                                    \
                /* Check if this process is included in the new PSet */                                                                                  \
                MPI_Info_get(info, "mpi_included", 6, boolean_string, &flag);                                                                            \
                MPI_Info_free(&info);                                                                                                                    \
                /* Create new communicator */                                                                                                            \
                MPI_Comm_disconnect(&DMR_INTERCOMM);                                                                                                     \
                /*printf("Before New Intercomm [%d/%d] %d: %s(%s,%d)\n", DMR_comm_rank, DMR_comm_size, getpid(), __FILE__, __func__, __LINE__);       */ \
                if (0 != strcmp(boolean_string, "False"))                                                                                                \
                {                                                                                                                                        \
                    strcpy(final_pset, main_pset);                                                                                                       \
                    MPI_Group_from_session_pset(DMR_session, main_pset, &wgroup);                                                                        \
                    MPI_Comm_create_from_group(wgroup, "mpi.forum.example", MPI_INFO_NULL, MPI_ERRORS_RETURN, &DMR_INTERCOMM);                           \
                    MPI_Group_free(&wgroup);                                                                                                             \
                    MPI_Comm_rank(DMR_INTERCOMM, &DMR_comm_rank);                                                                                        \
                    MPI_Comm_size(DMR_INTERCOMM, &DMR_comm_size);                                                                                        \
                    DMR_COMM_NEW = DMR_INTERCOMM;                                                                                                        \
                }                                                                                                                                        \
                free_string_array(q_output_psets, noutput2);                                                                                             \
                /* Finalize PSetOp*/                                                                                                                     \
                if (primary_proc)                                                                                                                        \
                {                                                                                                                                        \
                    MPI_Session_dyn_finalize_psetop(DMR_session, old_main_pset);                                                                         \
                }                                                                                                                                        \
                MPI_Comm_disconnect(&DMR_COMM_OLD);                                                                                                      \
                if (0 == strcmp(boolean_string, "False"))                                                                                                \
                {                                                                                                                                        \
                    /* Leaving Processes */                                                                                                              \
                    finalize_flag = 1;                                                                                                                   \
                }                                                                                                                                        \
            }                                                                                                                                            \
    }

#define DMR_FINALIZE1(FUNC_FINALIZE)                                                                                           \
    {                                                                                                                          \
        FUNC_FINALIZE;                                                                                                         \
        /* We need to cancel our last pset operation. If it was already executed we need to tell any new procs to terminate */ \
        if (primary_proc)                                                                                                      \
        {                                                                                                                      \
            op_req = MPI_PSETOP_CANCEL;                                                                                        \
            input_psets = (char **)malloc(1 * sizeof(char *));                                                                 \
            input_psets[0] = strdup(main_pset);                                                                                \
            noutput = 0;                                                                                                       \
            MPI_Session_dyn_v2a_psetop(DMR_session, &op_req, input_psets, 1, &output_psets, &noutput, MPI_INFO_NULL);          \
            free_string_array(input_psets, 1);                                                                                 \
            if (MPI_PSETOP_NULL == op_req)                                                                                     \
            {                                                                                                                  \
                strcpy(delta_pset, "mpi://SELF");                                                                              \
                MPI_Session_dyn_v2a_query_psetop(DMR_session, delta_pset, main_pset, &op_query, &q_output_psets, &noutput);    \
                if (MPI_PSETOP_NULL != op_query)                                                                               \
                {                                                                                                              \
                    if (0 != strcmp("", q_output_psets[1]))                                                                    \
                    {                                                                                                          \
                        MPI_Info_create(&info);                                                                                \
                        MPI_Info_set(info, "inter_pset", "dmr://finalize");                                                    \
                        MPI_Session_set_pset_data(DMR_session, q_output_psets[1], info);                                       \
                        MPI_Info_free(&info);                                                                                  \
                    }                                                                                                          \
                    MPI_Session_dyn_finalize_psetop(DMR_session, main_pset);                                                   \
                    MPI_Session_get_pset_data(DMR_session, delta_pset, q_output_psets[1], (char **)&keys[4], 1, true, &info);  \
                    free_string_array(q_output_psets, noutput);                                                                \
                    MPI_Info_free(&info);                                                                                      \
                }                                                                                                              \
            }                                                                                                                  \
        }                                                                                                                      \
        input_psets = (char **)malloc(1 * sizeof(char *));                                                                     \
        input_psets[0] = strdup(final_pset);                                                                                   \
        MPI_Session_pset_barrier(DMR_session, input_psets, 1, NULL);                                                           \
        free_string_array(input_psets, 1);                                                                                     \
        if (MPI_COMM_NULL != DMR_INTERCOMM)                                                                                    \
        {                                                                                                                      \
            MPI_Comm_disconnect(&DMR_INTERCOMM);                                                                               \
        }                                                                                                                      \
        MPI_Session_finalize(&DMR_session);                                                                                    \
    }

#define DMR_FINALIZE0()                                                                                                        \
    {                                                                                                                          \
        /* printf("test start pid: %d\n", getpid()); */                                                                             \
        /* We need to cancel our last pset operation. If it was already executed we need to tell any new procs to terminate */ \
        if (primary_proc)                                                                                                      \
        {                                                                                                                      \
            /* printf("test if primary pid: %d\n", getpid());      */                                                               \
            op_req = MPI_PSETOP_CANCEL;                                                                                        \
            input_psets = (char **)malloc(1 * sizeof(char *));                                                                 \
            input_psets[0] = strdup(main_pset);                                                                                \
            noutput = 0;                                                                                                       \
            MPI_Session_dyn_v2a_psetop(DMR_session, &op_req, input_psets, 1, &output_psets, &noutput, MPI_INFO_NULL);          \
            free_string_array(input_psets, 1);                                                                                 \
            if (MPI_PSETOP_NULL == op_req)                                                                                     \
            {                                                                                                                  \
                strcpy(delta_pset, "mpi://SELF");                                                                              \
                MPI_Session_dyn_v2a_query_psetop(DMR_session, delta_pset, main_pset, &op_query, &q_output_psets, &noutput);    \
                if (MPI_PSETOP_NULL != op_query)                                                                               \
                {                                                                                                              \
                    if (0 != strcmp("", q_output_psets[1]))                                                                    \
                    {                                                                                                          \
                        MPI_Info_create(&info);                                                                                \
                        MPI_Info_set(info, "inter_pset", "dmr://finalize");                                                    \
                        MPI_Session_set_pset_data(DMR_session, q_output_psets[1], info);                                       \
                        MPI_Info_free(&info);                                                                                  \
                        MPI_Session_dyn_finalize_psetop(DMR_session, main_pset);                                                   \
                        MPI_Session_get_pset_data(DMR_session, delta_pset, q_output_psets[1], (char **)&keys[4], 1, true, &info);  \
                        free_string_array(q_output_psets, noutput);                                                                \
                        MPI_Info_free(&info);                                                                                      \
                    }                                                                                                          \
                }                                                                                                              \
            }                                                                                                                  \
        }                                                                                                                      \
        input_psets = (char **)malloc(1 * sizeof(char *));                                                                     \
        input_psets[0] = strdup(final_pset);                                                                                   \
        /*printf("test pid: %d\n", getpid());    */                                                                                \
        fflush(NULL);                                                                                                          \
        MPI_Session_pset_barrier(DMR_session, input_psets, 1, NULL);                                                           \
        /*printf("test2 pid: %d\n", getpid());    */                                                                               \
        free_string_array(input_psets, 1);                                                                                     \
        if (MPI_COMM_NULL != DMR_INTERCOMM)                                                                                    \
        {                                                                                                                      \
            MPI_Comm_disconnect(&DMR_INTERCOMM);                                                                               \
        }                                                                                                                      \
        MPI_Session_finalize(&DMR_session);                                                                                    \
    }

void free_string_array(char **array, int size);

int DMR_Reconfiguration(char *argv[], MPI_Comm *DMR_INTERCOMM, int min, int max, int step, int pref);
void DMR_Send_expand(double *data, int size, MPI_Comm DMR_INTERCOM);
void DMR_Recv_expand(double **data, int *size, MPI_Comm DMR_INTERCOM);
void DMR_Send_shrink(double *data, int size, MPI_Comm DMR_INTERCOM);
void DMR_Recv_shrink(double **data, int *size, MPI_Comm DMR_INTERCOM);
void DMR_Set_parameters(MPI_Info mpi_info);

#endif