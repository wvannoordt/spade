#pragma once

#ifndef MPI_ENABLE
#define MPI_ENABLE 1
#endif

#if(MPI_ENABLE)
#include "mpi.h"
#endif

#if(MPI_ENABLE)
typedef MPI_Comm      mpi_comm_t;
typedef MPI_Datatype  mpi_data_t;
typedef MPI_Op        mpi_op_t;
typedef MPI_Request   request_t;
typedef MPI_Status    status_t;
typedef MPI_File      mpi_native_file_t;
#define MPI_CHECK(mycode) {mycode ;}

namespace spade::parallel
{
    static const mpi_op_t par_sum = MPI_SUM;
    static const mpi_op_t par_max = MPI_MAX;
}
#else
typedef unsigned long mpi_comm_t;
typedef unsigned long mpi_data_t;
typedef unsigned long mpi_op_t;
typedef unsigned long request_t;
typedef unsigned long status_t;
typedef unsigned long mpi_native_file_t;

static mpi_data_t MPI_DOUBLE    = 0;
static mpi_data_t MPI_CHAR      = 1;
static mpi_data_t MPI_FLOAT     = 2;
static mpi_data_t MPI_INT       = 3;
static mpi_data_t MPI_UINT64_T  = 4;
#define MPI_CHECK(mycode) ;

namespace spade::parallel
{
    static const mpi_op_t par_sum = 0;
    static const mpi_op_t par_max = 1;
}
#endif