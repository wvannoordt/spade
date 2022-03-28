#pragma once

#ifndef MPI_ENABLE
#define MPI_ENABLE 0
#else
#include "mpi.h"
#endif

#if(MPI_ENABLE)
typedef MPI_Comm      mpi_comm_t;
typedef MPI_Datatype  mpi_data_t;
typedef MPI_Op        mpi_op_t;
typedef MPI_Request   request_t;
typedef MPI_Status    status_t;
#define MPI_CHECK(mycode) {mycode ;}
#else
typedef unsigned long mpi_comm_t;
typedef unsigned long mpi_data_t;
typedef unsigned long mpi_op_t;
typedef unsigned long request_t;
typedef unsigned long status_t;
#define MPI_CHECK(mycode) ;
#endif