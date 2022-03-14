#pragma once

#ifndef MPI_ENABLE
#define MPI_ENABLE 0
#else
#include "mpi.h"
#endif

#if(MPI_ENABLE)
typedef MPI_Comm mpi_comm_t;
#define MPI_CHECK(mycode) {mycode ;}
#else
typedef int mpi_comm_t;
#define MPI_CHECK(mycode) ;
#endif