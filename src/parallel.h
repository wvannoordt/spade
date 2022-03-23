#pragma once

#include <concepts>
#include <type_traits>

#include "mpi_flags.h"

namespace cvdf::parallel
{
    template <class T> concept parallel_group = requires(T t)
    {
        t.rank();
        t.size();
        t.sync();
    };
    
    class mpi_t
    {
        public:
            mpi_t(int* argc, char*** argv)
            {
                MPI_CHECK(MPI_Init(argc, argv));
#if(MPI_ENABLE)
                this->channel = MPI_COMM_WORLD;
#else
                this->channel = 0;
#endif
                this->g_rank = 0;
                this->g_size = 1;
                MPI_CHECK(MPI_Comm_rank(this->channel, &this->g_rank));
                MPI_CHECK(MPI_Comm_size(this->channel, &this->g_size));
            }
            ~mpi_t(void)
            {
#if(MPI_ENABLE)
                MPI_Finalize();
#endif
            }
            
            int    rank(void) const {return g_rank;}
            int    size(void) const {return g_size;}
            bool isroot(void) const {return g_rank==0;}
            mpi_t& sync(void)
            {
                MPI_CHECK(MPI_Barrier(this->channel));
                return *this;
            }
        private:
            int g_rank, g_size;
            mpi_comm_t channel;
    };
}