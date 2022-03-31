#pragma once

#include <concepts>
#include <type_traits>

#include "core/mpi_flags.h"

namespace cvdf::parallel
{
    template <class T> concept parallel_group = requires(T t)
    {
        t.rank();
        t.size();
        t.sync();
    };
    
    static mpi_data_t get_data_type(double d)      { return MPI_DOUBLE; }
    static mpi_data_t get_data_type(char d)        { return MPI_CHAR; }
    static mpi_data_t get_data_type(float d)       { return MPI_FLOAT; }
    static mpi_data_t get_data_type(int d)         { return MPI_INT; }
    static mpi_data_t get_data_type(std::size_t d) { return MPI_UINT64_T; }
    
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
            const mpi_t& sync(void) const
            {
                MPI_CHECK(MPI_Barrier(this->channel));
                return *this;
            }
            
            template <typename data_t> request_t async_recv(data_t* buf, int count, int source)
            {
                auto dtype = get_data_type(data_t());
                request_t p;
                MPI_CHECK(MPI_Irecv(buf, count, dtype, source, 1, this->channel, &p));
                return p;
            }
            
            template <typename data_t> void sync_send(const data_t* buf, int count, int dest)
            {
                auto dtype = get_data_type(data_t());
                MPI_CHECK(MPI_Ssend(buf, count, dtype, dest, 1, this->channel));
            }
            
            void await_all(int count, request_t reqs[], status_t stats[])
            {
                MPI_CHECK(MPI_Waitall(count, reqs, stats));
            }
            
            template <typename... data_t> auto sum(data_t... datas) const
            {
                auto sum_loc = (... + datas);
                decltype(sum_loc) sum_glob;
                auto dtype = get_data_type(sum_loc);
                MPI_CHECK(MPI_Allreduce(&sum_loc, &sum_glob, 1, dtype, MPI_SUM, this->channel));
                return sum_glob;
            }
            
            template <typename data_t> auto reduce(const data_t& data, mpi_op_t op) const
            {
                data_t data_loc = data;
                auto dtype = get_data_type(data);
                data_t red_glob;
                MPI_CHECK(MPI_Allreduce(&data_loc, &red_glob, 1, dtype, op, this->channel));
                return red_glob;
            }
            
        private:
            int g_rank, g_size;
            mpi_comm_t channel;
    };
}