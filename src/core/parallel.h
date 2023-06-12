#pragma once

#include <type_traits>
#include <vector>

#include "core/config.h"
#include "core/config.h"
#include "core/print.h"
#include "core/mpi_flags.h"
#include "core/aliases.h"

#if (USE_SOURCE_LOCATION)
#include <source_location>
#endif

namespace spade::parallel
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
    
    const int default_root_rank = 0;
    
    class mpi_t
    {
        public:
            mpi_t(int* argc, char*** argv, bool wrap_init_in = true)
            {
                wrap_init = wrap_init_in;
                if (wrap_init)
                {
                    MPI_CHECK(MPI_Init(argc, argv));
                }
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
                if (wrap_init)
                {
                    MPI_Finalize();
                }
#endif
            }
            
            int    rank(void) const {return g_rank;}
            int    size(void) const {return g_size;}
            bool isroot(void) const {return g_rank==default_root_rank;}
            const mpi_t& sync(void) const
            {
                MPI_CHECK(MPI_Barrier(this->channel));
                return *this;
            }

#if(USE_SOURCE_LOCATION)
            const mpi_t& pause(const std::source_location location = std::source_location::current()) const
#else
            const mpi_t& pause(void) const
#endif
            {
                if (isroot())
                {
#if(USE_SOURCE_LOCATION)
                    print("Parallel group", (void*)this, " paused at:");
                    print("File:", location.file_name());
                    print("Line:", location.line());
#endif
                    print("Awaiting input...");
                    std::cin.get();
                }
                sync();
                return (*this);
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

            template <typename data_t> auto sum(const std::vector<data_t>& data) const
            {
                data_t sumval = data_t();
                for (const auto& e: data) sumval += e;
                return this->sum(sumval);
            }
            
            template <typename data_t> void append_root(aliases::vector<data_t>& root_result, const data_t& local_data) const
            {
                if (this->isroot()) root_result.resize(this->g_size);
                auto dtype = get_data_type(data_t());
                MPI_CHECK(MPI_Gather(&local_data, 1, dtype, &root_result[0], 1, dtype, default_root_rank, this->channel));
            }
            
            template <typename data_t> void append_root(aliases::vector<data_t>& root_result, const aliases::vector<data_t>& local_data) const
            {
                std::size_t output_size = this->sum(local_data.size());
                if (this->isroot()) root_result.resize(output_size);
                aliases::vector<int> sizes;
                this->append_root(sizes, (int)local_data.size());
                aliases::vector<int> offset;
                if (this->isroot())
                {
                    offset.resize(sizes.size());
                    offset[0] = 0;
                    for (int i = 1; i < this->g_size; ++i)
                    {
                        offset[i] = offset[i-1] + sizes[i-1];
                    }
                }
                auto dtype = get_data_type(data_t());
                MPI_CHECK(MPI_Gatherv(&local_data[0], local_data.size(), dtype, &root_result[0], &sizes[0], &offset[0], dtype, default_root_rank, this->channel));
            }
            
            template <typename data_t> auto reduce(const data_t& data, mpi_op_t op) const
            {
                data_t data_loc = data;
                auto dtype = get_data_type(data);
                data_t red_glob;
                MPI_CHECK(MPI_Allreduce(&data_loc, &red_glob, 1, dtype, op, this->channel));
                return red_glob;
            }
            
            mpi_comm_t get_channel(void) const {return channel;}
            
        private:
            int g_rank, g_size;
            mpi_comm_t channel;
            bool wrap_init; 
    };
    
    struct par_buf_t
    {
        aliases::vector<void*>       datas;
        aliases::vector<std::size_t> sizes;
        aliases::vector<std::size_t> offst;
        void add(void* base, const std::size_t& size, const std::size_t offset)
        {
            sizes.push_back(size);
            datas.push_back(base);
            offst.push_back(offset);
        }
        
        void clear(void)
        {
            datas.clear();
            sizes.clear();
            offst.clear();
        }
    };
    
    class mpi_file_t
    {
        public:
            mpi_file_t(const mpi_t& group_in, const std::string& filename)
            {
                this->open(filename);
                this->group = &group_in;
            }
            
            mpi_file_t(const mpi_t& group_in)
            {
                this->group = &group_in;
            }
            
            ~mpi_file_t(void)
            {
                if (this->is_open)
                {
                    this->close();
                }
            }
            
            mpi_file_t& open(const std::string& filename)
            {
                is_open = true;
                // mpi_comm_t chan = group->get_channel();
                auto cstr = filename.c_str();
                MPI_CHECK(MPI_File_open(MPI_COMM_WORLD, cstr, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &file_handle));
                return *this;
            }
            
            mpi_file_t& write_buf(par_buf_t& buf)
            {
                for (std::size_t i = 0; i < buf.offst.size(); ++i)
                {
                    MPI_CHECK(MPI_File_write_at(file_handle, buf.offst[i], buf.datas[i], buf.sizes[i], MPI_CHAR, &file_status));
                }
                return *this;
            }
            
            mpi_file_t& read_buf(par_buf_t& buf)
            {
                for (std::size_t i = 0; i < buf.offst.size(); ++i)
                {
                    MPI_CHECK(MPI_File_read_at(file_handle, buf.offst[i], buf.datas[i], buf.sizes[i], MPI_CHAR, &file_status));
                }
                return *this;
            }
            
            mpi_file_t& close(void)
            {
                is_open = false;
                MPI_CHECK(MPI_File_close(&file_handle));
                return *this;
            }
        private:        
            bool is_open = false;
            mpi_native_file_t file_handle;
            status_t file_status;
            const mpi_t* group;
    };
}
