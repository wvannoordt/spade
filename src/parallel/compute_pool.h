#pragma once

#include <memory>
#include <algorithm>
#include <mutex>
#include <variant>
#include <thread>
#include <barrier>

#include "mpi.h"
#include "core/except.h"
#include "core/vec_image.h"

namespace spade::parallel
{
    using mpi_comm_t = decltype(MPI_COMM_WORLD);
    inline void mpi_check(const int code)
    {
        if (code != MPI_SUCCESS)
        {
            char* buf;
            int len;
            if (MPI_Error_string(code, buf, &len) == MPI_SUCCESS)
            {
                std::string msg_str;
                msg_str.assign(buf, len);
                std::string fail_message = "MPI error occured: " + msg_str;
                throw except::sp_exception(fail_message);
            }
            else
            {
                std::string fail_message = "MPI error occured: unspecified error";
                throw except::sp_exception(fail_message);
            }
        }
    }
    
    struct proc_id_t
    {
        int node = 9999, thread = 9999, num_node = 9999, num_thread = 9999, rank = 9999, num_rank = 9999;
    };
    
    struct compute_env_t;
    
    using communicable_t = std::variant<
        double,      double*,
        char,        char*,
        float,       float*,
        int,         int*,
        std::size_t, std::size_t*
        >;
    
    using message_buf_t = std::variant<
        utils::vec_image_t<double>,      utils::const_vec_image_t<double>,
        utils::vec_image_t<char>,        utils::const_vec_image_t<char>,
        utils::vec_image_t<float>,       utils::const_vec_image_t<float>,
        utils::vec_image_t<int>,         utils::const_vec_image_t<int>,
        utils::vec_image_t<std::size_t>, utils::const_vec_image_t<std::size_t>
        >;
    
    using barrier_t = std::barrier<void(*)()>;
    //Triple vector, nice
    using mailbox_t = std::vector<std::vector<std::vector<message_buf_t>>>;
    struct common_thread_data_t
    {
        mutable std::mutex mut;
        mutable barrier_t barrier;
        mutable std::vector<communicable_t> buffer;
        
        mutable mailbox_t outbox;
        mutable mailbox_t inbox;
        
        std::vector<proc_id_t> rank_pids;
        bool p2p_enabled;
    };
    
    using  mpi_data_t = MPI_Datatype;
    static mpi_data_t get_mpi_type(const double&)      { return MPI_DOUBLE; }
    static mpi_data_t get_mpi_type(const char&)        { return MPI_CHAR; }
    static mpi_data_t get_mpi_type(const float&)       { return MPI_FLOAT; }
    static mpi_data_t get_mpi_type(const int&)         { return MPI_INT; }
    static mpi_data_t get_mpi_type(const std::size_t&) { return MPI_UINT64_T; }
    
    struct pool_t
    {
        proc_id_t pool_pid;
        mpi_comm_t comm;
        const common_thread_data_t& env;
        
        template <typename buffer_type>
        void post_send(const buffer_type& buf, int recv)
        {
            message_buf_t bf = buf;
            env.outbox[this->rank()][recv].push_back(bf);
        }
        
        template <typename buffer_type>
        void post_recv(const buffer_type& buf, int send)
        {
            message_buf_t bf = buf;
            env.inbox[this->rank()][send].push_back(bf);
        }
        
        template <typename buf_t>
        void send_all()
        {
            // Need to clear messages after we finish
            this->sync();
            
            // Receivers are the ones responsible for copying.
            int myrank = this->rank();
            auto& inbox = env.inbox[myrank];
            for (int donor = 0; donor < this->size(); ++donor)
            {
                auto my_id    = this->pid();
                auto their_id = env.rank_pids[donor];
                if (my_id.node == their_id.node)
                {
                    auto& recv_list = inbox[donor];
                    auto& send_list = env.outbox[donor][myrank];
                    std::string msg0 = "invalid communication: list size for sender and receiver do not match";
                    std::string msg1 = "invalid communication: message size mismatch";
                    if (recv_list.size() != send_list.size()) throw except::sp_exception(msg0);
                    for (int i_msg = 0; i_msg < recv_list.size(); ++i_msg)
                    {
                        auto recv_buf = std::get<buf_t>(recv_list[i_msg]);
                        auto send_buf = std::get<buf_t>(send_list[i_msg]);
                        if (recv_buf.size() != send_buf.size()) throw except::sp_exception(msg1);
                        std::copy(&send_buf[0], &send_buf[0] + send_buf.size(), &recv_buf[0]);
                    }
                }
                else
                {
                    std::string msg = "invalid communication: inter-process communication not yet working";
                    throw except::sp_exception(msg);
                }
            }
            
            auto& outbox = env.outbox[myrank];
            
            this->sync();
            
            for (auto& list:  inbox) list.clear();
            for (auto& list: outbox) list.clear();
        }
        
        bool p2p_enabled() const { return env.p2p_enabled; }
        
        void set_buf_size(const std::size_t& bsize) const
        {
            if (this->isbase())
            {
                std::lock_guard<std::mutex> guard(env.mut);
                env.buffer.resize(bsize);
            }
            env.barrier.arrive_and_wait();
        }
        
        template <typename data_t>
        requires (spade::utils::variant_contains<data_t, communicable_t>)
        void post(const data_t& data) const
        {
            this->set_buf_size(pool_pid.num_thread);
            int idx = this->thread();
            env.buffer[idx] = data;
            this->sync();
        }
        
        template <typename data_t>
        requires (spade::utils::variant_contains<data_t, communicable_t>)
        auto thread_sum(const data_t& local_val) const
        {
            this->post(local_val);
            data_t sum = 0;
            for (const auto& v: env.buffer)
            {
                sum += std::get<data_t>(v);
            }
            return sum;
        }
        
        template <typename data_t>
        data_t thread_broadcast(const int thread_id, const data_t value) const
        {
            this->post(value);
            return std::get<data_t>(env.buffer[thread_id]);
        }
        
        template <typename data_t>
        auto sum(const data_t& val) const
        {
            data_t sum_thr  = this->thread_sum(val);
            data_t sum_glob = 0;
            auto dtype = get_mpi_type(sum_thr);
            auto base = this->base_id();
            if (this->isbase())
            {
                mpi_check(MPI_Allreduce(&sum_thr, &sum_glob, 1, dtype, MPI_SUM, this->comm));
            } 
            return this->thread_broadcast(base, sum_glob);
        }
        
        const proc_id_t& pid() const { return pool_pid; }
        int    rank()  const { return this->pid().rank; }
        int    size()  const { return this->pid().num_rank; }
        int  thread()  const { return this->pid().thread; }
        int  base_id() const { return 0; }
        int  root_id() const { return 0; }
        bool isroot()  const { return this->rank()   == this->root_id(); }
        bool isbase()  const { return this->thread() == this->base_id(); }
        bool is_root() const { return this->isroot(); } //This is why we have proper naming conventions
        bool is_base() const { return this->isbase(); } //This is *also* why we have proper naming conventions
        
        const pool_t& pause() const
        {
            if (this->is_root()) std::cin.get();
            this->sync();
            return *this;
        }
        
        const pool_t& sync() const
        {
            // MPI node barrier
            if (this->isbase())
            {
                mpi_check(MPI_Barrier(this->comm));
            }
            
            // std::thread barrier
            env.barrier.arrive_and_wait();
            return *this;
        }
        
        template <typename func_t>
        const pool_t& sequential(const func_t& func) const
        {
            for (int i = 0; i < size(); ++i)
            {
                if (i == rank())
                {
                    func();
                }
                this->sync();
            }
            return *this;
        }
    };
    
    struct mpi_context_t
    {
        mpi_comm_t default_comm;
        mpi_context_t(int* argc, char*** argv)
        {
            int required = MPI_THREAD_FUNNELED;
            int provided;
            mpi_check(MPI_Init_thread(argc, argv, required, &provided));
            if (provided != required)
            {
                std::string fail_message = "MPI was unable to initialize with threading support. Please compile MPI with threading support.";
                throw except::sp_exception(fail_message);
            }
            default_comm = MPI_COMM_WORLD;
        }
        
        ~mpi_context_t()
        {            
            MPI_Finalize();
        }
    };
    
    std::ostream& operator << (std::ostream& os, const proc_id_t& pid)
    {
        os << "[pid] thread: {" << pid.thread << "/" << pid.num_thread
            << "} node: {" << pid.node << "/" << pid.num_node
            << "} rank: {" << pid.rank << "/" << pid.num_rank << "}";
        return os;
    }
    
    void nil() noexcept {}
    struct compute_env_t
    {
        std::shared_ptr<mpi_context_t> context;
        
        std::vector<std::thread> others;
        std::vector<proc_id_t> ids;
        proc_id_t root;
        
        std::vector<int> devices;
        
        //Thread shared memory
        common_thread_data_t data;
        
        compute_env_t(int* argc, char*** argv, const std::vector<int>& devices_in = {0})
        : devices{devices_in},
        data{
            std::mutex(),
            barrier_t(devices_in.size(), nil),
            std::vector<communicable_t>(),
            mailbox_t(),
            mailbox_t(),
            std::vector<proc_id_t>(),
            false}
        {
            std::sort(devices.begin(), devices.end());
            if (std::unique(devices.begin(), devices.end()) != devices.end())
            {
                throw except::sp_exception("Duplicate device in launch configuration!");
            }
            int num_threads = devices_in.size();
            context = std::make_shared<mpi_context_t>(argc, argv);
            mpi_check(MPI_Comm_rank(context->default_comm, &root.node));
            mpi_check(MPI_Comm_size(context->default_comm, &root.num_node));
            
            if (num_threads < 1)
            {
                throw except::sp_exception("Cannot execute with less than 1 thread!");
            }
            others.resize(num_threads-1);
            ids.resize(num_threads-1);
            root.thread     = 0;
            root.num_thread = 1 + others.size();
            int ithr = 1;
            for (auto& id: ids)
            {
                id.node       = root.node;
                id.num_node   = root.num_node;
                id.thread     = ithr++;
                id.num_thread = root.num_thread;
            }
            std::vector<int> each_threads(root.num_node, 0);
            each_threads[root.node] = num_threads;
            mpi_check(MPI_Allgather(&num_threads, 1, MPI_INT, &each_threads[0], 1, MPI_INT, context->default_comm));
            int total_threads = 0;
            for (auto& t: each_threads) total_threads += t;
            int threads_before = 0;
            for (int j = 0; j < root.node; ++j)
            {
                threads_before += each_threads[j];
            }
            root.rank     = threads_before;
            root.num_rank = total_threads;
            for (auto& id: ids)
            {
                id.rank     = threads_before + id.thread;
                id.num_rank = total_threads;
            }
            
            data.inbox.resize(total_threads);
            data.outbox.resize(total_threads);
            
            data.rank_pids.resize(total_threads);
            
            int rank_id = 0;
            for (int node = 0; node < root.num_node; ++node)
            {
                int nthreads_here = each_threads[node];
                for (int ithread = 0; ithread < nthreads_here; ++ithread)
                {
                    proc_id_t pid_loc;
                    pid_loc.node       = node;
                    pid_loc.num_node   = root.num_node;
                    pid_loc.thread     = ithread;
                    pid_loc.num_thread = root.num_thread;
                    pid_loc.rank       = rank_id;
                    pid_loc.num_rank   = root.num_rank;
                    data.rank_pids[rank_id] = pid_loc;
                    ++rank_id;
                }
            }
        }
        
        template <typename func_t>
        void exec(const func_t& func)
        {
            const auto wrapper = [&](spade::parallel::pool_t pool)
            {
                device::set_device(devices[pool.thread()]);
                data.p2p_enabled = device::enable_p2p(pool.thread(), devices);
                
                data.inbox [pool.rank()].resize(pool.size());
                data.outbox[pool.rank()].resize(pool.size());
                
                func(pool);
            };
            for (int i = 0; i < others.size(); ++i) others[i] = std::thread(wrapper, pool_t{ids[i], context->default_comm, data});
            wrapper(pool_t{root, context->default_comm, data});
            for (auto& t: others) t.join();
        }
    };
}