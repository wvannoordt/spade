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
#include "dispatch/copy.h"
#include "dispatch/device_type.h"
#include "dispatch/device_vector.h"

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
        int node = 9999, thread = 9999, num_node = 9999, num_thread = 9999, rank = 9999, num_rank = 9999, device_id = 9999;
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
    
    enum message_channel_spec
    {
        gpu_peer_messg,
        gpu_self_messg,
        gpu_long_messg,
        cpu_self_messg,
        cpu_long_messg
    };
    
    template <typename send_t, typename recv_t>
    static void send_buffer(const proc_id_t& sender, const proc_id_t& receiver, const send_t& send_buf, recv_t& recv_buf, const message_channel_spec& channel)
    {
        if (sender.node != receiver.node)
        {
            throw except::sp_exception("send_buffer cannot be used for inter-node communication");
        }
        
        switch (channel)
        {
            case cpu_self_messg:
            {
                std::copy(&send_buf[0], &send_buf[0] + send_buf.size(), &recv_buf[0]);
                break;
            }
            case gpu_self_messg:
            {
                device::copy(send_buf, recv_buf);
                break;
            }
            case gpu_peer_messg:
            {
                device::peer_copy(send_buf, sender.device_id, recv_buf, receiver.device_id);
                break;
            }
            default:
            {
                throw except::sp_exception("fatal: invalid send_buffer channel!");
            }
        }
    }
    
    template <device::is_device dev_t>
    static message_channel_spec channel_from_device(const proc_id_t& sender, const proc_id_t& receiver, const dev_t&)
    {
        bool same_node = sender.node == receiver.node;
        if constexpr (device::is_cpu<dev_t>)                       return same_node ? cpu_self_messg : cpu_long_messg;
        if (same_node && (sender.device_id == receiver.device_id)) return gpu_self_messg;
        if (same_node && (sender.device_id != receiver.device_id)) return gpu_peer_messg;
        return gpu_long_messg;
    }
    
    //Triple vector, nice
    using mailbox_t  = std::vector<std::vector<std::vector<message_buf_t>>>;
    using channels_t = std::vector<std::vector<std::vector<message_channel_spec>>>;
    struct common_thread_data_t
    {
        mutable std::mutex mut;
        mutable barrier_t barrier;
        mutable std::vector<communicable_t> buffer;
        
        mutable mailbox_t  outbox;
        mutable mailbox_t  inbox;
        mutable channels_t mesg_channel;
        
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
        void post_send(const buffer_type& buf, int recv, const message_channel_spec& chan = cpu_self_messg) const
        {
            message_buf_t bf = buf;
            env.outbox[this->rank()][recv].push_back(bf);
            env.mesg_channel[this->rank()][recv].push_back(chan);
        }
        
        template <typename buffer_type>
        void post_recv(const buffer_type& buf, int send) const
        {
            message_buf_t bf = buf;
            env.inbox[this->rank()][send].push_back(bf);
        }
        
        template <typename buf_t>
        void send_all() const
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
                    auto& chan_list = env.mesg_channel[donor][myrank];
                    std::string msg0 = "invalid communication: list size for sender and receiver do not match";
                    std::string msg1 = "invalid communication: message size mismatch";
                    if (recv_list.size() != send_list.size()) throw except::sp_exception(msg0);
                    for (int i_msg = 0; i_msg < recv_list.size(); ++i_msg)
                    {
                        auto recv_buf = std::get<buf_t>(recv_list[i_msg]);
                        auto send_buf = std::get<buf_t>(send_list[i_msg]);
                        if (recv_buf.size() != send_buf.size()) throw except::sp_exception(msg1);
                        const auto chan = chan_list[i_msg];
                        
                        send_buffer(their_id, my_id, send_buf, recv_buf, chan);
                        
                        // if (chan == cpu_messg)     std::copy(&send_buf[0], &send_buf[0] + send_buf.size(), &recv_buf[0]);
                        // if (chan == gpu_p2p_messg)
                        // {
                        //     auto er_code = cudaMemcpy(recv_buf.ptr, send_buf.ptr, send_buf.size()*sizeof(typename buf_t::value_type), cudaMemcpyDeviceToDevice);
                        //     if ((er_code != cudaSuccess))
                        //     {
                        //         print(recv_buf.ptr);
                        //         print(send_buf.ptr);
                        //         print(send_buf.size()*sizeof(buf_t));
                        //         print(cudaGetErrorString(er_code));
                        //         print("that sucks");
                        //     }
                        // }
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
            
            for (auto& list:  inbox)                   list.clear();
            for (auto& list: outbox)                   list.clear();
            for (auto& list: env.mesg_channel[myrank]) list.clear();
            
            this->sync();
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
        
        template <typename data_t, typename binary_t>
        requires (spade::utils::variant_contains<data_t, communicable_t>)
        auto thread_reduce(const data_t& local_val, const binary_t& bin_op) const
        {
            this->post(local_val);
            data_t sum = std::get<data_t>(env.buffer[0]);
            for (int pp = 1; pp < this->num_threads(); ++pp)
            {
                sum = bin_op(sum, std::get<data_t>(env.buffer[pp]));
            }
            return sum;
        }
        
        template <typename data_t>
        data_t thread_broadcast(const int thread_id, const data_t value) const
        {
            this->post(value);
            return std::get<data_t>(env.buffer[thread_id]);
        }
        
        template <typename data_t, typename binary_t>
        auto reduce(const data_t& val, const binary_t& bin_op) const
        {
            data_t reduc_val = this->thread_reduce(val, bin_op);
            auto   dtype     = get_mpi_type(reduc_val);
            auto   base      = this->base_id();
            data_t output;
            if (this->isbase())
            {
                std::vector<data_t> alldata;
                alldata.resize(this->num_nodes());
                alldata[this->node()] = reduc_val;
                mpi_check(MPI_Allgather(&reduc_val, 1, dtype, &alldata[0], 1, dtype, this->comm));
                output = alldata[0];
                for (int pp = 1; pp < alldata.size(); ++pp) output = bin_op(output, alldata[pp]);
            } 
            return this->thread_broadcast(base, output);
        }
        
        template <typename data_t>
        auto sum(const data_t& val) const
        {
            return this->reduce(val, [](const data_t& d0, const data_t& d1) { return d0 + d1; });
        }
        
        template <device::is_device device_t>
        bool has_shared_memory(const int rank, const int dest, const device_t&)
        {
            if constexpr (device::is_gpu<device_t>) return rank == dest;
            else
            {
                auto dest_pid = env.rank_pids[dest];
                return this->pid().node == dest_pid.node;
            }
        }
        
        const proc_id_t& pid() const                { return pool_pid; }
        const proc_id_t& pid(const int other) const { return env.rank_pids[other]; }
        int    rank()      const { return this->pid().rank; }
        int    size()      const { return this->pid().num_rank; }
        int  thread()      const { return this->pid().thread; }
        int  num_nodes()   const { return this->pid().num_node; }
        int  num_threads() const { return this->pid().num_thread; }
        int  node()        const { return this->pid().node; }
        int  device_id()   const { return this->pid().device_id; }
        int  base_id()     const { return 0; }
        int  root_id()     const { return 0; }
        bool isroot()      const { return this->rank()   == this->root_id(); }
        bool isbase()      const { return this->thread() == this->base_id(); }
        bool is_root()     const { return this->isroot(); } //This is why we have proper naming conventions
        bool is_base()     const { return this->isbase(); } //This is *also* why we have proper naming conventions
        
        const pool_t& pause() const
        {
            if (this->is_root())
            {
                //Note: we can implement some kind of logging here at some point!
                std::cout << "Paused. Press \"Enter\" to continue..." << std::endl;
                std::cin.get();
            }
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
        
        template <typename func_t>
        const pool_t& exclusive(const func_t& func) const
        {
            std::lock_guard<std::mutex> guard(env.mut);
            func();
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
            channels_t(),
            std::vector<proc_id_t>(),
            false}
        {
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
            root.device_id = devices[0];
            
            data.inbox.resize(total_threads);
            data.outbox.resize(total_threads);
            data.mesg_channel.resize(total_threads);
            
            data.rank_pids.resize(total_threads);
            
            int idx = 1;
            for (auto& id: ids)
            {
                id.device_id = devices[idx];
                ++idx;
            }
            
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
                    pid_loc.device_id  = devices[ithread];
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
                device::set_device(pool.device_id());
                
                data.p2p_enabled = device::enable_p2p(pool.thread(), devices);
                
                data.inbox       [pool.rank()].resize(pool.size());
                data.outbox      [pool.rank()].resize(pool.size());
                data.mesg_channel[pool.rank()].resize(pool.size());
                
                func(pool);
            };
            for (int i = 0; i < others.size(); ++i) others[i] = std::thread(wrapper, pool_t{ids[i], context->default_comm, data});
            wrapper(pool_t{root, context->default_comm, data});
            for (auto& t: others) t.join();
        }
    };
}