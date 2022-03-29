#pragma once

#include "core/print.h"
#include "core/ctrs.h"
#include "core/parallel.h"

namespace cvdf::partition
{
    class block_partition_t
    {
        public:
            block_partition_t(void)
            {
                
            }
            
            block_partition_t(const ctrs::array<std::size_t, 3>& num_blocks_in, parallel::mpi_t* group_in)
            {
                ctrs::copy_array(num_blocks_in, num_blocks);
                std::size_t total_blocks = num_blocks[0]*num_blocks[1]*num_blocks[2];
                num_local_blocks = total_blocks/group_in->size();
                if (group_in->rank() < (total_blocks%group_in->size())) num_local_blocks++;
                local_block_to_global_block.resize(num_local_blocks);
                global_block_to_local_block.resize(total_blocks,-1);
                global_block_to_rank.resize(total_blocks);
                std::size_t current_rank = 0;
                std::size_t loc_count = 0;
                for (auto lb: range(0, total_blocks))
                {
                    global_block_to_rank[lb[0]] = current_rank;
                    if (current_rank == group_in->rank())
                    {
                        local_block_to_global_block[loc_count] = lb[0];
                        global_block_to_local_block[lb[0]] = loc_count;
                        ++loc_count;
                    }
                    ++current_rank;
                    current_rank %= group_in->size();
                }
            }
            
            std::size_t get_num_local_blocks(void) const { return num_local_blocks; }
            std::size_t get_global_block(const std::size_t& lb_loc)   const { return local_block_to_global_block[lb_loc]; }
            std::size_t get_local_block (const std::size_t& lb_glob)  const { return global_block_to_local_block[lb_glob]; }
            
            std::size_t get_global_rank (const std::size_t& lb_glob) const { return global_block_to_rank[lb_glob]; }
            
        private:
            ctrs::array<std::size_t, 3> num_blocks;
            std::size_t num_local_blocks;
            std::vector<std::size_t> local_block_to_global_block;
            std::vector<std::size_t> global_block_to_local_block;
            std::vector<std::size_t> global_block_to_rank;
    };
}