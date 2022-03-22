#pragma once

#include "print.h"
#include "ctrs.h"
#include "parallel.h"

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
            }
        private:
            ctrs::array<std::size_t, 3> num_blocks;
            std::size_t num_local_blocks;
    };
}