#pragma once

#include "core/print.h"
#include "core/ctrs.h"
#include "core/parallel.h"
#include "core/tag_val.h"

#include "core/except.h"

namespace spade::partition
{
    static struct global_tag_t{} global;
    static struct local_tag_t{}  local;
    
    template <typename T> concept partition_tagged = utils::has_tag<T, global_tag_t> || utils::has_tag<T, local_tag_t>;
    class block_partition_t
    {
        public:
            constexpr static std::size_t no_value = std::numeric_limits<std::size_t>::max();
            block_partition_t()
            {
                
            }
            
            block_partition_t(const auto& block_config_in, const parallel::mpi_t& group_in)
            {
                num_global_blocks = block_config_in.total_num_blocks();
                num_local_blocks = num_global_blocks/group_in.size();
                if (group_in.rank() < (num_global_blocks%group_in.size())) num_local_blocks++;
                local_block_to_global_block.resize(num_local_blocks);

                
                global_block_to_local_block.resize(num_global_blocks, no_value);
                global_block_to_rank.resize(num_global_blocks);
                std::size_t current_rank = 0;
                std::size_t loc_count = 0;
                for (auto lb: range(0, num_global_blocks))
                {
                    global_block_to_rank[lb] = current_rank;
                    if (current_rank == group_in.rank())
                    {
                        local_block_to_global_block[loc_count] = lb;
                        global_block_to_local_block[lb] = loc_count;
                        ++loc_count;
                    }
                    ++current_rank;
                    current_rank %= group_in.size();
                }
                if (num_global_blocks != group_in.sum(num_local_blocks))
                    throw except::sp_exception("grid partition failed consistency check");
            }
            
            std::size_t get_num_global_blocks() const { return num_global_blocks; }
            std::size_t get_num_local_blocks () const { return num_local_blocks; }
            std::size_t get_global_block(const std::size_t& lb_loc)   const { return local_block_to_global_block[lb_loc]; }
            std::size_t get_local_block (const std::size_t& lb_glob)  const { return global_block_to_local_block[lb_glob]; }
            std::size_t get_global_rank (const std::size_t& lb_glob)  const { return global_block_to_rank[lb_glob]; }
            
            std::size_t get_rank (const partition_tagged auto& lb)  const { return global_block_to_rank[to_global(lb).value]; }
            
            template <partition_tagged idx_t> auto to_global(const idx_t& lb) const
            {
                if constexpr (utils::has_tag<idx_t, global_tag_t>) return lb;
                else return utils::tag[global](local_block_to_global_block[lb.value]);
            }
            
            // template <std::integral integ_t> auto to_local(const utils::tagged_value_t<global_tag_t, integ_t>& lb_glob) const
            // template <typename integ_t> auto to_local(const utils::tagged_value_t<global_tag_t, integ_t>& lb_glob) const
            template <partition_tagged idx_t> auto to_local(const idx_t& lb) const
            {
                if constexpr (utils::has_tag<idx_t, local_tag_t>) return lb;
                else return utils::tag[local](global_block_to_local_block[lb.value]);
            }
            
        private:
            std::size_t num_local_blocks, num_global_blocks;
            std::vector<std::size_t> local_block_to_global_block;
            std::vector<std::size_t> global_block_to_local_block;
            std::vector<std::size_t> global_block_to_rank;
    };
}