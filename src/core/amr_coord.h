#pragma once

#include <cstdint>

namespace spade::amr
{
    template <typename partition_t, typename bit_repr_t> struct amr_exact_coord_t
    {
        using int_t = int;
        bit_repr_t  bits;
        partition_t partition, num_partitions;
        
        amr_exact_coord_t()
        {
            partition = 0;
            bits = 0;
        }
        
        amr_exact_coord_t(const partition_t& partition_in, const partition_t& num_partitions_in, const bit_repr_t& bits_in)
        {
            partition = partition_in;
            num_partitions = num_partitions_in;
            bits = bits_in;
        }
        
        amr_exact_coord_t(const amr_exact_coord_t& rhs)
        {
            num_partitions = rhs.num_partitions;
            partition = rhs.partition;
            bits = rhs.bits;
        }
        
        static constexpr std::size_t get_byte_size()
        {
            return 8;
        }
        
        static constexpr std::size_t get_num_bits()
        {
            return sizeof(bit_repr_t)*get_byte_size();
        }
        
        char get_bit(const std::size_t& position) const
        {
            auto offset = get_num_bits() - 1 - position;
            bit_repr_t unit = 1;
            bit_repr_t mask = (unit<<offset);
            return (char)((bits&mask)>>offset);
        }
        
        void set_bit(const std::size_t& position, const int_t& value)
        {
            auto offset = get_num_bits() - 1 - position;
            bit_repr_t unit = 1;
            bit_repr_t mask = (unit<<offset);
            bit_repr_t bval = value;
            bit_repr_t val = ((bval&1)<<offset);
            bits = (bits&~mask)|(val&mask);
        }
        
        bit_repr_t get_interval_size(const int_t& level) const
        {
            auto offset = get_num_bits() - level;
            bit_repr_t unit = 1;
            return (unit<<offset);
        }
        
        amr_exact_coord_t& operator =(const amr_exact_coord_t& rhs)
        {
            bits = rhs.bits;
            partition = rhs.partition;
            num_partitions = rhs.num_partitions;
            return *this;
        }
        
        bool operator == (const amr_exact_coord_t & rhs) const
        {
            return ((bits==rhs.bits)&&(partition==rhs.partition));
        }
        
        bool operator < (const amr_exact_coord_t & rhs) const
        {
            if (partition<rhs.partition) return true;
            if (partition>rhs.partition) return false;
            return bits<rhs.bits;
        }
        
        bool operator > (const amr_exact_coord_t & rhs) const
        {
            if (partition<rhs.partition) return false;
            if (partition>rhs.partition) return true;
            return bits>rhs.bits;
        }
        
        bool operator <= (const amr_exact_coord_t& rhs) const
        {
            if (partition<rhs.partition) return true;
            if (partition>rhs.partition) return false;
            return bits<=rhs.bits;
        }
        
        bool operator >= (const amr_exact_coord_t & rhs) const
        {
            if (partition<rhs.partition) return false;
            if (partition>rhs.partition) return true;
            return bits>=rhs.bits;
        }

        template <typename base_coord_t>
        base_coord_t convert_to_coordinate(const base_coord_t& grid_min, const base_coord_t& partition_size) const
        {
            base_coord_t output = grid_min + partition*partition_size;
            base_coord_t dx = 0.5*partition_size;
            for (auto i: range(0,get_num_bits()))
            {
                output += dx*get_bit(i);
                dx *= 0.5;
            }
            return output;
        }
    };
    
    template <typename T1, typename T2> static std::ostream & operator<<(std::ostream & os, const amr_exact_coord_t<T1, T2> & pos)
    {
        os << "P:" << pos.partition << "[";
        T2 mask = 1;
        for (int i = 0; i < pos.get_num_bits(); i++)
        {
            os<<(int)(pos.get_bit(i));
        }
        os << "]";
        return os;
    }
    
    using amr_coord_t = amr_exact_coord_t<int,uint64_t>;
}