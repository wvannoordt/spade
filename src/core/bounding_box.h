#pragma once

#include <iostream>
#include <concepts>
#include "core/range.h"

namespace spade
{
    template <typename dtype, const std::size_t ar_size> struct bound_box_t
    {
        ctrs::array<dtype, 2*ar_size> bnds;
        
        using index_type = ctrs::array<dtype, ar_size>;
        
        _sp_hybrid const auto& data() const {return bnds;}
        
        _sp_hybrid bound_box_t(){}
        _sp_hybrid bound_box_t(const dtype& v)
        {
            for (int i = 0; i < 2*ar_size; ++i) bnds[i] = v;
        }
        
        template <typename... rs_t>
        requires ((sizeof...(rs_t) == 2*ar_size) && (std::same_as<rs_t, dtype> && ...))
        _sp_hybrid bound_box_t(const rs_t&... rs)
        {
            bnds = ctrs::array<dtype, 2*ar_size>(rs...);
        }
        
        _sp_hybrid bound_box_t inflate(const dtype& infl) const
        {
            bound_box_t output = *this;
            for (int d = 0; d < ar_size; ++d)
            {
                dtype hh = dtype(0.5)*this->size(d);
                dtype cc = dtype(0.5)*(this->min(d) + this->max(d));
                output.min(d) = cc - infl*hh;
                output.max(d) = cc + infl*hh;
            }
            return output;
        }
        
        template <ctrs::basic_array arr_t>
        requires(arr_t::size() == ar_size)
        _sp_hybrid bool contains(const arr_t& x) const
        {
            for (int d = 0; d < ar_size; ++d)
            {
                if ((x[d]<min(d)) || (x[d]>max(d))) return false;
            }
            return true;
        }
        
        _sp_hybrid bool intersects(const bound_box_t& rhs) const
        {
            for (int i = 0; i < ar_size; ++i)
            {
                if (!((this->min(i) <= rhs.max(i)) && (this->max(i) >= rhs.min(i)))) return false;
            }
            return true;
        }
        
        _sp_hybrid const dtype& min(size_t idx) const {return bnds[2*idx+0];}
        _sp_hybrid const dtype& max(size_t idx) const {return bnds[2*idx+1];}
        
        _sp_hybrid dtype& min(size_t idx) {return bnds[2*idx+0];}
        _sp_hybrid dtype& max(size_t idx) {return bnds[2*idx+1];}
        
        _sp_hybrid dtype size(size_t idx) const {return max(idx)-min(idx);}
        _sp_hybrid dtype volume() const
        {
            dtype output = 1;
            for (std::size_t i = 0; i < ar_size; i++)
            {
                output *= this->size(i);
            }
            return output;
        }
        
        _sp_hybrid dtype& operator() (const std::size_t& dim, const std::size_t& min_max) {return bnds[2*dim+min_max];}
        _sp_hybrid const dtype& operator() (const std::size_t& dim, const std::size_t& min_max) const {return bnds[2*dim+min_max];}

        _sp_hybrid bound_box_t operator ! () const
        {
            static_assert(std::same_as<dtype, bool>, "cannot apply unary operator ! to non-boolean bound box");
            bound_box_t output = (*this);
            for (auto& b: output) b = !b;
            return output;
        }
    };
    
    template <typename dtype, const size_t ar_size> static std::ostream & operator<<(std::ostream & os, const bound_box_t<dtype, ar_size> & pos)
    {
        os << "{";
        for (auto i: range(0, ar_size))
        {
            os << "(" << pos.min(i) << ", " << pos.max(i) << ")";
            if (i < ar_size-1) os << ", ";
        }
        os << "}";
        return os;
    }
}