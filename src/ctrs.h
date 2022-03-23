#pragma once

#include <type_traits>
#include <concepts>
#include "utils.h"
namespace cvdf::ctrs
{
    template <class T> concept basic_array = requires(T t, const int& i)
    {
        t.size();
        t[i];
    };
    
    template <class T, const size_t size, typename real_type> concept vec_nd = (sizeof(T) == size*sizeof(real_type)) && requires(T t, size_t i)
    {
        basic_array<T>;
        { t[i] } -> std::common_with<real_type>;
        typename T::value_type;
    };
    
    template <class T> concept integral_type = std::is_integral<T>::value;


    template <basic_array a1_t, basic_array a2_t> void copy_array(const a1_t& src, a2_t& dest)
    {
        std::size_t tsize = utils::min(src.size(), dest.size());
        for (std::size_t i = 0; i < tsize; ++i) dest[i] = src[i];
    }
    
    
    template<typename dtype, const size_t ar_size> struct array
    {
        typedef dtype value_type;
        dtype data[ar_size];
        dtype& operator [] (size_t idx) {return data[idx];}
        dtype* begin() noexcept {return &data[0];}
        dtype* end()   noexcept {return &data[0]+ar_size;}
        constexpr size_t size(void) const noexcept {return ar_size;}
        void fill(const dtype& val)
        {
            for (size_t i = 0; i < this->size(); i++) data[i] = val;
        }
        template <class param> void set_r(const size_t i, const param& p)
        {
            data[i] = p;
        }
        template <class param, class... params> void set_r(const size_t i, const param& p, params... ps)
        {
            data[i] = p;
            set_r(i+1, ps...);
        }
        template <class... params> array(params... ps)
        {
            static_assert(sizeof...(ps)==ar_size, "Incorrect number of arguments passed to array constructor");
            set_r(0, ps...);
        }
        array(const dtype& val) {fill(val);}
        array(void) {fill(dtype());}
        template <typename rhs_t> array<dtype, ar_size>& operator = (const rhs_t& rhs)
        {
            this->fill(rhs);
            return *this;
        }
        template <integral_type idx_t> const dtype& operator[] (const idx_t& idx) const noexcept { return data[idx]; }
        template <integral_type idx_t>       dtype& operator[] (const idx_t& idx)       noexcept { return data[idx]; }
    };
    
    template <basic_array arr_l_t, basic_array arr_r_t> auto collapse_index(const arr_l_t& idx, const arr_r_t& dims)
    {
        static_assert(std::is_integral<typename arr_l_t::value_type>::value, "collapse_index requires integral arrays");
        static_assert(std::is_integral<typename arr_r_t::value_type>::value, "collapse_index requires integral arrays");
        typename arr_l_t::value_type coeff = 1;
        typename arr_l_t::value_type output = 0;
        for (auto i: range(0, dims.size()))
        {
            output += coeff*idx[i[0]];
            coeff*=dims[i[0]];
        }
        return output;
    }
    
    template <typename idx_t, basic_array arr_t> auto expand_index(const idx_t& idx, const arr_t& dims)
    {
        static_assert(std::is_integral<typename arr_t::value_type>::value, "collapse_index requires integral arrays");
        typename arr_t::value_type coeff = 1;
        typename arr_t::value_type modcoeff = 1;
        arr_t output;
        idx_t temp_idx = idx;
        idx_t sum = 0;
        for (auto i: range(0, dims.size()))
        {
            modcoeff *= dims[i[0]];
            output[i[0]] = ((temp_idx-sum) % modcoeff) / coeff;
            sum += output[i[0]]*coeff;
            coeff *= dims[i[0]];
        }
        return output;
    }
}