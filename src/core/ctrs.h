#pragma once

#include <type_traits>
#include <concepts>
#include "core/utils.h"
namespace cvdf::ctrs
{
    template <class T> concept basic_array = requires(T t, const int& i)
    {
        t.size();
        t[i];
    };
    
    template <class T> concept not_basic_array = !basic_array<T>;
    
    template <class T, typename data_t> concept basic_array_of_type = basic_array<T> && requires(T t, const int& i)
    {
        {t[i]}-> std::convertible_to<data_t>;
    };
    
    template <class T, const size_t size, typename real_type> concept vec_nd = basic_array<T>
    && (sizeof(T) == size*sizeof(real_type))
    && requires(T t, size_t i)
    {
        { t[i] } -> std::common_with<real_type>;
        typename T::value_type;
    };
    
    template <class T> concept integral_type = std::is_integral<T>::value;


    template <basic_array a1_t, basic_array a2_t> static void copy_array(const a1_t& src, a2_t& dest)
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
        template <typename ftype> void fill(const ftype& val)
        {
            for (size_t i = 0; i < this->size(); i++) data[i] = val;
        }
        template <not_basic_array param> void set_r(const size_t i, const param& p)
        {
            data[i] = p;
        }
        template <not_basic_array param, class... params> void set_r(const size_t i, const param& p, params... ps)
        {
            data[i] = p;
            set_r(i+1, ps...);
        }
        template <basic_array param, class... params> void set_r(const size_t i, const param& p)
        {
            for (std::size_t i = 0; i < p.size(); i++) data[i] = p[i];
        }
        template <class... params> array(params... ps)
        {
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
        
        template <basic_array arr_t> auto& operator += (const arr_t& rhs)
        {
            for (std::size_t i = 0; i < this->size(); i++) data[i] += rhs[i];
            return *this;
        }
        
        template <basic_array arr_t> auto& operator -= (const arr_t& rhs)
        {
            for (std::size_t i = 0; i < this->size(); i++) data[i] -= rhs[i];
            return *this;
        }
        
        template <typename data_t> auto& operator *= (const data_t& rhs)
        {
            for (std::size_t i = 0; i < this->size(); i++) data[i] *= rhs;
            return *this;
        }
        
        template <typename data_t> auto& operator /= (const data_t& rhs)
        {
            for (std::size_t i = 0; i < this->size(); i++) data[i] /= rhs;
            return *this;
        }
        
        auto& operator += (const dtype& rhs)
        {
            for (std::size_t i = 0; i < this->size(); i++) data[i] += rhs;
            return *this;
        }
        
        auto& operator -= (const dtype& rhs)
        {
            for (std::size_t i = 0; i < this->size(); i++) data[i] -= rhs;
            return *this;
        }
        
        template <basic_array arr_t> auto& operator %= (const arr_t& rhs)
        {
            for (std::size_t i = 0; i < this->size(); i++) data[i] %= rhs[i];
            return *this;
        }
        
        bool operator == (const array<dtype,ar_size>& rhs) const
        {
            bool output = true;
            for (std::size_t i = 0; i < this->size(); i++) output = (output&&(data[i]==rhs[i]));
            return output;
        }
    };
    
    template <typename dtype, const size_t ar_size> static std::ostream & operator<<(std::ostream & os, const array<dtype, ar_size>& arr)
    {
       os << "[";
       for (size_t i = 0; i < arr.size(); i++)
       {
           os << arr.data[i];
           if (i< arr.size()-1) os << ", ";
       }
       os << "]";
       return os;
    }
    
    template <basic_array arr_l_t, basic_array arr_r_t> static auto collapse_index(const arr_l_t& idx, const arr_r_t& dims)
    {
        static_assert(std::is_integral<typename arr_l_t::value_type>::value, "collapse_index requires integral arrays");
        static_assert(std::is_integral<typename arr_r_t::value_type>::value, "collapse_index requires integral arrays");
        typename arr_l_t::value_type coeff = 1;
        typename arr_l_t::value_type output = 0;
        for (auto i: range(0, dims.size()))
        {
            output += coeff*idx[i];
            coeff*=dims[i];
        }
        return output;
    }
    
    template <typename idx_t, basic_array arr_t> static auto expand_index(const idx_t& idx, const arr_t& dims)
    {
        static_assert(std::is_integral<typename arr_t::value_type>::value, "collapse_index requires integral arrays");
        typename arr_t::value_type coeff = 1;
        typename arr_t::value_type modcoeff = 1;
        arr_t output;
        idx_t temp_idx = idx;
        idx_t sum = 0;
        for (auto i: range(0, dims.size()))
        {
            modcoeff *= dims[i];
            output[i] = ((temp_idx-sum) % modcoeff) / coeff;
            sum += output[i]*coeff;
            coeff *= dims[i];
        }
        return output;
    }
    
    template <typename data_t, const std::size_t ar_size> struct static_stack
    {
        data_t data[ar_size];
        int next = 0;
        void push(const data_t& new_data)
        {
            data[next++] = new_data;
            next -= (next==ar_size)*ar_size;
        }
        const data_t& operator[] (const uint i) const
        {
            return data[(i+next)%ar_size];
        }
        data_t& operator[] (const uint i)
        {
            return data[(i+next)%ar_size];
        }
    };
}