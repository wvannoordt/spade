#pragma once

#include <type_traits>
#include <cmath>
#include <utility>
#include <initializer_list>

#include "core/config.h"
#include "core/utils.h"
#include "core/range.h"
#include "core/const_index.h"
#include "core/unit_vector.h"

namespace spade::ctrs
{
    template <class T> concept basic_array = requires(T t, const int& i)
    {
        typename T::value_type;
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
    
    template <typename left_t, typename right_t> concept convertible_array =
    basic_array<left_t> && basic_array<right_t> &&
    std::convertible_to<typename left_t::value_type, typename right_t::value_type>;
    
    template <typename left_t, typename right_t> concept same_size = (left_t::size() == right_t::size());

    template <basic_array arr_t>
    requires(std::floating_point<typename arr_t::value_type>)
    _sp_hybrid constexpr static auto array_norm(const arr_t& arr)
    {
        using data_t = typename arr_t::value_type;
        data_t output = data_t();
        for (int i = 0; i < arr.size(); ++i) output += arr[i]*arr[i];
        return sqrt(output);
    }
    
    template <basic_array arr_t>
    requires (std::floating_point<typename arr_t::value_type> && arr_t::size() == 3)
    _sp_hybrid constexpr static arr_t cross_prod(const arr_t& v0, const arr_t& v1)
    {
        return {v0[1]*v1[2] - v0[2]*v1[1], v0[2]*v1[0] - v0[0]*v1[2], v0[0]*v1[1] - v0[1]*v1[0]};
    }
    
    template <basic_array arr1_t, basic_array arr2_t>
    requires (arr1_t::size() == arr2_t::size())
    _sp_hybrid constexpr static auto dot_prod(const arr1_t& arr1, const arr2_t& arr2)
    {
        using f1_t = typename arr1_t::value_type;
        using f2_t = typename arr2_t::value_type;
        using data_t = decltype(f1_t()*f2_t());
        data_t output = arr1[0]*arr2[0];
        for (int i = 1; i < arr1.size(); ++i) output += arr1[i]*arr2[i];
        return output;
    }
    
    template <basic_array a1_t, basic_array a2_t>
    _sp_hybrid constexpr static void copy_array(const a1_t& src, a2_t& dest, const typename a2_t::value_type& default_val)
    {
        std::size_t tsize = utils::min(src.size(), dest.size());
        for (std::size_t i = 0; i < dest.size(); ++i) dest[i] = default_val;
        for (std::size_t i = 0; i < tsize;       ++i) dest[i] = src[i];
    }
    
    template <basic_array a1_t, basic_array a2_t>
    requires(a1_t::size() == a2_t::size())
    _sp_hybrid constexpr static void copy_array(const a1_t& src, a2_t& dest)
    {
        for (std::size_t i = 0; i < a2_t::size(); ++i) dest[i] = src[i];
    }
    
    template <basic_array arr_t, typename rhs_t>
    _sp_hybrid static void fill_array(arr_t& arr, const rhs_t& val)
    {
        for (std::size_t i = 0; i < arr.size(); ++i) arr[i] = val;
    }
    
    namespace detail{struct no_type_t{};}
    
    template<typename dtype, const std::size_t ar_size, typename derived_t = detail::no_type_t>
    struct arithmetic_array_t
    {
        using self_type = std::conditional<
            std::is_same<derived_t, detail::no_type_t>::value,
            arithmetic_array_t<dtype, ar_size, derived_t>,
            derived_t>::type;
        
        using value_type = dtype;
        using index_type = int;
        
        constexpr static bool has_derived = std::same_as<derived_t, detail::no_type_t>;
        
        dtype data[ar_size];
        _sp_hybrid dtype& operator [] (index_type idx) {return data[idx];}
        _sp_hybrid const dtype& operator [] (index_type idx) const {return data[idx];}
        
        template <udci::integral_t ii>
        requires(ii < ar_size)
        _sp_hybrid dtype& operator[] (const udci::idx_const_t<ii>& idx) {return data[ii];}
        
        template <udci::integral_t ii>
        requires(ii < ar_size)
        _sp_hybrid const dtype& operator[] (const udci::idx_const_t<ii>& idx) const {return data[ii];}
        
        template <udci::integral_t ii>
        _sp_hybrid const dtype& operator[] (const udci::idx_const_t<ii>& idx) const
        {
            static_assert(ii<ar_size, "constant literal index must be less than array size"); return data[ii];
        }
        
        _sp_hybrid dtype* begin() noexcept {return &data[0];}
        _sp_hybrid dtype* end()   noexcept {return &data[0]+ar_size;}

        _sp_hybrid const dtype* begin() const noexcept {return &data[0];}
        _sp_hybrid const dtype* end()   const noexcept {return &data[0]+ar_size;}

        _sp_hybrid constexpr static index_type size() {return ar_size;}
        
        _sp_hybrid self_type& self() { return static_cast<self_type&>(*this); }
        
        template <typename ftype> _sp_hybrid void fill(const ftype& val, const index_type imin, const index_type imax)
        {
            for (index_type i = imin; i < imax; i++) data[i] = val;
        }
        
        template <typename ftype> _sp_hybrid void fill(const ftype& val)
        {
            fill(val, 0, this->size());
        }
        
        template <not_basic_array param> _sp_hybrid  void set_r(const index_type i, const param& p)
        {
            data[i] = p;
        }
        template <not_basic_array param, class... params> _sp_hybrid  void set_r(const index_type i, const param& p, params... ps)
        {
            data[i] = p;
            set_r(i+1, ps...);
        }
        
        template <not_basic_array... params>
        requires(sizeof...(params) == ar_size)
        _sp_hybrid  arithmetic_array_t(params... ps)
        {
            set_r(0, ps...);
        }
        
        _sp_hybrid arithmetic_array_t(std::initializer_list<dtype>& llist)
        {
            std::copy(this->begin(), this->end(), llist.begin());
        }
        
       _sp_hybrid arithmetic_array_t(const dtype& val) {fill(val);}
        
        _sp_hybrid arithmetic_array_t() {}
        _sp_hybrid constexpr arithmetic_array_t(const arithmetic_array_t& rhs) = default;
        
        template <typename rhs_t>
        requires (!basic_array<rhs_t>)
        _sp_hybrid auto& operator = (const rhs_t& rhs)
        {
            this->fill(rhs);
            return self();
        }
        
        template <basic_array rhs_t>
        _sp_hybrid constexpr self_type& operator = (const rhs_t& rhs)
        {
            copy_array(rhs, *this);
            return self();
        }
        
        _sp_hybrid constexpr self_type& operator = (const dtype& rhs)
        {
            this->fill(rhs);
            return self();
        }
        
        template <integral_type idx_t> _sp_hybrid const dtype& operator[] (const idx_t& idx) const noexcept { return data[idx]; }
        template <integral_type idx_t> _sp_hybrid       dtype& operator[] (const idx_t& idx)       noexcept { return data[idx]; }
        
        template <basic_array arr_t>
        _sp_hybrid self_type& operator += (const arr_t& rhs)
        {
            #pragma unroll
            for (index_type i = 0; i < this->size(); i++) data[i] += rhs[i];
            return self();
        }
        
        template <typename rhs_t>
        _sp_hybrid auto operator + (const arithmetic_array_t<rhs_t, ar_size, derived_t>& rhs) const
        {
            arithmetic_array_t<decltype(dtype()+rhs_t()), ar_size, derived_t> output;
            #pragma unroll
            for (index_type i = 0; i < this->size(); i++) output[i] = data[i] + rhs[i];
            return output.self();
        }
        
        template <typename rhs_t>
        _sp_hybrid auto operator - (const arithmetic_array_t<rhs_t, ar_size, derived_t>& rhs) const
        {
            arithmetic_array_t<decltype(dtype()+rhs_t()), ar_size, derived_t> output;
            #pragma unroll
            for (index_type i = 0; i < this->size(); i++) output[i] = data[i] - rhs[i];
            return output.self();
        }
        
        template <basic_array arr_t>
        _sp_hybrid self_type& operator -= (const arr_t& rhs)
        {
            #pragma unroll
            for (index_type i = 0; i < this->size(); i++) data[i] -= rhs[i];
            return self();
        }
        
        template <typename data_t>
        _sp_hybrid self_type& operator *= (const data_t& rhs)
        {
            #pragma unroll
            for (index_type i = 0; i < this->size(); i++) data[i] *= rhs;
            return self();
        }
        
        template <typename data_t>
        _sp_hybrid self_type& operator /= (const data_t& rhs)
        {
            #pragma unroll
            for (index_type i = 0; i < this->size(); i++) data[i] /= rhs;
            return self();
        }
        
        _sp_hybrid auto& operator += (const dtype& rhs)
        {
            #pragma unroll
            for (index_type i = 0; i < this->size(); i++) data[i] += rhs;
            return self();
        }
        
        _sp_hybrid auto& operator -= (const dtype& rhs)
        {
            #pragma unroll
            for (index_type i = 0; i < this->size(); i++) data[i] -= rhs;
            return self();
        }
        
        template <basic_array arr_t> auto& operator %= (const arr_t& rhs)
        {
            #pragma unroll
            for (index_type i = 0; i < this->size(); i++) data[i] %= rhs[i];
            return self();
        }
        
        _sp_hybrid bool operator == (const arithmetic_array_t<dtype, ar_size, derived_t>& rhs) const
        {
            bool output = true;
            #pragma unroll
            for (index_type i = 0; i < this->size(); i++) output = (output && (data[i]==rhs[i]));
            return output;
        }
    };
    
    template<typename scal_type, typename dtype, const std::size_t ar_size, typename derived_t>
    _sp_hybrid auto operator*(const scal_type& lhs, const arithmetic_array_t<dtype, ar_size, derived_t>& rhs)
    {
        using out_t = dtype;
        arithmetic_array_t<out_t, ar_size, derived_t> output;
        #pragma unroll
        for (int i = 0; i < rhs.size(); ++i) output[i] = (rhs[i] * lhs);
        return output.self();
    }
    
    template<typename div_t, typename dtype, const std::size_t ar_size, typename derived_t>
    _sp_hybrid auto operator/(const arithmetic_array_t<dtype, ar_size, derived_t>& rhs, const div_t& lhs)
    {
        using out_t = dtype;
        arithmetic_array_t<out_t, ar_size, derived_t> output;
        #pragma unroll
        for (int i = 0; i < rhs.size(); ++i) output[i] = (rhs[i] / lhs);
        return output.self();
    }
    
    template<typename dtype, const std::size_t ar_size, typename derived_t>
    _sp_hybrid auto operator+(const dtype& lhs, const arithmetic_array_t<dtype, ar_size, derived_t>& rhs)
    {
        arithmetic_array_t<dtype, ar_size, derived_t> output = rhs;
        for (auto& i: output) i+= lhs;
        return output.self();
    }
    
    template <typename data_t, const std::size_t ar_size> using array = arithmetic_array_t<data_t, ar_size>;
    template <typename data_t> using v3d = array<data_t,3>;
    
    
    template <basic_array rhs_t>
    _sp_hybrid constexpr static array<typename rhs_t::value_type, rhs_t::size()> to_array(const rhs_t& rhs)
    {
        array<typename rhs_t::value_type, rhs_t::size()> output;
        copy_array(rhs, output);
        return output;
    }
    
    template <basic_array a1_t, basic_array a2_t>
    _sp_hybrid constexpr static void copy_array(const a1_t& src, a2_t& dest)
    {
        std::size_t tsize = utils::min(src.size(), dest.size());
        for (std::size_t i = 0; i < tsize; ++i) dest[i] = src[i];
    }
    
    template <typename dtype, const size_t ar_size, typename derived_t>
    static std::ostream &
    operator<<(std::ostream & os, const arithmetic_array_t<dtype, ar_size, derived_t>& arr)
    {
       os << "[";
       for (size_t i = 0; i < arr.size(); i++)
       {
           os << arr.data[i];
           if (i < arr.size()-1) os << ", ";
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
    
    //only for use with external data-owning applications
    template <typename data_t> struct unsafe_vector_alias_t
    {
        typedef data_t value_type;
        data_t* base;
        std::size_t safe_size;
        const data_t& operator [] (const std::size_t& idx) const {return base[idx];}
        data_t& operator [] (const std::size_t& idx) {return base[idx];}
        std::size_t size() const {return safe_size;}
        unsafe_vector_alias_t(){}
        unsafe_vector_alias_t(data_t* base_in, const std::size_t& safe_size_in){base = base_in; safe_size = safe_size_in;}
    };
    
    template <typename target_t, typename data_t, const std::size_t ar_size, typename derived_t>
    _sp_hybrid inline target_t array_cast(const arithmetic_array_t<data_t, ar_size, derived_t>& arr)
    {
        //Kind of stupid but quite useful.
        target_t output;
        copy_array(arr, output);
        return output;
    }
    
    template <typename val_t, typename data_t, const std::size_t ar_size, typename derived_t>
    _sp_hybrid inline auto array_val_cast(const arithmetic_array_t<data_t, ar_size, derived_t>& arr)
    {
        //Kind of stupid but quite useful.
        arithmetic_array_t<val_t, ar_size, derived_t> output;
        copy_array(arr, output);
        return output;
    }
    
    template <typename... vals_t>
    _sp_hybrid inline auto make_array(const vals_t&... vals)
    {
        using val_t    = typename std::common_type<vals_t...>::type;
        using output_t = array<val_t, sizeof...(vals_t)>;
        return output_t{vals...};
    }
}