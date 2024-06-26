#pragma once

#include <vector>
#include "core/cuda_incl.h"
#include "core/ctrs.h"

namespace spade::utils
{
    template <typename data_t>
    struct vec_image_t
    {
        using value_type = data_t;
        
        data_t* ptr;
        std::size_t csize;
        
        _sp_hybrid std::size_t size() const { return csize; }
        _sp_hybrid data_t& back() { return ptr[csize-1]; }
        _sp_hybrid data_t& operator [] (const std::size_t& idx) { return ptr[idx]; }
        _sp_hybrid const data_t& operator [] (const std::size_t& idx) const { return ptr[idx]; }
        
        vec_image_t& swap (vec_image_t& rhs)
        {
            std::swap(rhs.ptr,   ptr);
            std::swap(rhs.csize, csize);
            return *this;
        }
    };
    
    template <typename data_t>
    struct const_vec_image_t
    {
        using value_type = data_t;
        
        const data_t* ptr;
        std::size_t csize;
        
        _sp_hybrid std::size_t size() const { return csize; }
        _sp_hybrid const data_t& back() const { return ptr[csize-1]; }
        _sp_hybrid const data_t& operator [] (const std::size_t& idx) const { return ptr[idx]; }
        
        const_vec_image_t& swap (const_vec_image_t& rhs)
        {
            std::swap(rhs.ptr,   ptr);
            std::swap(rhs.csize, csize);
            return *this;
        }
    };
    
    namespace detail
    {
        template <const int idx, const std::size_t rank, std::integral coeff_t, std::integral i_t>
        _sp_hybrid inline coeff_t r_md_offst(const ctrs::array<coeff_t, rank>& sizes, const i_t& i)
        {
            return i;
        }
        
        template <const int idx, const std::size_t rank, std::integral coeff_t, std::integral i_t, std::integral... is_t>
        requires(sizeof...(is_t) > 0) // stupid CUDA compiler
        _sp_hybrid inline coeff_t r_md_offst(const ctrs::array<coeff_t, rank>& sizes, const i_t& i, const is_t&... is)
        {
            return i + sizes[idx]*(r_md_offst<idx+1>(sizes, is...));
        }
    }
    
    template <typename data_t, const std::size_t rank, std::integral coeff_t>
    struct md_vec_image_t
    {
        using value_type = data_t;
        
        data_t* ptr;
        ctrs::array<coeff_t, rank> sizes;
        std::size_t csize;
        
        _sp_hybrid std::size_t size() const { return csize; }
        
        template <std::integral... is_t>
        requires (sizeof...(is_t) == rank)
        _sp_hybrid data_t& operator () (const is_t&... is)
        {
            coeff_t offst = detail::r_md_offst<0>(sizes, is...);
            return ptr[offst];
        }
        
        _sp_hybrid data_t* end()
        {
            std::size_t offs = 1;
            for (int i = 0; i < rank; ++i) offs *= sizes[i];
            return ptr + offs;
        }
        
        template <std::integral... is_t>
        requires (sizeof...(is_t) == rank)
        _sp_hybrid const data_t& operator () (const is_t&... is) const 
        {
            coeff_t offst = detail::r_md_offst<0>(sizes, is...);
            return ptr[offst];
        }
    };
    
    template <typename data_t, const std::size_t rank, std::integral coeff_t>
    struct const_md_vec_image_t
    {
        using value_type = data_t;
        
        const data_t* ptr;
        ctrs::array<coeff_t, rank> sizes;
        std::size_t csize;
        
        _sp_hybrid std::size_t size() const { return csize; }
        
        _sp_hybrid const data_t* end() const
        {
            std::size_t offs = 1;
            for (int i = 0; i < rank; ++i) offs *= sizes[i];
            return ptr + offs;
        }

        template <std::integral... is_t>
        requires (sizeof...(is_t) == rank)
        _sp_hybrid const data_t& operator () (const is_t&... is) const 
        {
            coeff_t offst = detail::r_md_offst<0>(sizes, is...);
            return ptr[offst];
        }
    };
    
    template <typename vec_t>
    requires (!(requires(vec_t v) { v.transfer(); v.itransfer(); }))
    _sp_hybrid const auto make_vec_image(const vec_t& v) { return const_vec_image_t{&v[0], v.size()}; }
    
    template <typename vec_t>
    requires (!(requires(vec_t v) { v.transfer(); v.itransfer(); }))
    _sp_hybrid auto make_vec_image(vec_t& v) { return vec_image_t{&v[0], v.size()}; }
    
    template <typename vec_t, std::integral coeff_t, const std::size_t rank>
    requires (!(requires(vec_t v) { v.transfer(); v.itransfer(); }))
    _sp_hybrid const auto make_vec_image(const vec_t& v, const ctrs::array<coeff_t, rank>& sizes)
    {
        using output_t = const_md_vec_image_t<typename utils::remove_all<decltype(v[0])>::type, rank, coeff_t>;
        return output_t{&v[0], sizes, v.size()};
    }
    
    template <typename vec_t, std::integral coeff_t, const std::size_t rank>
    requires (!(requires(vec_t v) { v.transfer(); v.itransfer(); }))
    _sp_hybrid auto make_vec_image(vec_t& v, const ctrs::array<coeff_t, rank>& sizes)
    {
        using output_t = md_vec_image_t<typename utils::remove_all<decltype(v[0])>::type, rank, coeff_t>;
        return output_t{&v[0], sizes, v.size()};
    }
    
    template <typename vec_t, std::integral... is_t>
    requires (!(requires(vec_t v) { v.transfer(); v.itransfer(); }) && sizeof...(is_t) > 0)
    _sp_hybrid const auto make_vec_image(const vec_t& v, const is_t&... sizes)
    {
        using val_t    = typename utils::remove_all<decltype(v[0])>::type;
        using coeff_t  = typename std::common_type<is_t...>::type;
        using output_t = const_md_vec_image_t<val_t, sizeof...(sizes), coeff_t>;
        ctrs::array<coeff_t, sizeof...(is_t)> arr{sizes...};
        return output_t{&v[0], arr, v.size()};
    }
    
    template <typename vec_t, std::integral... is_t>
    requires (!(requires(vec_t v) { v.transfer(); v.itransfer(); }) && sizeof...(is_t) > 0)
    _sp_hybrid auto make_vec_image(vec_t& v, const is_t&... sizes)
    {
        using val_t    = typename utils::remove_all<decltype(v[0])>::type;
        using coeff_t  = typename std::common_type<is_t...>::type;
        using output_t = md_vec_image_t<val_t, sizeof...(sizes), coeff_t>;
        ctrs::array<coeff_t, sizeof...(is_t)> arr{sizes...};
        return output_t{&v[0], arr, v.size()};
    }
    
    template <typename new_t, typename old_t, const std::size_t rank, std::integral coeff_t>
    requires(sizeof(new_t) == sizeof(old_t))
    _sp_hybrid inline auto vec_img_cast(md_vec_image_t<old_t, rank, coeff_t>& old_data)
    {
        using output_t = md_vec_image_t<new_t, rank, coeff_t>;
        return output_t{(new_t*)(old_data.ptr), old_data.sizes, old_data.csize};
    }
}