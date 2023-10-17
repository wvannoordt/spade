#pragma once

#include <vector>
#include "core/cuda_incl.h"

namespace spade::utils
{
    template <typename data_t>
    struct vec_image_t
    {
        using value_type = data_t;
        data_t* ptr;
        std::size_t csize;
        _sp_hybrid std::size_t size() const { return csize; }
        _sp_hybrid data_t& operator [] (const std::size_t& idx) { return ptr[idx]; }
        _sp_hybrid const data_t& operator [] (const std::size_t& idx) const { return ptr[idx]; }
    };
    
    template <typename data_t>
    struct const_vec_image_t
    {
        const data_t* ptr;
        std::size_t csize;
        _sp_hybrid std::size_t size() const { return csize; }
        _sp_hybrid const data_t& operator [] (const std::size_t& idx) const { return ptr[idx]; }
    };
    
    template <typename vec_t>
    const auto make_vec_image(const vec_t& v) { return const_vec_image_t{&v[0], v.size()}; }
    
    template <typename vec_t>
    auto make_vec_image(vec_t& v) { return vec_image_t{&v[0], v.size()}; }
}