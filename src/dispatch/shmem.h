#pragma once

#include "core/ctrs.h"
#include "core/vec_image.h"
#include "core/const_index.h"

namespace spade::dispatch::shmem
{
    struct empty_t
    {
        _sp_hybrid constexpr static bool is_empty() { return true; }
        _sp_hybrid constexpr static std::size_t bytes() { return 0; }
    };
    
    template <typename data_t>
    using vec_t = utils::vec_image_t<data_t>;
    
    template <typename data_t>
    vec_t<data_t> vec(const std::size_t& size)
    {
        return vec_t<data_t>{nullptr, size};
    }
    
    template <typename... datas_t> struct shmem_t;
    
    template<typename data_t, typename... datas_t>
    struct shmem_t<data_t, datas_t...>
    {
        _sp_hybrid constexpr static bool is_empty() { return false; }
        _sp_hybrid std::size_t bytes() const { return sizeof(data_t) + next.bytes(); }
        volatile data_t* data = nullptr;
        shmem_t<datas_t...> next;
        
        shmem_t(const data_t&, const datas_t&... datas) : next(datas...) {}
        _sp_hybrid void bind_ptr(volatile char* ptr)
        {
            data = (data_t*)ptr;
            next.bind_ptr(ptr + sizeof(data_t));
        }
        
        template <const udci::integral_t idx>
        _sp_hybrid auto& operator [] (const udci::idx_const_t<idx>&)
        {
            if constexpr (idx == 0) return *(data);
            else return next[udci::idx_const_t<idx-1>()];
        }
        
        template <const int idx> using type_at = typename std::conditional<idx==0, data_t, typename shmem_t<datas_t...>::type_at<idx-1>>::type;
    };
    
    template <typename data_t, typename... datas_t>
    struct shmem_t<vec_t<data_t>, datas_t...>
    {
        _sp_hybrid constexpr static bool is_empty() { return false; }
        _sp_hybrid std::size_t bytes() const { return data.csize*sizeof(data_t) + next.bytes(); }
        vec_t<data_t> data;
        shmem_t<datas_t...> next;
        
        shmem_t(const vec_t<data_t>& data_in, const datas_t&... datas) : data{nullptr, data_in.csize}, next(datas...) {}
        _sp_hybrid void bind_ptr(volatile char* ptr)
        {
            data.ptr = (data_t*)ptr;
            next.bind_ptr(ptr + data.csize*sizeof(data_t));
        }
        
        template <const udci::integral_t idx>
        _sp_hybrid auto& operator [] (const udci::idx_const_t<idx>&)
        {
            if constexpr (idx == 0) return data;
            else return next[udci::idx_const_t<idx-1>()];
        }
        
        template <const int idx> using type_at = typename std::conditional<idx==0, vec_t<data_t>, typename shmem_t<datas_t...>::type_at<idx-1>>::type;
    };
    
    template <typename data_t>
    struct shmem_t<vec_t<data_t>>
    {
        _sp_hybrid constexpr static bool is_empty() { return false; }
        _sp_hybrid std::size_t bytes() const { return data.csize*sizeof(data_t); }
        vec_t<data_t> data;
        
        shmem_t(const vec_t<data_t>& data_in) : data{nullptr, data_in.csize} {}
        _sp_hybrid void bind_ptr(volatile char* ptr)
        {
            data.ptr = (data_t*)ptr;
        }
        
        _sp_hybrid auto& operator [] (const udci::idx_const_t<0>&)
        {
            return data;
        }
        
        template <const int idx> using type_at = vec_t<data_t>;
    };
    
    template <typename data_t>
    struct shmem_t<data_t>
    {
        _sp_hybrid constexpr static bool is_empty() { return false; }
        _sp_hybrid std::size_t bytes() const { return sizeof(data_t); }
        volatile data_t* data = nullptr;
        
        shmem_t(const data_t&) {}
        _sp_hybrid void bind_ptr(volatile char* ptr)
        {
            data = (data_t*)ptr;
        }
        
        _sp_hybrid auto& operator [] (const udci::idx_const_t<0>&)
        {
            return *data;
        }
        
        template <const int idx> using type_at = data_t;
    };
    
    
    template <typename... datas_t>
    auto make_shmem(const datas_t&... datas)
    {
        return shmem_t<datas_t...>(datas...);
    }
}