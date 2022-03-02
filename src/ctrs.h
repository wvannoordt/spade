#pragma once

#include <type_traits>
#include <concepts>

namespace cvdf::ctrs
{
    template <class T, const size_t size, typename real_type> concept vec_nd = (sizeof(T) == size*sizeof(real_type)) && requires(T t, size_t i)
    {
        { t[i] } -> std::common_with<real_type>;
    };
    
    template<typename dtype, const size_t ar_size> struct array
    {
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
    };
}