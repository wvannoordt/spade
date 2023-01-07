#pragma once

#include <vector>
#include <algorithm>

#include "array-containers/ac_base.h"

namespace spade::array_containers
{
    template <typename data_t> struct vector_wrapper : public array_container<data_t, vector_wrapper<data_t>>
    {
        using base_t = array_container<data_t, vector_wrapper<data_t>>;
        using size_type = base_t::size_type;
        using value_type = base_t::value_type;
        std::vector<data_t> data;
        virtual void resize(const size_type& new_size, const value_type val = data_t()) {data.resize(new_size, val);};
        virtual size_type size() const {return data.size();};
        virtual void fill(const data_t& val) {std::fill(data.begin(), data.end(), val);};
        virtual value_type& operator [] (const size_type i) {return data[i];};
        virtual const value_type& operator [] (const size_type i) const {return data[i];};
    };
}