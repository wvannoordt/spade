#pragma once

#include <concepts>

//TODO

namespace spade::block_config
{
    template <typename T> concept block_config = requires(T t)
    {
        t;
    };
}