#pragma once

#include <concepts>

//TODO

namespace cvdf::block_config
{
    template <typename T> concept block_config = requires(T t)
    {
        t;
    };
}