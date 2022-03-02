#pragma once

namespace cvdf::output
{
    template <class T> concept exportable_type = requires(T t)
    {
        //todo: write this
        t;
    };
}