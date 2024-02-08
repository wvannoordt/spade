#pragma once

namespace spade::utils
{
    // This whole type is here to get aroun the stupid std::vector<bool> issue:
    // https://en.cppreference.com/w/cpp/container/vector_bool
    struct sbool
    {
        bool val;
        operator bool() const       { return val; }
        sbool& operator = (const bool rhs) { val = rhs; return *this; }
    };
}