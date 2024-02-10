#pragma once

namespace spade::utils
{
    // This whole type is here to get around the stupid std::vector<bool> issue:
    // https://en.cppreference.com/w/cpp/container/vector_bool
    struct sbool
    {
        bool val;
        operator bool() const       { return val; }
        sbool() = default;
        sbool(const bool d) {val = d;}
        sbool& operator = (const bool rhs) { val = rhs; return *this; }
    };
}