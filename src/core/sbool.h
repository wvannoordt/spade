#pragma once

namespace spade::utils
{
    // This whole type is here to get around the stupid std::vector<bool> issue:
    // https://en.cppreference.com/w/cpp/container/vector_bool
    struct sbool
    {
        bool val;
        _sp_hybrid operator bool() const       { return val; }
        _sp_hybrid sbool() = default;
        _sp_hybrid sbool(const bool d) {val = d;}
        _sp_hybrid sbool& operator = (const bool rhs) { val = rhs; return *this; }
    };
}