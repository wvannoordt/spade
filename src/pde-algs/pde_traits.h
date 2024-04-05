#pragma once

#include "sym/sym.h"

namespace spade::pde_algs
{
    struct pde_increment_base_t
    {
        constexpr static auto trait_label() { using namespace sym::literals; return "pde_increment"_sym; }
    };
    
    static struct toverwrite_t : public pde_increment_base_t
    {
        constexpr static bool increment_mode = false;
    } overwrite;
    
    static struct tincrement_t : public pde_increment_base_t
    {
        constexpr static bool increment_mode = true;
    } increment;
}