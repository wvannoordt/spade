#pragma once

#include <concepts>

#include "core/cuda_incl.h"
#include "sym/symbol.h"

namespace spade::ode
{
    template <sym::is_symbol sym_t, typename func_t, typename specifier_t>
    struct expression_t
    {
        using symbol_type = sym_t;
        const func_t func;
        const specifier_t specifier;
        expression_t(const sym_t&, const func_t& f_in, const specifier_t& spec)
        : func{f_in}, specifier{spec} {}
    };
    
    template <typename val_t>
    struct tridiag_values_t
    {
        val_t lower, diag, upper, rhs;
    };
    
    template <typename val_t>
    _sp_hybrid static tridiag_values_t<val_t> tridiag(const val_t& l, const val_t& d, const val_t& u, const val_t& r)
    {
        return tridiag_values_t{l,d,u,r};
    }
    
    template <typename... exprs_t> struct expression_list_t;
    template <typename expr_t, typename... exprs_t>
    struct expression_list_t<expr_t, exprs_t...>
    {
        using top_sym   = typename expr_t::symbol_type;
        using next_type = expression_list_t<exprs_t...>;
        using syms_type = typename next_type::syms_type::extend<typename expr_t::symbol_type>;
        expr_t expr;
        next_type next;
    };
    
    template <> struct expression_list_t<>
    {
        using top_sym   = int;
        using syms_type = sym::vector_t<>;
    };
    
    static constexpr auto make_expressions()
    {
        return expression_list_t<>{};
    }
    
    template <typename expr_t, typename... exprs_t>
    static constexpr auto make_expressions(const expr_t& expr, const exprs_t&... exprs)
    {
        return expression_list_t<expr_t, exprs_t...>{expr, make_expressions(exprs...)};
    }
    
    template <sym::is_symbol sym_t, typename... exprs_t>
    // requires (something something)
    _sp_hybrid static constexpr auto get_expression(const sym_t& sym, const expression_list_t<exprs_t...>& expressions)
    {
        if constexpr(std::same_as<sym_t, typename expression_list_t<exprs_t...>::top_sym>)
        {
            return expressions.expr;
        }
        else
        {
            return get_expression(sym, expressions.next);
        }
    }
}