#pragma once
#include <type_traits>
#include <concepts>
namespace cvdf::grid
{
    template <typename derived_t, typename arith_t>
    struct arithmetic_wrapper_t
    {
        arith_t val;
        arithmetic_wrapper_t(void) {}
        arithmetic_wrapper_t(const arithmetic_wrapper_t& rhs) {val = rhs.val;}
        arithmetic_wrapper_t(const arith_t& rhs) {val = rhs;}
        
        derived_t& operator = (const arith_t& rhs)
        {
            val = rhs;
            return static_cast<derived_t&>(*this);
        }
        
        derived_t& operator /= (const arith_t& rhs)
        {
            val /= rhs;
            return static_cast<derived_t&>(*this);
        }
        derived_t& operator *= (const arith_t& rhs)
        {
            val *= rhs;
            return static_cast<derived_t&>(*this);
        }
        
        derived_t& operator += (const arith_t& rhs)
        {
            val += rhs;
            return static_cast<derived_t&>(*this);
        }
        derived_t& operator -= (const arith_t& rhs)
        {
            val -= rhs;
            return static_cast<derived_t&>(*this);
        }
        operator arith_t() const { return val; }
    };
    
    template <typename arith_t> struct cell_t : public arithmetic_wrapper_t<cell_t<arith_t>, arith_t>
    {
        typedef arithmetic_wrapper_t<cell_t<arith_t>, arith_t> base_t;
        typedef arith_t value_type;
        using base_t::base_t;
    };
    template <typename arith_t> struct node_t : public arithmetic_wrapper_t<node_t<arith_t>, arith_t>
    {
        typedef arithmetic_wrapper_t<node_t<arith_t>, arith_t> base_t;
        typedef arith_t value_type;
        using base_t::base_t;
    };
    template <typename arith_t> struct face_t : public arithmetic_wrapper_t<face_t<arith_t>, arith_t>
    {
        typedef arithmetic_wrapper_t<face_t<arith_t>, arith_t> base_t;
        typedef arith_t value_type;
        using base_t::base_t;
    };
}