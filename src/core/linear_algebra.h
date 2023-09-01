#pragma once

#include <initializer_list>

#include "core/cuda_incl.h"
#include "core/ctrs.h"
#include "core/print.h"

namespace spade::linear_algebra
{
    template <class T, const size_t size, typename real_type>
    concept mat_nd = requires(T t, size_t i, size_t j, typename T::vector_type rhs)
    {
        { t(i,j) } -> std::common_with<real_type>;
        typename T::value_type;
        // { t*rhs } -> ctrs::basic_array;
    };

    template <typename data_t, const std::size_t mat_size> struct dense_mat
    {
        typedef ctrs::array<data_t, mat_size> vector_type;
        typedef data_t value_type;
        constexpr static std::size_t size() {return mat_size;}
        data_t m[mat_size*mat_size];
        _sp_hybrid dense_mat(){}
        _sp_hybrid dense_mat(const value_type& init) {for (int iii = 0; iii < mat_size*mat_size; ++iii) m[iii] = init;}
        _sp_hybrid dense_mat(const std::initializer_list<std::initializer_list<data_t>>& init_list)
        {
            std::size_t i = 0;
            for (const auto& outer: init_list)
            {
                for (const auto& inner: outer)
                {
                    m[i] = inner;
                    ++i;
                }
            }
        }
        _sp_hybrid data_t& operator() (const std::size_t& row, const std::size_t& col) {return m[col+row*mat_size];}
        _sp_hybrid const data_t& operator() (const std::size_t& row, const std::size_t& col) const {return m[col+row*mat_size];}
        
        template <ctrs::basic_array vec_t>
        _sp_hybrid dense_mat<data_t, mat_size>& set_row(const std::size_t i, const vec_t& v)
        {
            for (std::size_t j = 0; j < mat_size; j++) (*this)(i,j)=v[j];
            return *this;
        }
        
        template <ctrs::basic_array vec_t>
        _sp_hybrid dense_mat<data_t, mat_size>& set_col(const std::size_t j, const vec_t& v)
        {
            for (std::size_t i = 0; i < mat_size; i++) (*this)(i,j)=v[i];
            return *this;
        }
        
        _sp_hybrid vector_type operator * (const vector_type& rhs) const
        {
            vector_type output;
            for (std::size_t i = 0; i < mat_size; ++i)
            {
                value_type sum = 0.0;
                for (std::size_t j = 0; j < mat_size; ++j)
                {
                    sum += (*this)(i,j)*rhs[j];
                }
                output[i] = sum;
            }
            return output;
        }
        
        _sp_hybrid dense_mat<data_t, mat_size> operator - (const dense_mat<data_t, mat_size>& rhs) const
        {
            dense_mat<data_t, mat_size> output;
            for (std::size_t i = 0; i < mat_size*mat_size; ++i) output.m[i] = m[i]-rhs.m[i];
            return output;
        }
        
        _sp_hybrid dense_mat<data_t, mat_size> operator + (const dense_mat<data_t, mat_size>& rhs) const
        {
            dense_mat<data_t, mat_size> output;
            for (std::size_t i = 0; i < mat_size*mat_size; ++i) output.m[i] = m[i]+rhs.m[i];
            return output;
        }
        
        _sp_hybrid  data_t det() const
        {
            static_assert(mat_size==3, "Determinant calculation not yet implemented for matrices other thatn 3x3");
            return (*this)(0,0)*((*this)(1,1)*(*this)(2,2) - (*this)(1,2)*(*this)(2,1))
                +  (*this)(0,1)*((*this)(1,2)*(*this)(2,0) - (*this)(1,0)*(*this)(2,2))
                +  (*this)(0,2)*((*this)(1,0)*(*this)(2,1) - (*this)(1,1)*(*this)(2,0));
        }

        template <typename kern_t>
        _sp_hybrid dense_mat& fill(const kern_t& kern)
        {
            for (int j = 0; j < mat_size; ++j)
            {
                for (int i = 0; i < mat_size; ++i)
                {
                    (*this)(i,j) = kern(i,j);
                }
            }
            return *this;
        }
    };
    
    template <typename rhs_t, typename data_t, const std::size_t mat_size>
    _sp_hybrid  dense_mat<data_t, mat_size> operator * (const rhs_t& rhs, const dense_mat<data_t, mat_size>& mat)
    {
        dense_mat<data_t, mat_size> output;
        for (std::size_t i = 0; i < mat_size*mat_size; ++i) output.m[i] = mat.m[i]*rhs;
        return output;
    }
    
    template <typename data_t, const std::size_t sys_size, ctrs::basic_array rhs_t>
    void brute_solve_inplace(dense_mat<data_t, sys_size>& a, rhs_t& y)
    {
        for (int k = 0; k < sys_size; ++k)
        {
            auto diag = a(k,k);
            for (int r = 0; r < sys_size; ++r) a(k,r)/=diag;
            y[k]/=diag;
            for (int r = k+1; r < sys_size; ++r)
            {
                auto root_elem = a(r, k);
                for (int c = k; c < sys_size; ++c)
                {
                    a(r,c) -= root_elem*a(k,c);
                }
                y[k] *= root_elem;
                y[r] -= y[k];
                y[k] /= root_elem;
            }
        }
        for (int k = sys_size-1; k>=0; --k)
        {
            for (int r = k-1; r >= 0; --r)
            {
                auto root_elem = a(r,k);
                for (int c = sys_size-1; c >= k; --c)
                {
                    a(r,c) -= root_elem*a(k,c);
                }
                y[k] *= root_elem;
                y[r]-=y[k];
                y[k] /= root_elem;
            }
        }
    }
}
