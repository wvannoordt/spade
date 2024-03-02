#pragma once

#include "core/ctrs.h"

namespace spade::sampling
{
    const static struct multilinear_t
    {
        constexpr static int stencil_size() { return 8; }
        
        template <
            ctrs::basic_array                              indices_t,
            ctrs::basic_array                              coeffs_t,
            grid::multiblock_grid                          grid_t,
            ctrs::basic_array                              reduced_t,
            ctrs::basic_array                              delta_t,
            std::invocable<typename indices_t::value_type> exclude_t
            >
        requires (
            (indices_t::size() == coeffs_t::size()) &&
            (indices_t::size() >= stencil_size())
            )
        bool try_make_cloud(
            indices_t& indices,
            coeffs_t& coeffs,
            const grid_t& grid,
            const typename indices_t::value_type& landed_cell,
            const typename grid_t::coord_point_type& x_sample,
            const reduced_t& reduced_idx,
            const delta_t& deltai,
            const exclude_t& exclude_crit) const
        {
          using coeff_t = typename coeffs_t::value_type;
          using pnt_t   = typename grid_t::coord_point_type;
          int nmx = 4;
          if (grid_t::dim()==3) nmx = 8;

          grid::cell_idx_t lower_left = landed_cell;
          for (int d = 0; d < 3; ++d) lower_left.i(d) += (deltai[d]<0)?deltai[d]:0;

          indices = landed_cell;
          coeffs  = 0.0;

          bool success = true;
          for (int ct = 0; ct < nmx; ++ct)
            {
              ctrs::array<int, 3> di = 0;
              di[0] = ct%2;
              di[1] = (ct/2)%2;
              di[2] = (ct/4)%2;

                indices[ct].i()  = lower_left.i() + di[0];
                indices[ct].j()  = lower_left.j() + di[1];
                indices[ct].k()  = lower_left.k() + di[2];
                indices[ct].lb() = lower_left.lb();

                if (exclude_crit(indices[ct]))
                  {
                    success = false;
                  }


                coeff_t coeff = 1.0;
                for (int d = 0; d < grid_t::dim(); ++d)
                {
                    coeff_t dir_coeff0 = reduced_idx[d] - lower_left[d];
                    coeff_t dir_coeff1 = 1.0 - dir_coeff0;
                    coeff_t dir_coeff  = di[d]*dir_coeff0 + (1.0-di[d])*dir_coeff1;
                    coeff *= dir_coeff;
                }
                coeffs[ct] = coeff;
            }


          if (!success)
            {
              success=false;
              coeffs = 0.0;
              for (int ct = 0; ct < nmx; ++ct)
                {
                  if (!exclude_crit(indices[ct]))
                    {
                      coeffs[ct] = 1.0;
                      success = true;
                    }
                }
            }

            return success;
        }
    } multilinear;
    
    const static struct nearest_t
    {
        constexpr static int stencil_size() { return 1; }
    
        template <
            ctrs::basic_array                              indices_t,
            ctrs::basic_array                              coeffs_t,
            grid::multiblock_grid                          grid_t,
            ctrs::basic_array                              reduced_t,
            ctrs::basic_array                              delta_t,
            std::invocable<typename indices_t::value_type> exclude_t
            >
        requires (
            (indices_t::size() == coeffs_t::size()) &&
            (indices_t::size() >= stencil_size())
            )
        bool try_make_cloud(
            indices_t& indices,
            coeffs_t& coeffs,
            const grid_t& grid,
            const typename indices_t::value_type& landed_cell,
            const typename grid_t::coord_point_type& x_sample,
            const reduced_t& reduced_idx,
            const delta_t& deltai,
            const exclude_t& exclude_crit) const
        {
            indices   = landed_cell;
            coeffs    = 0.0;
            coeffs[0] = 1.0;
            
            return !exclude_crit(landed_cell);
        }
    } nearest;
    
    const static struct wlsqr_t
    {
        //Need to check this!
        constexpr static int stencil_size() { return 20; }
    
        template <
            ctrs::basic_array                              indices_t,
            ctrs::basic_array                              coeffs_t,
            grid::multiblock_grid                          grid_t,
            ctrs::basic_array                              reduced_t,
            ctrs::basic_array                              delta_t,
            std::invocable<typename indices_t::value_type> exclude_t
            >
        requires (
            (indices_t::size() == coeffs_t::size()) &&
            (indices_t::size() >= stencil_size())
            )
        bool try_make_cloud(
            indices_t& indices,
            coeffs_t& coeffs,
            const grid_t& grid,
            const typename indices_t::value_type& landed_cell,
            const typename grid_t::coord_point_type& x_sample,
            const reduced_t& reduced_idx,
            const delta_t& deltai,
            const exclude_t& exclude_crit) const
        {
          using real_t  = typename coeffs_t::value_type;
          using pnt_t   = typename grid_t::coord_point_type;
          //int nmx = 8;
          //if (grid_t::dim()==3)
          const int nmx = 12;

          std::vector<int> nx_lyr{8,};

          grid::cell_idx_t lower_left = landed_cell;
          for (int d = 0; d < 3; ++d) lower_left.i(d) += (deltai[d]<0)?deltai[d]:0;

          indices     = landed_cell;
          coeffs      = 0.0;
          int point_count      = 0;
          const int mx_lyr     = 3;
          const int order      = 1;
          const int nequations = 4;
          // create matric and rhs
          spade::linear_algebra::dense_mat<real_t, nequations> mat_a(0.);
          spade::ctrs::array<real_t, nequations> coef_mat(0.);
          std::vector<std::vector<real_t> > mat_b(nequations, std::vector<real_t>(nmx));
          spade::ctrs::array<real_t, nmx> rhs(0.);
          //
          bool enough_pts = false;
          //
          auto grid_img = grid.image(spade::device::cpu);
          //
          int num_exchange = 2;//cbrehm: needs to be adjusted
          //
          int il = -num_exchange;
          int iu =  grid.get_num_cells(0)+num_exchange;
          int jl = -num_exchange;
          int ju =  grid.get_num_cells(1)+num_exchange;
          int kl = -num_exchange;
          int ku =  grid.get_num_cells(2)+num_exchange;
          //
          for (int lyr= 0; lyr < mx_lyr; ++lyr)
            {
              const int nmx_lyr = pow(2+lyr*2,3);
              for (int ct = 0; ct < nmx_lyr; ++ct)
                {
                  ctrs::array<int, 3> di = 0;
                  const int dd = 2*(lyr+1);
                  const int dd2= dd*dd;
                  const int dd3= dd2*dd;

                  di[0] = ct     %dd-lyr;
                  di[1] =(ct/dd )%dd-lyr;
                  di[2] =(ct/dd2)%dd-lyr;
                  //
                  grid::cell_idx_t test_ipt = lower_left;
                  test_ipt.i() += di[0];
                  test_ipt.j() += di[1];
                  test_ipt.k() += di[2];
                  //
                  auto test_xpt = grid_img.get_coords(test_ipt);

                  //
                  bool in_box = true;
                  in_box = in_box && (test_ipt.i() >= il);
                  in_box = in_box && (test_ipt.j() >= jl);
                  in_box = in_box && (test_ipt.k() >= kl);
                  in_box = in_box && (test_ipt.i() <  iu);
                  in_box = in_box && (test_ipt.j() <  ju);
                  in_box = in_box && (test_ipt.k() <  ku);
                  //
                  if (!exclude_crit(test_ipt) &&
                      in_box &&
                      ((di[0]<=-lyr)||(di[0]>lyr)||
                       (di[1]<=-lyr)||(di[1]>lyr)||
                       (di[2]<=-lyr)||(di[2]>lyr)))
                    {
                      //
                      const real_t dx = test_ipt.i()-reduced_idx[0];
                      const real_t dy = test_ipt.j()-reduced_idx[1];
                      const real_t dz = test_ipt.k()-reduced_idx[2];
                      //
                      indices[point_count] = test_ipt;
                      //
                      const real_t eps = dx*1e-4;
                      real_t weight    = 1./(dx*dx + dy*dy + dz*dz + eps*eps);
                      //
                      if constexpr (order==1)
                                     {
                                       coef_mat[0] = 1.0;
                                       coef_mat[1] = dx;
                                       coef_mat[2] = dy;
                                       coef_mat[3] = dz;
                                     }
                      else if constexpr (order==2)
                                          {
                                            coef_mat[0] =  1.0;
                                            coef_mat[1] =  dx;
                                            coef_mat[2] =  dy;
                                            coef_mat[3] =  dz;
                                            coef_mat[4] =  dx*dx;
                                            coef_mat[5] =2*dy*dx;
                                            coef_mat[6] =2*dz*dx;
                                            coef_mat[7] =  dy*dy;
                                            coef_mat[8] =2*dz*dy;
                                            coef_mat[9] =  dz*dz;
                                          }
                      for (int i_row=0; i_row<nequations; ++i_row)
                        {
                          for (int i_col=0; i_col<nequations; ++i_col)
                            {
                              mat_a(i_row,i_col) += coef_mat[i_row]*coef_mat[i_col]*weight;
                            }
                          mat_b[i_row][point_count] = weight*coef_mat[i_row];
                        }
                      rhs[point_count] = weight;
                      ++point_count;
                      if (point_count==nmx)
                        {
                          enough_pts = true;
                          break;
                        }
                    }
                }
              if (enough_pts) break;
            }
          int errorflag = spade::linear_algebra::brute_invert(mat_a);
          //
          coeffs    = 0.0;
          //
          real_t sum_check =0.;
          for (int c=0;c<nmx;++c)
            {
              for (int k=0;k<nequations;++k)
                {
                  coeffs[c]+=mat_a(0,k)*mat_b[k][c];
                }
              sum_check += coeffs[c];
            }
          //
          if (!((errorflag==0)&&(std::abs(1.-sum_check)<1e-4)))
            {
              print("check     = ",errorflag);
              print("sum_check = ",(1.-sum_check));
              print("xsample   = ",x_sample);
              print("indices   = ",indices);
              print("coeffs    = ",coeffs);
              //              std::cin.get();
            }
          //
          return ((errorflag==0)&&(std::abs(1.-sum_check)<1e-4));
        }
    } wlsqr;
    
    template <typename... strategies_t>
    struct sample_strategy_cascade_t;
    
    template <typename strategy_t, typename... strategies_t>
    struct sample_strategy_cascade_t<strategy_t, strategies_t...>
    {
        using next_type = sample_strategy_cascade_t<strategies_t...>;
        constexpr static int stencil_size()
        {
            return (strategy_t::stencil_size() > next_type::stencil_size()) ? strategy_t::stencil_size() : next_type::stencil_size();
        }
        
        const strategy_t here;
        const next_type  next;
        
        template <
            ctrs::basic_array                              indices_t,
            ctrs::basic_array                              coeffs_t,
            grid::multiblock_grid                          grid_t,
            ctrs::basic_array                              reduced_t,
            ctrs::basic_array                              delta_t,
            std::invocable<typename indices_t::value_type> exclude_t
            >
        requires (
            (indices_t::size() == coeffs_t::size()) &&
            (indices_t::size() >= stencil_size())
            )
        bool try_make_cloud(
            indices_t& indices,
            coeffs_t& coeffs,
            const grid_t& grid,
            const typename indices_t::value_type& landed_cell,
            const typename grid_t::coord_point_type& x_sample,
            const reduced_t& reduced_idx,
            const delta_t& deltai,
            const exclude_t& exclude_crit) const
        {
            if (!here.try_make_cloud(indices, coeffs, grid, landed_cell, x_sample, reduced_idx, deltai, exclude_crit))
            {
                return next.try_make_cloud(indices, coeffs, grid, landed_cell, x_sample, reduced_idx, deltai, exclude_crit);
            }
            return true;
        }
    };
    
    template <typename strategy_t>
    struct sample_strategy_cascade_t<strategy_t>
    {
        constexpr static int stencil_size()
        {
            return strategy_t::stencil_size();
        }
        
        const strategy_t here;
        
        template <
            ctrs::basic_array                              indices_t,
            ctrs::basic_array                              coeffs_t,
            grid::multiblock_grid                          grid_t,
            ctrs::basic_array                              reduced_t,
            ctrs::basic_array                              delta_t,
            std::invocable<typename indices_t::value_type> exclude_t
            >
        requires (
            (indices_t::size() == coeffs_t::size()) &&
            (indices_t::size() >= stencil_size())
            )
        bool try_make_cloud(
            indices_t& indices,
            coeffs_t& coeffs,
            const grid_t& grid,
            const typename indices_t::value_type& landed_cell,
            const typename grid_t::coord_point_type& x_sample,
            const reduced_t& reduced_idx,
            const delta_t& deltai,
            const exclude_t& exclude_crit) const
        {
            return here.try_make_cloud(indices, coeffs, grid, landed_cell, x_sample, reduced_idx, deltai, exclude_crit);
        }
    };
    
    template <typename... strategies_t>
    inline auto prioritize(const strategies_t&... strategies)
    {
        return sample_strategy_cascade_t<strategies_t...>{strategies...};
    }
}
