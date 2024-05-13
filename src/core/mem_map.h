#pragma once

#include <type_traits>

#include "core/config.h"
#include "core/aliases.h" 
#include "core/static_for.h"
#include "core/cuda_incl.h"
#include "core/permute.h"

namespace spade::mem_map
{
    template <const int i0, const int i1> struct static_dim_t
    {
        static constexpr std::size_t size() {return i1-i0;}
        static constexpr int start()        {return i0;}
        static constexpr int end()          {return i1;}
        static constexpr int rank()         {return 1;}
        static constexpr int num_coeffs()   {return 1;}
    };

    template <typename integral_t=int> struct dynamic_dim_t
    {
        integral_t i0, i1;
        std::size_t size()  const           {return i1-i0;}
        integral_t  start() const           {return i0;}
        integral_t  end()   const           {return i1;}
        static constexpr int rank()        {return 1;}
        static constexpr int num_coeffs()  {return 1;}
        dynamic_dim_t(){}
        dynamic_dim_t(integral_t i0_in, integral_t i1_in):i0{i0_in},i1{i1_in} {}
    };

    namespace detail
    {
        template <typename dim_t, typename... dims_t> struct rank_sum_t
        {
            constexpr static int value = dim_t::rank() + rank_sum_t<dims_t...>::value;
        };
        
        template <typename dim_t> struct rank_sum_t<dim_t> {constexpr static int value = dim_t::rank();};
        
        template <typename dim_t, typename... dims_t> struct coeff_count_sum_t
        {
            constexpr static int value = dim_t::num_coeffs() + coeff_count_sum_t<dims_t...>::value;
        };
        
        template <typename dim_t> struct coeff_count_sum_t<dim_t> {constexpr static int value = dim_t::num_coeffs();};
    }

    using offset_type = long int;

    struct empty_view_t
    {
        std::tuple<> views;
        _sp_hybrid constexpr static int rank()        {return 0;}
        _sp_hybrid constexpr static int num_views()   {return 0;}
        _sp_hybrid constexpr static int num_coeffs()  {return 0;}
        using coeff_array_t = ctrs::array<int, 0>;
        coeff_array_t get_extent_array() const
        {
            return coeff_array_t();
        }
        coeff_array_t get_lower_bound_array() const
        {
            return coeff_array_t();
        }
    };

    template <typename... dims_t> struct recti_view_t
    {
        _sp_hybrid constexpr static int rank()        {return detail::rank_sum_t<dims_t...>::value;}
        _sp_hybrid constexpr static int num_views()   {return sizeof...(dims_t);}
        _sp_hybrid constexpr static int num_coeffs()  {return detail::coeff_count_sum_t<dims_t...>::value;}
        std::tuple<dims_t...> views;
        recti_view_t(){}
        recti_view_t(dims_t... views_in) : views(std::make_tuple(views_in...)) {}
        
        using coeff_array_t = ctrs::array<int, num_coeffs()>;
        
        coeff_array_t get_extent_array() const
        {
            coeff_array_t output;
            int idx = 0;
            algs::static_for<0,num_views()>([&](const auto& ii) -> void
            {
                const int i = ii.value;
                const auto& view = std::get<i>(views);
                if constexpr (requires{view.size();})
                {
                    output[idx++] = view.size();
                }
                else
                {
                    auto sub_coeffs = view.get_extent_array();
                    for (const auto& c: sub_coeffs) output[idx++] = c;
                }
            });
            return output;
        }
        
        //TODO: remove duplication
        coeff_array_t get_lower_bound_array() const
        {
            coeff_array_t output;
            int idx = 0;
            algs::static_for<0,num_views()>([&](const auto& ii) -> void
            {
                const int i = ii.value;
                const auto& view = std::get<i>(views);
                if constexpr (requires{view.size();})
                {
                    output[idx++] = -view.start();
                }
                else
                {
                    auto sub_coeffs = view.get_lower_bound_array();
                    for (const auto& c: sub_coeffs) output[idx++] = c;
                }
            });
            return output;
        }
    };
    
    template <const int tsize, typename dim_t> struct tile_uni_dim_t
    {
        using dim_type = dim_t;
        dim_t dim;
        constexpr static int tile_size() {return static_math::pow<2,tsize>::value;}
        tile_uni_dim_t(const dim_t& dim_in) : dim{dim_in} {}
        std::size_t size()  const          {return dim.size();}
        _sp_hybrid constexpr static int rank()        {return 1;}
        _sp_hybrid constexpr static int num_coeffs()  {return 2;}
        _sp_hybrid constexpr static int mask()        {return tile_size()-1;}
        _sp_hybrid constexpr static int exp()         {return tsize;}
    };
    
    template <typename... uni_dims_t>
    requires (((uni_dims_t::dim_type::rank()==1) && (uni_dims_t::dim_type::num_coeffs()==1)) && ...)
    struct tiled_dim_t
    {
        constexpr static bool tile_multi_dim_tag = true;
        template <const int idx> constexpr static int tile_size = utils::get_pack_type<idx, uni_dims_t...>::type::tile_size();
        std::tuple<uni_dims_t...> views;
        tiled_dim_t(const uni_dims_t&... views_in) : views(std::make_tuple(views_in...)) {}
        _sp_hybrid constexpr static int num_views()   {return sizeof...(uni_dims_t);}
        _sp_hybrid constexpr static int rank()        {return sizeof...(uni_dims_t);}
        _sp_hybrid constexpr static int num_coeffs()  {return 2*rank();}
        using coeff_array_t = ctrs::array<int, num_coeffs()>;
        
        coeff_array_t get_extent_array() const { return get_extent_array([](const auto& vw) -> auto {return vw.size();});}
        template <typename accessor_t> coeff_array_t get_extent_array(const accessor_t& accessor) const
        {
            coeff_array_t output;
            int idx = 0;
            algs::static_for<0,rank()>([&](const auto& ii) -> void
            {
                const int i = ii.value;
                const auto& view = std::get<i>(views);

                output[i]        = view.tile_size();
                output[i+rank()] = view.dim.size()/view.tile_size();
            });
            return output;
        }
        coeff_array_t get_lower_bound_array() const
        {
            coeff_array_t output;
            int idx = 0;
            algs::static_for<0,rank()>([&](const auto& ii) -> void
            {
                const int i = ii.value;
                const auto& view = std::get<i>(views);

                output[i]        = 0;
                output[i+rank()] = 0;
            });
            return output;
        }
    };

    namespace detail
    {
        template <const int idx, const int view_idx, const int cur_pos, typename views_t>
        _sp_hybrid constexpr const auto& recurse_dim_retrieve(const views_t& view_collection)
        {
            const auto& view = std::get<view_idx>(view_collection.views);
            
            if constexpr(view.rank()==0)
            {
                return recurse_dim_retrieve<idx, view_idx+1, cur_pos>(view_collection);
            }
            else
            {
                //  0   1   2   3   4   5   6   7
                //  +   +   +   + [ +   +   + ] +
                //                  0   1   2
                //  0   1   2   3   4           5
                
                //Here we have hit a match
                if constexpr (idx == cur_pos)
                {
                    //if the current dimension is singular, just return it
                    //otherwise, it's a multi-dimension and we simply return
                    //the first element
                    if constexpr (requires{view.size();}) { return view; }
                    else { return std::get<0>(view.views); }
                }
                else
                {
                    //Simply increment to the next view
                    if constexpr (requires{view.size();})
                    {
                        return recurse_dim_retrieve<idx, view_idx+1, cur_pos+1>(view_collection);
                    }
                    else
                    {
                        //check if the index in question is inside the current multiview
                        const int num_views = view.num_views();
                        if constexpr (idx >= cur_pos+num_views)
                        {
                            //if not, simply continue pase the end of the multiview
                            return recurse_dim_retrieve<idx, view_idx+1, cur_pos+num_views>(view_collection);
                        }
                        else
                        {
                            //if so, then we index into the multiview but relative to its
                            //starting position (which is the current position)
                            return recurse_dim_retrieve<idx-cur_pos, 0, 0>(view);
                        }
                    }
                }
            }
        }
    }
    
    template <const int idx, typename views_t>
    _sp_hybrid constexpr const auto& get_dim_from_map(const views_t& views)
    {
        return detail::recurse_dim_retrieve<idx, 0, 0>(views);
    }
    
    namespace detail
    {
        template <typename idx_t>          struct ar_idx_size        {constexpr static int value = 1;};
        template <ctrs::basic_array idx_t> struct ar_idx_size<idx_t> {constexpr static int value = idx_t::size();};
        
        template <typename idx_t, typename... idxs_t> struct size_sum
        {
            constexpr static int value = ar_idx_size<idx_t>::value + size_sum<idxs_t...>::value;
        };
        template <typename idx_t> struct size_sum<idx_t> {constexpr static int value = ar_idx_size<idx_t>::value;};
    }
    
    template <typename map_t> struct mem_map_t
    {
        map_t mmap;
        ctrs::array<std::size_t, map_t::num_coeffs()> i_coeff;
        int offset_base;
        mem_map_t(){}
        mem_map_t(const map_t& map_in):mmap{map_in}
        {
            compute_coeffs();
            compute_offset_base();
        }
        
        void compute_coeffs()
        {
            auto sizes = mmap.get_extent_array();
            for (int i = 0; i < i_coeff.size(); ++i)
            {
                i_coeff[i] = 1;
                for (int j = 0; j < i; ++j) i_coeff[i] *= sizes[j];
            }
            permute_tile_coeffs();
        }
        
        void compute_coeffs(const ctrs::array<int, map_t::num_coeffs()> permutation)
        {
            auto sizes = mmap.get_extent_array();
            sizes = utils::permute(sizes, permutation);
            
            for (int i = 0; i < i_coeff.size(); ++i)
            {
                i_coeff[i] = 1;
                for (int j = 0; j < i; ++j) i_coeff[i] *= sizes[j];
            }
            i_coeff = utils::ipermute(i_coeff, permutation);
            permute_tile_coeffs();
        }
        
        void permute_tile_coeffs()
        {
            int base_index = 0;
            algs::static_for<0, map_t::num_views()>([&](const auto& ii) -> void
            {
                const int i = ii.value;
                const auto& dim = std::get<i>(mmap.views);
                if constexpr (requires{dim.tile_multi_dim_tag;})
                {
                    const int nc = dim.rank();
                    ctrs::array<std::size_t, 2*nc> to_insert = 0;
                    algs::static_for<0,nc>([&](const auto& jj) -> void
                    {
                        const int j = jj.value;
                        to_insert[2*j]   = i_coeff[base_index+j];
                        to_insert[2*j+1] = i_coeff[base_index+j+nc];
                    });
                    for (int k = 0; k < to_insert.size(); ++k)
                    {
                        i_coeff[base_index+k] = to_insert[k];
                    }
                    base_index+=2*nc;
                }
                else
                {
                    base_index += dim.num_coeffs();
                }
            });
        }
        
        void compute_offset_base()
        {
            offset_base = 0;
            auto sizes = mmap.get_lower_bound_array();
            int idx = 0;
            for (const auto& ii: sizes) offset_base += i_coeff[idx++]*ii;
        }
        
        template <
            const int ar_idx,    //for the current provided index, what element are we currently concerned with?
            const int coeff_idx, //What is the index of the coefficient to be multiplied?
            const int dim_idx,   //Which dimension are we considering?
            typename offset_t,
            typename idx_t, typename... idxs_t>
        _sp_hybrid void rec_off_calc(offset_t& offset, const idx_t& idx, const idxs_t&... idxs) const
        {
            constexpr bool skip_index = [&]()
            {
                if constexpr (ctrs::basic_array<idx_t>) { return (idx_t::size()==0); }
                else { return false; }
            }();
            if constexpr (skip_index)
            {
                rec_off_calc<0, coeff_idx, dim_idx>(offset, idxs...);
            }
            else
            {
                const auto idx_val = [&]()->auto
                {
                    if constexpr (ctrs::basic_array<idx_t>) return idx[ar_idx];
                    else return idx;
                }();
                
                const auto& dim = get_dim_from_map<dim_idx>(mmap);
                
                //need to detect the current dimension somehow.
                const int next_coeff = dim.num_coeffs();
                
                if constexpr(next_coeff == 1)
                {
                    //increment the offset according to a conventional index
                    offset += idx_val*i_coeff[coeff_idx];
                }
                else
                {
                    const int idx_val_maj = (idx_val-dim.dim.start()) >> dim.exp();
                    const int idx_val_min = (idx_val-dim.dim.start()) & dim.mask();
                    //tiling logic
                    offset += idx_val_min*i_coeff[coeff_idx];
                    offset += idx_val_maj*i_coeff[coeff_idx+1];
                }
                
                //logic for incrementing the
                //coefficient index...
                if constexpr (ctrs::basic_array<idx_t>)
                {
                    if constexpr(ar_idx == idx_t::size()-1)
                    {
                        if constexpr (sizeof...(idxs) > 0)
                        {
                            rec_off_calc<0, coeff_idx+next_coeff, dim_idx+1>(offset, idxs...);
                        }
                    }
                    else
                    {
                        rec_off_calc<ar_idx+1, coeff_idx+next_coeff, dim_idx+1>(offset, idx, idxs...);
                    }
                }
                else
                {
                    if constexpr (sizeof...(idxs) > 0)
                    {
                        rec_off_calc<0,coeff_idx+next_coeff, dim_idx+1>(offset, idxs...);
                    }
                }
            }
        }
        
        template <typename... idxs_t>
        requires (detail::size_sum<idxs_t...>::value == map_t::rank())
        _sp_hybrid offset_type compute_offset(const idxs_t&... idxs) const
        {
            offset_type output = 0;
            rec_off_calc<0,0,0>(output, idxs...);
            return output + offset_base;
        }
    };

    template <typename view_t> static constexpr std::size_t map_size(const view_t& view)
    {
        std::size_t out = 1;
        const auto sizes = view.get_extent_array();
        for (const auto p:sizes) out *= p;
        return out;
    }
    
    template <typename T> concept integral_array           = ctrs::basic_array<T> && std::integral<typename T::value_type>;
    template <typename T> concept integral_scalar_or_array = std::integral<T> || integral_array<T>;
    
    template <const int idx, typename i_t, integral_scalar_or_array... is_t>
    _sp_hybrid _sp_inline static constexpr int get_val(const i_t& i)
    {
        if constexpr (ctrs::basic_array<i_t>)
        {
            return i[idx];
        }
        else
        {
            return i;
        }
    }
    
    template <const int idx, integral_scalar_or_array i_t, integral_scalar_or_array... is_t>
    requires(sizeof...(is_t) > 0)
    _sp_hybrid _sp_inline static constexpr int get_val(const i_t& i, const is_t&... is)
    {
        if constexpr (ctrs::basic_array<i_t>)
        {
            if constexpr (idx < i_t::size())
            {
                return i[idx];
            }
            else
            {
                return get_val<idx - i_t::size()>(is...);
            }
        }
        else
        {
            if constexpr (idx == 0)
            {
                return i;
            }
            else
            {
                return get_val<idx - 1>(is...);
            }
        }
    }
    
    // Re-implementing memory maps with a slightly less cumbersome syntax, etc
    // Memory map tags
    
    inline struct tlinear_t {} linear;
    template <const int rank>
    struct linear_t
    {
        using tag_type = tlinear_t;
        ctrs::array<int, rank> min, max;
        _sp_hybrid constexpr int size(int i) const { return max[i] - min[i]; }
        _sp_hybrid constexpr std::size_t volume() const
        {
            std::size_t output = 1;
            algs::static_for<0, rank>([&](const auto& ii)
            {
                constexpr int i = ii.value;
                output *= size(i);
            });
            return output;
        }
        
        template <integral_scalar_or_array... is_t>
        _sp_hybrid _sp_inline constexpr std::size_t compute_offset(const is_t&... is) const noexcept
        {
            // off = v + i*nv + j*nv*ni + k*nv*ni*nj + l*nv*ni*nj*nk;
            std::size_t output = 0;
            algs::static_for<0, rank-1>([&](const auto& ii)
            {
                constexpr int i = ii.value;
                output += get_val<rank - i - 1>(is...) - min[rank - i - 1];
                output *= size(rank - i - 2);
            });
            output += get_val<0>(is...) - min[0];
            return output;
        }
    };
    
    inline struct ttiled_t {} tiled;
    //Todo: make this more generic
    struct tiled_t
    {
        using tag_type = ttiled_t;
        ctrs::array<int, 5> min, max;
        _sp_hybrid int size(int i) const { return max[i] - min[i]; }
        _sp_hybrid std::size_t volume() const
        {
            std::size_t output = 1;
            algs::static_for<0, 5>([&](const auto& ii)
            {
                constexpr int i = ii.value;
                output *= size(i);
            });
            return output;
        }
        template <integral_scalar_or_array... is_t>
        _sp_hybrid _sp_inline std::size_t compute_offset(const is_t&... is) const noexcept
        {
            constexpr int tpow = 2;
            constexpr int tsz  = 1 << tpow;
            constexpr int mask = tsz - 1;
            
            int ni  = size(1);
            int nj  = size(2);
            int nk  = size(3);
            int nlb = size(4);
            int nv  = size(0);
            
            ni = ni >> tpow;
            nj = nj >> tpow;
            nk = nk >> tpow;
            
            int v  = get_val<0>(is...) - min[0];
            int i  = get_val<1>(is...) - min[1];
            int j  = get_val<2>(is...) - min[2];
            int k  = get_val<3>(is...) - min[3];
            int lb = get_val<4>(is...) - min[4];
            
            int ii = i & mask;
            i = i >> tpow;
            
            int jj = j & mask;
            j = j >> tpow;
            
            int kk = k & mask;
            k = k >> tpow;
            
            // std::size_t output = v;
            // output *= nlb;
            // output += lb;
            // output *= nk;
            // output += k;
            // output *= nj;
            // output += j;
            // output *= ni;
            // output += i;
            // output *= tsz;
            // output += kk;
            // output *= tsz;
            // output += jj;
            // output *= tsz;
            // output += ii;
            
            std::size_t output = lb;
            output *= nk;
            output += k;
            output *= nj;
            output += j;
            output *= ni;
            output += i;
            output *= nv;
            output += v;
            output *= tsz;
            output += kk;
            output *= tsz;
            output += jj;
            output *= tsz;
            output += ii;
            
            return output;
        }
    };
    
    inline struct ttiled_small_t {} tiled_small;
    //Todo: make this more generic
    struct tiled_small_t
    {
        using tag_type = ttiled_small_t;
        ctrs::array<int, 5> min, max;
        _sp_hybrid int size(int i) const { return max[i] - min[i]; }
        _sp_hybrid std::size_t volume() const
        {
            std::size_t output = 1;
            algs::static_for<0, 5>([&](const auto& ii)
            {
                constexpr int i = ii.value;
                output *= size(i);
            });
            return output;
        }
        template <integral_scalar_or_array... is_t>
        _sp_hybrid _sp_inline std::size_t compute_offset(const is_t&... is) const noexcept
        {
            constexpr int tpow = 1;
            constexpr int tsz  = 1 << tpow;
            constexpr int mask = tsz - 1;
            
            int ni  = size(1);
            int nj  = size(2);
            int nk  = size(3);
            int nlb = size(4);
            int nv  = size(0);
            
            ni = ni >> tpow;
            nj = nj >> tpow;
            nk = nk >> tpow;
            
            int v  = get_val<0>(is...) - min[0];
            int i  = get_val<1>(is...) - min[1];
            int j  = get_val<2>(is...) - min[2];
            int k  = get_val<3>(is...) - min[3];
            int lb = get_val<4>(is...) - min[4];
            
            int ii = i & mask;
            i = i >> tpow;
            
            int jj = j & mask;
            j = j >> tpow;
            
            int kk = k & mask;
            k = k >> tpow;
            
            std::size_t  output = lb;
            output *= nk;
            output += k;
            output *= nj;
            output += j;
            output *= ni;
            output += i;
            output *= nv;
            output += v;
            output *= tsz;
            output += kk;
            output *= tsz;
            output += jj;
            output *= tsz;
            output += ii;
            
            return output;
        }
    };
    
    template <typename alias_t, typename device_t>
    inline auto make_grid_map(const tlinear_t&, const alias_t&, const device_t&,
        const ctrs::array<int, 2>& is,
        const ctrs::array<int, 2>& js,
        const ctrs::array<int, 2>& ks,
        const ctrs::array<int, 2>& ls)
    {
        if constexpr (ctrs::basic_array<alias_t>)
        {
            if constexpr (device::is_gpu<device_t>)
            {
                // return linear_t<5>{{is[0], js[0], ks[0], ls[0], 0}, {is[1], js[1], ks[1], ls[1], alias_t::size()}};
                return linear_t<5>{{0, is[0], js[0], ks[0], ls[0]}, {alias_t::size(), is[1], js[1], ks[1], ls[1]}};
            }
            else
            {
                return linear_t<5>{{0, is[0], js[0], ks[0], ls[0]}, {alias_t::size(), is[1], js[1], ks[1], ls[1]}};
            }
        }
        else
        {
            return linear_t<4>{{is[0], js[0], ks[0], ls[0]}, {is[1], js[1], ks[1], ls[1]}};
        }
    }
    
    template <typename alias_t, typename device_t>
    inline auto make_grid_map(const ttiled_t&, const alias_t&, const device_t&,
        const ctrs::array<int, 2>& is,
        const ctrs::array<int, 2>& js,
        const ctrs::array<int, 2>& ks,
        const ctrs::array<int, 2>& ls)
    {
        static_assert(ctrs::basic_array<alias_t>, "tiled memory map only available for array alias types");
        
        return tiled_t{{0, is[0], js[0], ks[0], ls[0]}, {alias_t::size(), is[1], js[1], ks[1], ls[1]}};
    }
    
    template <typename alias_t, typename device_t>
    inline auto make_grid_map(const ttiled_small_t&, const alias_t&, const device_t&,
        const ctrs::array<int, 2>& is,
        const ctrs::array<int, 2>& js,
        const ctrs::array<int, 2>& ks,
        const ctrs::array<int, 2>& ls)
    {
        static_assert(ctrs::basic_array<alias_t>, "tiled memory map only available for array alias types");
        
        return tiled_small_t{{0, is[0], js[0], ks[0], ls[0]}, {alias_t::size(), is[1], js[1], ks[1], ls[1]}};
    }
}