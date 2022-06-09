#pragma once
#include <iostream>
#include <array>

template <typename index_t, const std::size_t idx_dim> struct md_iterator
{
    constexpr size_t size(void) const noexcept {return idx_dim;}
    index_t idx_v, start_v, end_v;
    typedef std::size_t base_type_t;
    typename std::conditional<(idx_dim<=1), base_type_t, md_iterator<index_t,idx_dim-1>>::type next;
    template <typename... idxs_t> static void set_start_r(std::size_t& base){base = 0;}
    template <class callee_t> static void set_start_r(callee_t& arg)
    {
        arg.idx_v = arg.start_v;
        set_start_r(arg.next);
    }
    template <typename... idxs_t> static void set_end_r(std::size_t& base){base = 1;}
    template <class callee_t> static void set_end_r(callee_t& arg)
    {
        arg.idx_v = arg.start_v;
        set_end_r(arg.next);
    }
    md_iterator& set_start(void)
    {
        set_start_r(*this);
        return *this;
    }
    md_iterator& set_end(void)
    {
        set_end_r(*this);
        return *this;
    }
    template <typename... idxs_t> void call_set_val_r(std::size_t& base, idxs_t... idxs){base = 0;}
    template <class callee_t, typename... idxs_t> void call_set_val_r(callee_t& callee, idxs_t... idxs)
    {
        callee.set_val_r(idxs...);
    }
    template <typename idx_t, typename... idxs_t> void set_val_r(const idx_t& i0, const idx_t& i1, idxs_t... idxs)
    {
        static_assert(std::is_integral<idx_t>::value, "argumets to md_iterator constructor must be integral!");
        start_v = i0;
        idx_v   = i0;
        end_v   = i1;
        call_set_val_r(this->next, idxs...);
    }
    
    md_iterator(void) {}
    template <typename... idxs_t> md_iterator(idxs_t... idxs)
    {
        static_assert(sizeof...(idxs_t)==2*idx_dim, "md_iterator of size N requires 2N integral arguments!");
        set_val_r(idxs...);
    }
    
    md_iterator& operator ++ (void)
    {
        ++idx_v;
        if (idx_v>=end_v) { idx_v = start_v; ++next; }
        return *this;
    }
    
    bool operator != (const md_iterator& rhs)
    {
        return (this->idx_v != rhs.idx_v) || (this->next != rhs.next);
    }
    
    template <const std::size_t dim_i> static const md_iterator& deref_iter(const md_iterator<index_t,dim_i>* ptr, const index_t& val)
    {
        return *ptr;
    }
    static index_t deref_iter(const md_iterator<index_t,1>* ptr, const index_t& val)
    {
        return val;
    }
    
    auto operator * (void) const
    {
        return deref_iter(this, idx_v);
    }
    constexpr index_t brack_r(const typename md_iterator<index_t,0>::base_type_t& iter, std::size_t idx) const noexcept
    {
        return -1;
    }
    template <typename index_t_i, const std::size_t idx_dim_i>
    constexpr index_t brack_r(const md_iterator<index_t_i,idx_dim_i>& iter, std::size_t idx) const noexcept
    {
        return iter[idx];
    }
    constexpr index_t operator[] (const std::size_t idx) const noexcept
    {
        if (idx==0)
        {
            return idx_v;
        }
        else
        {
            return brack_r(next,idx-1);
        }
    }
};

template <typename index_t, const std::size_t ar_dim>
static void iter_ar_set_start_end(
    typename md_iterator<index_t, 0>::base_type_t& iter,
    const std::size_t idx,
    const std::array<index_t,ar_dim>& start_vs,
    const std::array<index_t,ar_dim>& end_vs)
{
    
}

template <typename index_t, const std::size_t idx_dim, const std::size_t ar_dim>
static void iter_ar_set_start_end(
    md_iterator<index_t, idx_dim>& iter,
    const std::size_t idx,
    const std::array<index_t,ar_dim>& start_vs,
    const std::array<index_t,ar_dim>& end_vs)
{
    iter.start_v = start_vs[idx];
    iter.end_v   = end_vs  [idx];
    iter_ar_set_start_end(iter.next, idx+1, start_vs, end_vs);
}

template <typename index_t, const std::size_t range_dim> struct md_range_t
{
    std::array<index_t,range_dim> start_v;
    std::array<index_t,range_dim> end_v;
    template <typename... init_t> md_range_t(void){}
    template <typename... init_t> md_range_t(init_t... inits)
    {
        set_r(0, inits...);
    }
    std::size_t size(void) const
    {
        std::size_t output = 1;
        for (std::size_t i = 0; i < range_dim; ++i)
        {
            output *= (end_v[i]>start_v[i])?(end_v[i]-start_v[i]):0;
        }
        return output;
    }
    template <typename init_t_1, typename init_t_2> void set_r(
        const std::size_t& idx, const init_t_1& imin, const init_t_2& imax)
    {
        static_assert(std::is_integral<init_t_1>::value && std::is_integral<init_t_2>::value, "Arguments to range must be integral");
        start_v[idx] = imin;
        end_v  [idx] = imax;
    }
    template <typename init_t_1, typename init_t_2, typename... init_t> void set_r(
        const std::size_t& idx, const init_t_1& imin, const init_t_2& imax,
        init_t... inits)
    {
        static_assert(std::is_integral<init_t_1>::value && std::is_integral<init_t_2>::value, "Arguments to range must be integral");
        start_v[idx] = imin;
        end_v  [idx] = imax;
        set_r(idx+1, inits...);
    }
    
    md_range_t<index_t, 1> subrange(const std::size_t& i) const
    {
        return md_range_t<index_t, 1>(this->start_v[i], this->end_v[i]);
    }
    
    template <typename index_r_t, const std::size_t range_dim_r>
    md_range_t<decltype(index_r_t()+index_t()), (range_dim+range_dim_r)>
    operator* (const md_range_t<index_r_t,range_dim_r>& rhs) const
    {
        md_range_t<decltype(index_r_t()+index_t()), (range_dim+range_dim_r)> output;
        for (std::size_t ii = 0; ii < range_dim; ++ii)
        {
            output.start_v[ii] = start_v[ii];
            output.end_v  [ii] = end_v  [ii];
        }
        for (std::size_t ii = 0; ii < range_dim_r; ++ii)
        {
            output.start_v[range_dim + ii] = rhs.start_v[ii];
            output.end_v  [range_dim + ii] = rhs.end_v  [ii];
        }
        return output;
    }
    
    auto begin(void) const
    {
        md_iterator<index_t,range_dim> output;
        iter_ar_set_start_end(output, 0, start_v, end_v);
        return output.set_start();
    }
    auto end(void) const
    {
        md_iterator<index_t,range_dim> output;
        iter_ar_set_start_end(output, 0, start_v, end_v);
        return output.set_end();
    }
};

template <typename prod> static auto gpt_r(const prod& p) { return p; }
template <typename prod, typename... prods> static  auto gpt_r(const prod& p, prods... ps) { return p*gpt_r(ps...); }
template <typename... prods> static auto get_prod(prods... ps) { return gpt_r(ps...); }
template <typename... idxs_t> static auto range(idxs_t... idxs)
{
    static_assert(2*(sizeof...(idxs_t)/2)==sizeof...(idxs_t), "mdrange requires an even number of integral arguments!");
    return md_range_t<decltype(get_prod(idxs...)),sizeof...(idxs_t)/2>(idxs...);
}
