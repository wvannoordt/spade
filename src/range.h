#pragma once

template <typename range_type> struct range_iterator_t
{
    range_type value;
    range_iterator_t(const range_type& value_in) { value = value_in; }
    bool operator != (const range_iterator_t& rhs)
    {
        return this->value < rhs.value;
    }
    range_iterator_t<range_type> operator ++ (void)
    {
        ++value;
        return *this;
    }
    range_type operator * (void) {return value;}
};

template <typename range_type> struct range_t
{
    range_type start_v, end_v;
    range_t(const range_type& start_in, const range_type& end_in)
    {
        start_v = start_in;
        end_v   = end_in;
    }
    range_iterator_t<range_type> begin(void) const noexcept
    {
        return range_iterator_t<range_type>(start_v);
    }
    range_iterator_t<range_type> end(void) const noexcept
    {
        return range_iterator_t<range_type>(end_v);
    }
};
template <typename range_type> range_t<range_type> range(const range_type& start, const range_type& end)
{
    return range_t<range_type>(start, end);
}