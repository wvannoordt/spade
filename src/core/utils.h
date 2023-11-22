#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <type_traits>

namespace spade::utils
{
    template<class T, class U=
        typename std::remove_cv<
        typename std::remove_pointer<
        typename std::remove_reference<
        typename std::remove_extent<
        T
        >::type
        >::type
        >::type
        >::type
        > struct remove_all : remove_all<U> {};
    template<class T> struct remove_all<T, T> { using type = T; };
    
    
    namespace detail
    {
        template <typename T> struct istdv_t { constexpr static bool value = false;};
        template <typename T, typename A> struct istdv_t<std::vector<T,A>> { constexpr static bool value = true;};
    }
    template <typename T> concept is_std_vector = detail::istdv_t<T>::value;

    template <class ltype, class rtype> static void sum_reduce(ltype& total, const rtype& r)
    {
        total += r;
    }
    
    template <class reduction_t, class invokable_t, class param_t>
    static auto reduce_over_params(const reduction_t& reduce, const invokable_t& func, const param_t& param)
    {
        return func(param);
    }
    
    template <class reduction_t, class invokable_t, class param_t, class... params_t>
    static auto reduce_over_params(const reduction_t& reduce_op, const invokable_t& func, const param_t& param, const params_t&... params)
    {
        auto total = func(param);
        reduce_op(total, reduce_over_params(reduce_op, func, params...));
        return total;
    }
    template <class invokable_t, class param_t, class... params_t>
    static auto reduce_over_params(const invokable_t& func, const param_t& param, const params_t&... params)
    {
        typedef decltype(func(param)) ret_type;
        return reduce_over_params(sum_reduce<ret_type,ret_type>, func, param, params...);
    }
    
    template <class invokable_t, class param_t>
    _sp_hybrid static void foreach_param(const invokable_t& func, const param_t& param)
    {
        func(param);
    }
    
    template <class invokable_t, class param_t, class... params_t>
    _sp_hybrid static void foreach_param(const invokable_t& func, const param_t& param, const params_t&... params)
    {
        func(param);
        foreach_param(func, params...);
    }
    
    template <typename tp1_t> _sp_hybrid [[nodiscard]] constexpr static auto max(const tp1_t& t1)
    {
        return t1;
    }
    
    template <typename tp1_t, typename tp2_t> _sp_hybrid [[nodiscard]] constexpr static auto max(const tp1_t& t1, const tp2_t& t2)
    {
        return t1<t2?t2:t1;
    }
    
    template <typename tp_t, typename... tps_t> _sp_hybrid [[nodiscard]] constexpr static auto max(const tp_t& t, const tps_t&... ts)
    {
        return max(t, max(ts...));
    }
    
    template <typename tp1_t> _sp_hybrid [[nodiscard]] constexpr static auto min(const tp1_t& t1)
    {
        return t1;
    }
    
    template <typename tp1_t, typename tp2_t> _sp_hybrid [[nodiscard]] constexpr static auto min(const tp1_t& t1, const tp2_t& t2)
    {
        return t1<t2?t1:t2;
    }
    
    template <typename tp_t, typename... tps_t> _sp_hybrid [[nodiscard]] constexpr static auto min(const tp_t& t, const tps_t&... ts)
    {
        return min(t, min(ts...));
    }
    
    template <typename rtype> _sp_hybrid [[nodiscard]] constexpr static auto abs(const rtype& val)
    {
        return val<0?-val:val;
    }
    
    template <typename rtype> _sp_hybrid [[nodiscard]] constexpr static auto sign(const rtype& val)
    {
        return val<0?-1:1;
    }
    
    static inline void get_format_substrings(std::vector<std::string>& subStrings, const std::string& templateStr)
    {
        std::string delimiter = "{}";
        std::string templateStrCopy = templateStr;
        size_t pos = 0;
        std::string token;
        while ((pos = templateStrCopy.find(delimiter)) != std::string::npos)
        {
            token = templateStrCopy.substr(0, pos);
            subStrings.push_back(token);
            templateStrCopy.erase(0, pos + delimiter.length());
        }
        subStrings.push_back(templateStrCopy);
    }

    static inline std::vector<std::string> string_split(std::string templateStr, std::string delimiter)
    {
        std::vector<std::string> subStrings;
        std::string templateStrCopy = templateStr;
        size_t pos = 0;
        std::string token;
        while ((pos = templateStrCopy.find(delimiter)) != std::string::npos)
        {
            token = templateStrCopy.substr(0, pos);
            subStrings.push_back(token);
            templateStrCopy.erase(0, pos + delimiter.length());
        }
        subStrings.push_back(templateStrCopy);
        return subStrings;
    }
    
    template <typename T> static inline void strformat_recursive (const std::string& templateStr, std::ostringstream& strstream, std::vector<std::string>& subStrings, int& lev, const T& t)
    {
        strstream << subStrings[lev] << t;
        lev++;
    }
    
    template <typename T, typename... Ts> static inline void strformat_recursive (const std::string& templateStr, std::ostringstream& strstream, std::vector<std::string>& subStrings, int& lev, const T& t, Ts... ts)
    {
        strstream << subStrings[lev] << t;
        lev++;
        strformat_recursive(templateStr, strstream, subStrings, lev, ts...);
    }
    
    static inline std::string strformat (std::string templateStr) { return templateStr; }
    
    template <typename T, typename... Ts> static inline std::string strformat (const std::string& templateStr, const T& t, Ts... ts)
    {
        std::ostringstream strstream;
        std::vector<std::string> subStrings;
        get_format_substrings(subStrings, templateStr);
        if (((sizeof...(Ts))+1)!=(subStrings.size()-1))
        {
            std::string err_string = "strformat invocation with template string \"" + templateStr +"\" and wrong number of arguments (" +
                std::to_string((sizeof...(Ts))+1) + " args instead of " + std::to_string(subStrings.size()-1) + ")";
            throw std::runtime_error(err_string.c_str());
        }
        int lev = 0;
        strformat_recursive(templateStr, strstream, subStrings, lev, t, ts...);
        strstream << subStrings[lev];
        return strstream.str();
    }
    
    template <class data_t> static std::string zfill(const data_t& data, const std::size_t num_zeros)
    {
        std::string output = std::to_string(data);
        while (output.length() < num_zeros) output = "0"+output;
        return output;
    }
    
    template <typename data_t>
    static std::string join(const data_t& data, const std::string& join_dlm)
    {
        std::string output = "";
        std::stringstream ss(output);
        int ct = 0;
        for (const auto& i: data)
        {
            if (ct != 0) ss << join_dlm;
            ss << i;
            ++ct;
        }
        return ss.str();
    }
    
    static inline bool is_big_endian()
    {
        int num = 1;
        return (! ( *(char *)&num == 1 ));
    }
    
    static inline void random_seed(unsigned int seed)
    {
        srand(seed);
    }
    template <typename out_t=double> static inline out_t unitary_random()
    {
        return ((out_t) rand() / (RAND_MAX));
    }
    
    namespace padding
    {
        enum padding
        {
            left, right
        };
    }
    
    static std::string to_string(const std::string& str)
    {
        return str;
    }
    
    template <typename str_t>
    static std::string to_string(const str_t& str)
    {
        return std::to_string(str);
    }
    
    template <typename data_t, const padding::padding pd = padding::right>
    static inline std::string pad_str(const data_t& data, const std::size_t& pad_size, const char& pad_char = ' ')
    {
        std::string output = to_string(data);
        if constexpr (pd == padding::left)
        {
            while (output.length() < pad_size) output = pad_char + output;
        }
        else
        {
            while (output.length() < pad_size) output = output + pad_char;
        }
        return output;
    }
    
    template <const std::size_t position, const std::size_t idx, typename type_t, typename... types_t>
    struct get_pack_type_helper 
    {
        typedef std::conditional<position==idx, type_t, typename get_pack_type_helper<position+1, idx, types_t...>::type>::type type;
    };
    
    template <const std::size_t position, const std::size_t idx, typename type_t>
    struct get_pack_type_helper<position, idx, type_t>
    {
        typedef type_t type;
    };
    
    
    template <const std::size_t idx, typename... types_t>
    requires (idx < sizeof...(types_t))
    struct get_pack_type
    {
        typedef get_pack_type_helper<0,idx,types_t...>::type type;
    };
    
    template <const std::size_t idx, typename... types_t>
    using pack_element = typename get_pack_type<idx, types_t...>::type;
    
    namespace detail
    {        
        template <typename T, typename... Ts>
        struct index_of_pack_impl;
        
        template <typename T, typename... Ts>
        struct index_of_pack_impl<T, T, Ts...>
        {
            constexpr static int value = 0;
        };
        
        template <typename T, typename U, typename... Ts>
        struct index_of_pack_impl<T, U, Ts...>
        {
            constexpr static int value = 1 + index_of_pack_impl<T, Ts...>::value;
        };
    }
    
    template <typename T, typename... Ts>
    constexpr static int index_of_pack = detail::index_of_pack_impl<T, Ts...>::value;
    
    // Used to store a reallocation-agnostic
    // reference to a vector element
    template <typename data_t>
    struct vector_location
    {
        std::vector<data_t>* container = nullptr;
        std::size_t offset;
        vector_location(std::vector<data_t>& container_in, const std::size_t offset_in)
        : container{&container_in}, offset{offset_in} {}
        vector_location(std::vector<data_t>* container_in, const std::size_t offset_in)
        : container{container_in},  offset{offset_in} {}
        vector_location& operator = (const vector_location& rhs)
        {
            offset = rhs.offset;
            container = rhs.container;
            return *this;
        }
        
        bool is_null()      const { return container == nullptr; }
        data_t& get()
        {
            auto& v = (*container);
            return v[offset];
        }
        
        const data_t& get() const
        {
            const auto& v = (*container);
            return v[offset];
        }
        
        constexpr static vector_location null() { return vector_location(nullptr, 0); }
    };
    
    template <typename T, typename... Ts>
    struct is_contained;
    
    template <typename T, typename U, typename... Vs>
    struct is_contained<T, U, Vs...>
    {
        constexpr static bool value = std::is_same<T, U>::value || is_contained<T, Vs...>::value;
    };
    
    template <typename T>
    struct is_contained<T>
    {
        constexpr static bool value = false;
    };
    
    template <typename... List>
    struct is_unique_impl;
    
    template <typename U, typename... Vs>
    struct is_unique_impl<U, Vs...>
    {
        constexpr static bool value = !is_contained<U, Vs...>::value && is_unique_impl<Vs...>::value;
    };
    
    template <>
    struct is_unique_impl<>
    {
        constexpr static bool value = true;
    };
    
    template <typename... Ts> concept is_unique = is_unique_impl<Ts...>::value;
    
    template <typename T> using raw_ptr_t       = T*;
    template <typename T> using const_raw_ptr_t = const T*;
}
