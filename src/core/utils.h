#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>
#include <iostream>
#include <cstdlib>

namespace spade::utils
{
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
    static auto reduce_over_params(const reduction_t& reduce_op, const invokable_t& func, const param_t& param, params_t... params)
    {
        auto total = func(param);
        reduce_op(total, reduce_over_params(reduce_op, func, params...));
        return total;
    }
    template <class invokable_t, class param_t, class... params_t>
    static auto reduce_over_params(const invokable_t& func, const param_t& param, params_t... params)
    {
        typedef decltype(func(param)) ret_type;
        return reduce_over_params(sum_reduce<ret_type,ret_type>, func, param, params...);
    }
    
    template <class invokable_t, class param_t>
    static void foreach_param(const invokable_t& func, const param_t& param)
    {
        func(param);
    }
    
    template <class invokable_t, class param_t, class... params_t>
    static void foreach_param(const invokable_t& func, const param_t& param, params_t... params)
    {
        func(param);
        foreach_param(params...);
    }
    
    
    
    template <typename tp1_t, typename tp2_t> static auto max(const tp1_t& t1, const tp2_t& t2)
    {
        return t1<t2?t2:t1;
    }
    
    template <typename tp_t, typename... tps_t> static auto max(const tp_t& t, tps_t... ts)
    {
        return max(t, max(ts...));
    }
    
    template <typename tp1_t, typename tp2_t> static auto min(const tp1_t& t1, const tp2_t& t2)
    {
        return t1<t2?t1:t2;
    }
    
    template <typename tp_t, typename... tps_t> static auto min(const tp_t& t, tps_t... ts)
    {
        return min(t, min(ts...));
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
    
    static inline bool is_big_endian()
    {
        int num = 1;
        return (! ( *(char *)&num == 1 ));
    }
    
    static inline void random_seed(unsigned int seed)
    {
        srand(seed);
    }
    
    template <typename out_t=double> static inline out_t unitary_random(void)
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
    template <typename data_t, const padding::padding pd = padding::right> static inline std::string pad_str(const data_t& data, const std::size_t& pad_size, const char& pad_char = ' ')
    {
        std::string output = std::to_string(data);
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
}
