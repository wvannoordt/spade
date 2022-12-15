#pragma once

#include <type_traits>
#include <iostream>
#include <ostream>

namespace color
{
    enum color
    {
        black   = 30,
        red     = 31,
        green   = 32,
        yellow  = 33,
        blue    = 34,
        magenta = 35,
        cyan    = 36,
        white   = 37,
        none    = 0
    };
}

static std::string color_format(const color::color& c)
{
    return std::string("\033[0;" + std::to_string((int)c) + "m");
}

static const std::string dlm = " ";

namespace detail
{
    template <typename stream_t>
    static void stream_content(stream_t& stm, const color::color& content)
    {
        stm << color_format(content);
    }
    
    template <typename stream_t, typename content_t>
    static void stream_content(stream_t& stm, const content_t& content)
    {
        stm << content;
    }
    
    template <typename T> struct use_delim
    {
        static constexpr bool value = true;
    };
    template <> struct use_delim<color::color>
    {
        static constexpr bool value = false;
    };
    
    template <typename T, typename... Ts> struct any_color
    {
        static constexpr bool value = std::is_same<T, color::color>::value || any_color<Ts...>::value;
    };
    template <typename T> struct any_color<T>
    {
        static constexpr bool value = std::is_same<T, color::color>::value;
    };
}

template <typename T> static void print_recursive(std::ostream& stm, const T& t)
{
    detail::stream_content(stm, t);
    stm << std::endl;
}

template <typename T, typename... Ts> static  void print_recursive(std::ostream& stm, const T& t, const Ts&... ts)
{
    detail::stream_content(stm, t);
    if constexpr (detail::use_delim<T>::value)
    {
        detail::stream_content(stm, dlm);
    }
    print_recursive(stm, ts...);
}

template <typename T, typename... Ts> static void print(const T& t, const Ts&... ts)
{
    print_recursive(std::cout, t, ts...);
    if constexpr (detail::any_color<T, Ts...>::value) detail::stream_content(std::cout, color::none);
}

static void print(void)
{
    std::cout << std::endl;
}

static inline std::string zfill(int num, int numZeros)
{
    std::string output = std::to_string(num);
    while(output.length() < numZeros) output = "0" + output;
    return output;
}