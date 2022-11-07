#pragma once

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
        white   = 37//,
        // default = 0
    };
}

// static std::string color_format(const int& )

static const std::string dlm = " ";
template <typename T> static void print_recursive(std::ostream& stm, T t)
{
    stm << t << std::endl;
}

template <typename T, typename... Ts> static  void print_recursive(std::ostream& stm, T t, Ts... ts)
{
    stm << t << dlm;
    print_recursive(stm, ts...);
}

template <typename T, typename... Ts> static void print(const T& t, Ts... ts)
{
    print_recursive(std::cout, t, ts...);
}

static void print(void)
{
    std::cout << std::endl;
}

template <typename... Ts> static void PrintToStream(std::ostream& stm, Ts... ts)
{
    print_recursive(stm, ts...);
}

static inline std::string zfill(int num, int numZeros)
{
    std::string output = std::to_string(num);
    while(output.length() < numZeros) output = "0" + output;
    return output;
}