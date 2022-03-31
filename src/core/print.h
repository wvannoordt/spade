#pragma once

#include <iostream>
#include <ostream>

template <typename T> static void print_recursive(std::ostream& stm, T t)
{
    stm << t << std::endl;
}

template <typename T, typename... Ts> static  void print_recursive(std::ostream& stm, T t, Ts... ts)
{
    stm << t << " ";
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