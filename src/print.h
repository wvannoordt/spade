#pragma once

#include <iostream>
#include <ostream>

template <typename T> void print_recursive(std::ostream& stm, T t)
{
    stm << t << std::endl;
}

template <typename T, typename... Ts> void print_recursive(std::ostream& stm, T t, Ts... ts)
{
    stm << t << " ";
    print_recursive(stm, ts...);
}

template <typename... Ts> void print(Ts... ts)
{
    print_recursive(std::cout, ts...);
}

template <typename... Ts> void PrintToStream(std::ostream& stm, Ts... ts)
{
    print_recursive(stm, ts...);
}

static inline std::string zfill(int num, int numZeros)
{
    std::string output = std::to_string(num);
    while(output.length() < numZeros) output = "0" + output;
    return output;
}