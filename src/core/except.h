#pragma once

#include <stdexcept>
#include <string>

#include "core/point.h"

namespace spade::except
{
    struct sp_exception : public std::exception
    {
        std::string message;
        sp_exception(const std::string& message_in) : message{message_in} {}
        const char* what() const throw()
        {
            return message.c_str();
        }
    };
    
    template <typename float_t>
    struct points_exception : public sp_exception
    {
        std::vector<coords::point_t<float_t>> data;
        points_exception(const std::string& message_in, const std::vector<coords::point_t<float_t>>& data_in) : sp_exception(message_in), data{data_in} {}
    };
}