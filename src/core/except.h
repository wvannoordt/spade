#pragma once

#include <stdexcept>
#include <string>

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
}