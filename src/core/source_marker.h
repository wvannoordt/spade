#pragma once

#include <string>
#include <sstream>
#include <source_location>

namespace spade::utils
{
    struct source_marker_t
    {
#if (USE_SOURCE_LOCATION)
        std::source_location data;
        source_marker_t(const std::source_location location = std::source_location::current()) : data{location} {}
#endif
        
        std::string short_str() const
        {
            std::stringstream ss;
#if (USE_SOURCE_LOCATION)
            ss << "File: " << data.file_name() << "\n";
            ss << "Line: " << data.line();
#else
            ss << "(source location unsupported)";
#endif
            return ss.str();
        }
        
        std::string long_str() const
        {
            std::stringstream ss;
#if (USE_SOURCE_LOCATION)
            ss << "File:     " << data.file_name() << "\n";
            ss << "Line:     " << data.line() << "\n";
            ss << "Function: " << data.function_name();
#else
            ss << "(source location unsupported)";
#endif
            return ss.str();
        }
    };
}