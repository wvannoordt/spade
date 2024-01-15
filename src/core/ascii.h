#pragma once

#include <sstream>
#include <fstream>
#include <string>

#include "core/except.h"

namespace spade::utils
{
    struct ascii_file_t
    {
        std::ifstream fh;
        std::string p_line{""};
        std::istringstream iss;
        std::string fn;
        std::size_t line_num = 0;
        ascii_file_t(const std::string& fname) : fh{fname}, fn{fname}
        {
            if (!fh) throw except::sp_exception("cannot open file \"" + fname + "\"");
        }
        
        const std::string& next_line()
        {
            ++line_num;
            if (!std::getline(fh, p_line)) throw except::sp_exception("unexpected end of file \"" + fn + "\", line " + std::to_string(line_num) + ".");
            iss.clear();
            iss.str(p_line);
            return p_line;
        }

        const std::string& line() const
        {
            return p_line;
        }
        
        void expect(const std::string& content)
        {
            if (content != p_line)
            {
                throw except::sp_exception(
                    "Error reading \"" + fn + "\": expecting \""
                    + content + "\" on line " + std::to_string(line_num) + ", but found \"" + p_line + "\"");
            }
        }
        
        template <typename data_t>
        void tparse_sing(std::istringstream& iss, data_t& data)
        {
            iss >> data;
            if (iss.fail())
            {
                throw except::sp_exception(
                    "Error reading \"" + fn + "\": attempting parse on line " + std::to_string(line_num) + ", but failed: \"" + p_line + "\"");
            }
        }
        
        template <typename data_t>
        void r_parse(std::istringstream& iss, data_t& data)
        {
            tparse_sing(iss, data);
        }
        
        template <ctrs::basic_array data_t>
        void r_parse(std::istringstream& iss, data_t& data)
        {
            for (int i = 0; i < data.size(); ++i)
            {
                tparse_sing(iss, data[i]);
            }
        }
        
        template <typename data_t, typename... datas_t>
        void r_parse(std::istringstream& iss, data_t& data, datas_t&... datas)
        {
            tparse_sing(iss, data);
            r_parse(iss, datas...);
        }
        
        template <typename... datas_t>
        void parse(datas_t&... datas)
        {
            r_parse(iss, datas...);
        }

        template <typename data_t>
        bool try_parse_sing(std::istringstream& iss, data_t& data)
        {
            iss >> data;
            return !iss.fail();
        }

        template <typename data_t>
        bool r_try_parse(std::istringstream& iss, data_t& data)
        {
            return try_parse_sing(iss, data);
        }
        
        template <typename data_t, typename... datas_t>
        bool r_try_parse(std::istringstream& iss, data_t& data, datas_t&... datas)
        {
            if (!try_parse_sing(iss, data)) return false;
            return r_try_parse(iss, datas...);
        }

        template <typename... datas_t>
        bool try_parse(datas_t&... datas)
        {
            return r_try_parse(iss, datas...);
        }
    };
}