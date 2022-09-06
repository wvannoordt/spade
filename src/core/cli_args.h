#pragma once
#include <string>
#include <vector>
#include <map>
#include <sstream>
namespace spade::cli_args
{
    struct arg_t
    {
        arg_t(const std::string& raw_in) {raw = raw_in;}
        template <typename fundamental_t> operator fundamental_t() const
        {
            fundamental_t p;
            std::stringstream iss(this->raw);
            iss >> p;
            return p;
        }
        
        std::string raw;
    };
    
    struct shortname_args_t
    {
        shortname_args_t(int argc, char** argv)
        {
            for (std::size_t i = 0; i < argc; ++i)
            {
                std::string p(argv[i]);
                if (p.length()>0) raw_args.push_back(p);
            }
            
            for (std::size_t i = 0; i < argc; ++i)
            {
                std::string key = raw_args[i];
                if ((key[0] == '-') && (i <= raw_args.size()-2))
                {
                    if (!this->has_arg(key)) data.insert({key, raw_args[i+1]});
                }
            }
        }
        bool has_arg(const std::string& key)
        {
            return !(data.find(key) == data.end());
        }
        std::string operator[] (const std::string& key)
        {
            return data[key];
        }
        std::vector<std::string> raw_args;
        std::map<std::string, std::string> data;
    };
}