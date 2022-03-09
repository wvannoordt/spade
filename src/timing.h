#pragma once
#include <chrono>
#include <string>
#include "cmf.h"
using cmf::print;
namespace timing
{
    class scoped_tmr_t
    {
        public:
            scoped_tmr_t(const std::string& name_in)
            {
                name = name_in;
                start = std::chrono::high_resolution_clock::now();
            }
            
            ~scoped_tmr_t(void)
            {
                end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
                print(name, "took", duration.count(), "ms.");
            }
            
        private:
            std::string name;
            decltype(std::chrono::high_resolution_clock::now()) start, end;
    };
}