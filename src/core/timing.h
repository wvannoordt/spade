#pragma once
#include <chrono>
#include <string>
namespace spade::timing
{
    class scoped_tmr_t
    {
        public:
            scoped_tmr_t(const std::string& name_in)
            {
                name = name_in;
                start = std::chrono::high_resolution_clock::now();
            }
            
            ~scoped_tmr_t()
            {
                end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                print(name, "took", duration.count()/1000.0, "ms.");
            }
            
        private:
            std::string name;
            decltype(std::chrono::high_resolution_clock::now()) start, end;
    };
    
    struct tmr_t
    {
        private:
            decltype(std::chrono::high_resolution_clock::now()) t_start, t_end;
        public:
            void start() { t_start = std::chrono::high_resolution_clock::now(); }
            void stop()  { t_end   = std::chrono::high_resolution_clock::now(); }
            
            //milliseconds
            double duration() const
            {
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start);
                return duration.count()/1000.0;
            }
    };
}