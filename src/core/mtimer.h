#pragma once

#include <chrono>

#include "typedef.h"

namespace spade::timing
{
    template <typename duration_t = real_t> struct mtimer_t
    {
        using time_t = decltype(std::chrono::steady_clock::now());
        
        std::vector<time_t> start_time, end_time;
        std::vector<duration_t> avg;
        std::vector<std::size_t> num_samples;
        std::vector<std::string> names;
        std::map<std::string, int> name_map;
        
        mtimer_t()
        {
            start_time.resize(1);
            end_time.resize(1);
            avg.resize(1, 0.0);
            num_samples.resize(1, 0);
            default_init_names();
            create_name_map();
        }
        
        mtimer_t(const std::size_t& size)
        {
            start_time.resize(size);
            end_time.resize(size);
            avg.resize(size, 0.0);
            num_samples.resize(size, 0);
            default_init_names();
            create_name_map();
        }
        
        void rsetnames(){}
        template <typename... strings_t> void rsetnames(const std::string& name_in, const strings_t&... names_in) {names.push_back(name_in); rsetnames(names_in...);}
        
        template <typename... strings_t> mtimer_t(const std::string& name_in, const strings_t&... names_in)
        {
            const std::size_t size = 1 + sizeof...(strings_t);
            start_time.resize(size);
            end_time.resize(size);
            avg.resize(size, 0.0);
            num_samples.resize(size, 0);
            names.clear();
            names.push_back(name_in);
            rsetnames(names_in...);
            create_name_map();
        }
        
        void default_init_names()
        {
            names.clear();
            for (auto i: range(0, start_time.size())) names.push_back(std::to_string(i));
        }
        
        void create_name_map()
        {
            name_map.clear();
            for (auto i: range(0, start_time.size()))
            {
                name_map.insert({this->names[i], i});
            }
        }
        
        void start() { this->start(0); }
        void stop () { this->stop(0); }
        
        void start(const std::string& name_in) { this->start(name_map[name_in]); }
        void stop (const std::string& name_in) { this->stop(name_map[name_in]); }
        
        void start(const std::size_t& idx)
        {
            start_time[idx] = std::chrono::steady_clock::now();
        }
        
        void stop(const std::size_t& idx)
        {
        	end_time[idx] = std::chrono::steady_clock::now();
            auto dur = this->duration(idx);
            auto alpha = duration_t(num_samples[idx])/(num_samples[idx]+1);
            avg[idx] = alpha*avg[idx] + (1.0-alpha)*dur;
            num_samples[idx]++;
        }
        
        
        duration_t duration() const { return this->duration(0); }
        duration_t average () const { return this->average(0); }
        
        duration_t duration(const std::size_t& idx) const
        {
            auto diff = end_time[idx] - start_time[idx];
            return (duration_t)(std::chrono::duration_cast<std::chrono::milliseconds>(diff).count())/duration_t(1000.0);
        }
        
        duration_t average(const std::size_t& idx) const { return avg[idx]; }
        
        std::size_t size() const {return num_samples.size();}
    };
    
    template <typename duration_t> static std::ostream& operator << (std::ostream& os, const mtimer_t<duration_t>& tmr)
    {
        os << "Timer: elapsed/average (seconds)\n";
        for (int i = 0; i < tmr.size(); ++i)
        {
            os << utils::pad_str(tmr.names[i], 15) << ":" << utils::pad_str(tmr.duration(i), 13) << "/" << utils::pad_str(tmr.average(i), 13) << "\n";
        }
        return os;
    }
}