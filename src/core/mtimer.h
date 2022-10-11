#pragma once

#include <chrono>

#include "typedef.h"

namespace spade::utils
{
    template <typename duration_t = real_t> struct mtimer_t
    {
        using time_t = decltype(std::chrono::steady_clock::now());
        
        std::vector<time_t> start_time, end_time;
        std::vector<duration_t> avg;
        std::vector<std::size_t> num_samples;
        
        mtimer_t()
        {
            start_time.resize(1);
            end_time.resize(1);
            avg.resize(1, 0.0);
            num_samples.resize(1, 0);
        }
        
        mtimer_t(const std::size_t& size)
        {
            start_time.resize(size);
            end_time.resize(size);
            avg.resize(size, 0.0);
            num_samples.resize(size, 0);
        }
        
        void start() { this->start(0); }
        void stop () { this->stop(0); }
        
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
            os << utils::pad_str(i, 5) << ":" << utils::pad_str(tmr.duration(i), 13) << "/" << utils::pad_str(tmr.average(i), 13) << "\n";
        }
        return os;
    }
}