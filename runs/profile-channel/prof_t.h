#pragma once
#include <iostream>
#include <string>
namespace postprocessing
{
    template <typename data_t> struct prof_t
    {
        prof_t(const int& ny, const data_t& def, const std::string& name_in, std::vector<prof_t<data_t>*>& registry_in)
        {
            inst.resize(ny, def);
            avg.resize(ny, def);
            registry = &registry_in;
            registry->push_back(this);
            num_avg = 0;
        }
        void aggregate(void)
        {
            data_t alpha = (data_t)num_avg/(num_avg+1);
            data_t beta  = 1.0/(num_avg+1);
            for (int i = 0; i < this->size(); ++i)
            {
                avg[i] = alpha*avg[i] + beta*inst[i];
            }
            num_avg++;
        }
        std::vector<data_t> inst, avg;
        std::size_t size(void) const { return inst.size(); }
        std::vector<prof_t<data_t>*>* registry;
        std::string name;
        int num_avg;
    };
}