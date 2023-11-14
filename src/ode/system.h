#pragma once

namespace spade::ode
{
    template <typename data_t, typename mesh_t, typename variables_t>
    struct system_t
    {
        using fundamental_type = data_t;
        using value_type       = data_t;
        using mesh_type        = mesh_t;
        using variable_list    = variables_t;
        
        mesh_type mesh;
        
        system_t(const data_t&, const mesh_t& mesh_in, const variables_t&)
        : mesh{mesh_in}
        {
            
        }
    };
}