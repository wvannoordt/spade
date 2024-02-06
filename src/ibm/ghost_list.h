#pragma once

#include "dispatch/shared_vector.h"
#include "grid/grid_index_types.h"

namespace spade::ibm
{
    template <typename float_t, template <typename> typename container_t = device::basic_shared_vector>
    struct ghost_list_t
    {
        using idx_t = grid::cell_idx_t;
        using pnt_t = coords::point_t<float_t>;
        using vec_t = ctrs::array<float_t, pnt_t::size()>;
        using image_type       = ghost_list_t<float_t, utils::vec_image_t>;
        using const_image_type = ghost_list_t<float_t, utils::const_vec_image_t>;
        
        //container_t is either std::vector<...> for CPU or device::device_vector<...> for GPU
        
        //Stores the points locally
        container_t<idx_t> indices;          // i, j, k, lb of the ghost points
        container_t<pnt_t> boundary_points;
        container_t<vec_t> boundary_normals;
        container_t<pnt_t> closest_points;
        container_t<vec_t> closest_normals;
        container_t<int>   directions;
        container_t<int>   signs;
        
        // Oh for to be able to use deducing this
        template <typename device_t>
        image_type image(const device_t& device)
        {
            image_type output;
            output.indices          = utils::make_vec_image(indices.data(device));
            output.boundary_points  = utils::make_vec_image(boundary_points.data(device));
            output.boundary_normals = utils::make_vec_image(boundary_normals.data(device));
            output.closest_points   = utils::make_vec_image(closest_points.data(device));
            output.closest_normals  = utils::make_vec_image(closest_normals.data(device));
            output.directions       = utils::make_vec_image(directions.data(device));
            output.signs            = utils::make_vec_image(signs.data(device));
            return output;
        }
        
        template <typename device_t>
        const_image_type image(const device_t& device) const
        {
            const_image_type output;
            output.indices          = utils::make_vec_image(indices.data(device));
            output.boundary_points  = utils::make_vec_image(boundary_points.data(device));
            output.boundary_normals = utils::make_vec_image(boundary_normals.data(device));
            output.closest_points   = utils::make_vec_image(closest_points.data(device));
            output.closest_normals  = utils::make_vec_image(closest_normals.data(device));
            output.directions       = utils::make_vec_image(directions.data(device));
            output.signs            = utils::make_vec_image(signs.data(device));
            return output;
        }
        
        void transfer()
        {
            //todo
            // if constexpr (device::has_device<container_t<int>>)
            // {
                indices.transfer();
                boundary_points.transfer();
                boundary_normals.transfer();
                closest_points.transfer();
                closest_normals.transfer();
                directions.transfer();
                signs.transfer();
            // }
        }
        
        void itransfer()
        {
            // if constexpr (device::has_device<container_t<int>>)
            // {
                indices.itransfer();
                boundary_points.itransfer();
                boundary_normals.itransfer();
                closest_points.itransfer();
                closest_normals.itransfer();
                directions.itransfer();
                signs.itransfer();
            // }
        }
    };
}