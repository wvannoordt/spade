#pragma once

#include <concepts>
#include <type_traits>

namespace spade::algs
{    
    namespace detail
    {
        template <typename item_t, typename... items_t> struct transform_list
        {
            const item_t* item;
            transform_list<items_t...> list;
            void set_items(const item_t& item_in, const items_t&... items_in)
            {
                item = &item_in;
                list.set_items(items_in...);
            }
            
            template <typename operand_t> void apply_forward_transforms(operand_t& q) const
            {
                item->transform_forward(q);
                list.apply_forward_transforms(q);
            }
            
            template <typename operand_t> void apply_inverse_transforms(operand_t& q) const
            {
                list.apply_inverse_transforms(q);
                item->transform_inverse(q);
            }
        };
        template <typename item_t> struct transform_list<item_t>
        {
            const item_t* item;
            void set_items(const item_t& item_in)
            {
                item = &item_in;
            }
            
            template <typename operand_t> void apply_forward_transforms(operand_t& q) const
            {
                item->transform_forward(q);
            }
            
            template <typename operand_t> void apply_inverse_transforms(operand_t& q) const
            {
                item->transform_inverse(q);
            }
        };
    }
    
    template <typename... transforms_t> struct composite_transform_t
    {
        detail::transform_list<transforms_t...> list;
        composite_transform_t(){}
        composite_transform_t(const transforms_t&... transforms)
        {
            list.set_items(transforms...);
        }
        
        template <typename operand_t> void transform_forward(operand_t& q) const
        {
            list.apply_forward_transforms(q);
        }
        
        template <typename operand_t> void transform_inverse(operand_t& q) const
        {
            list.apply_inverse_transforms(q);
        }
    };
}