#include "spade.h"

template <typename in_t> void print_addr(const in_t& in)
{
    print(&in);
}
int main(int argc, char** argv)
{
    struct trans1_t
    {
        void transform_forward(double& d) const {d = d*d;}
        void transform_inverse(double& d) const {d = std::sqrt(d);}
    } trans1;
    
    struct trans2_t
    {
        void transform_forward(double& d) const {d = 3.0*d;}
        void transform_inverse(double& d) const {d = d/3.0;}
    } trans2;
    
    double z = 2.0;
    print("z = ", z);
    trans1.transform_forward(z);
    trans2.transform_forward(z);
    print("z = ", z);
    trans2.transform_inverse(z);
    trans1.transform_inverse(z);
    print("z = ", z);
    
    print("=================");
    spade::algs::composite_transform_t composite_transform(trans1, trans2);
    print("z = ", z);
    composite_transform.transform_forward(z);
    print("z = ", z);
    composite_transform.transform_inverse(z);
    print("z = ", z);
}
