#include "spade.h"

int main(int argc, char** argv)
{
    print("hello from SPADE");

    print("Here's some random stuff");
    //spade provides some arithmetic statically-sized arrays
    using real_t = double;
    spade::ctrs::array<real_t, 2> ar0(0.1, 0.2);
    spade::ctrs::array<real_t, 2> ar1(0.2, 0.4);
    print(ar0);
    ar0 += ar1;
    print(ar1);

    //There isn't much to say here...

    return 0;
}
