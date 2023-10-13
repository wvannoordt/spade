#include <print.h>
#include "geolytical.h"

int main(int argc, char** argv)
{
    bbox bounds;
    bounds.xmin = -1.0;
    bounds.xmax =  1.0;
    bounds.ymin = -1.0;
    bounds.ymax =  1.0;
    bounds.zmin = -1.0;
    bounds.zmax =  1.0;
    geolytical::Sphere ball(256, bounds);
    ball.OutputToVtk("../sphere.vtk");
    return 0;
}
