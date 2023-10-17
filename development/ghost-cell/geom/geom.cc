#include <print.h>
#include "geolytical.h"

int main(int argc, char** argv)
{
	double rad = 0.05;
    bbox bounds;
    bounds.xmin = -rad;
    bounds.xmax =  rad;
    bounds.ymin = -rad;
    bounds.ymax =  rad;
    bounds.zmin = -rad;
    bounds.zmax =  rad;
    geolytical::Sphere ball(256, bounds);
    ball.OutputToVtk("sphere.vtk");
    return 0;
}
