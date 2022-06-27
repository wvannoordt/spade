#include "spade.h"
#include "prof_t.h"

int main(int argc, char** argv)
{    
    int ny = 100;
    int nt = 100000;
    
    std::vector<postprocessing::prof_t<double>*> reg;
    postprocessing::prof_t<double> yprof(ny, 0.0, "y", reg);
    postprocessing::prof_t<double> fprof(ny, 0.0, "f", reg);
    
    double y0 = 0.0;
    double y1 = 1.0;
    double t0 = 0.0;
    double t1 = 60.0;
    
    auto func1 = [&](const double& y, const double& t) -> double {return sin(2.0*spade::consts::pi*y*y*y + 2.0*spade::consts::pi*t);};
    auto func2 = [&](const double& y, const double& t) -> double {return sin(2.0*spade::consts::pi*t);};
    for (auto n: range(0,nt))
    {
        double t = t0 + n*(t1-t0)/(nt-1);
        for (auto i: range(0,ny))
        {
            double y = y0 + i*(y1-y0)/(ny-1);
            yprof.inst[i] = y;
            fprof.inst[i] = func1(y,t);
        }
        for (auto p:reg) p->aggregate();
    }
    
    std::ofstream myfile("soln.dat");
    for (auto i: range(0,ny))
    {
        myfile << yprof.avg[i] << "   " << fprof.avg[i] << "\n";
    }
    return 0;
}
