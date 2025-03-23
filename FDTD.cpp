#include "field.h"
#include "solver.h"
#include <iostream>

int main(int argc, char **argv)
{
    int dump_num = 10;
    Simu simulation(80, 60, 0.5, 100);
    std::cout << "simulation init" << std::endl;
    for (int dump_count = 0; simulation.field.getT() < simulation.total_step; simulation.field.push())
    {
        simulation.evol();
        dump_count++;
        std::cout << "step " << simulation.field.getT() << " done" << std::endl;
        if (dump_count == dump_num)
        {
            // TODO: simulation.field.writeToHDF5 not implement yet
            // simulation.field.writeToHDF5(const std::string &filename, double dx, double dy, double dt);
            dump_count = 0;
            std::cout << "step " << simulation.field.getT() << " dump done" << std::endl;
        }
    }
    return 0;
}
