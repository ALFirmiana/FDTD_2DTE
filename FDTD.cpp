#include "field.h"
#include "solver.h"

int main(int argc, char **argv)
{
    int dump_num = 10;
    Simu simulation(80, 60, 0.5, 100);
    for (int dump_count = 0; simulation.field.getT() < simulation.total_step; simulation.field.push())
    {
        simulation.evol();
        dump_count++;
        if (dump_count == dump_num)
        {
            // TODO: simulation.field.writeToHDF5 not implement yet
            // simulation.field.writeToHDF5(const std::string &filename, double dx, double dy, double dt);
            dump_count = 0;
        }
    }
    return 0;
}
