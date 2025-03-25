#include "field.h"
#include "solver.h"
#include <iostream>
#include <string>

int main(int argc, char **argv)
{
    int dump_num = 10;
    Simu simulation(80, 60, 0.5, 100);
    std::string const output_file = "test.h5";
    std::cout << "simulation init" << std::endl;
    for (int dump_count = 0; simulation.field.getT() < simulation.total_step; simulation.field.push())
    {
        simulation.evol();
        dump_count++;
        std::cout << "step " << simulation.field.getT() << " done" << std::endl;
        if (dump_count == dump_num)
        {
            simulation.field.writeToHDF5(output_file);
            dump_count = 0;
            std::cout << "step " << simulation.field.getT() << " dump done" << std::endl;
        }
    }
    return 0;
}
