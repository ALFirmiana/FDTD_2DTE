#include "field.h"
#include "params.h"
#include "solver.h"
#include <cmath>
#include <iostream>

int main(int argc, char **argv)
{
    Simu simulation(NX, NY, DT, TOTLE_STEP);
    std::cout << "simulation init" << std::endl;

    simulation.field.writeToText();
    while (simulation.field.getT() < simulation.total_step)
    {
        // simulation.field.setBz(75, 50, std::sin(simulation.field.getT()));
        simulation.evol();
        std::cout << "step " << simulation.field.getT() << " done" << std::endl;

        if (simulation.field.getT() % DUMP_NUM == 0)
        {
            simulation.field.writeToText();
            std::cout << "step " << simulation.field.getT() << " dump done" << std::endl;
        }
    }

    return 0;
}
