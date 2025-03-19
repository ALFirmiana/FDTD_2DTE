#include "solver.h"
#include "field.h"

Simu::Simu(int nx, int ny, double dx, double dy, double dt) : field(nx, ny), dx(dx), dy(dy), dt(dt)
{
}

void Simu::evol()
{
}
