#include "solver.h"
#include "field.h"
#include <cmath>

Simu::Simu(int nx, int ny, double dt, int total_step) : field(nx, ny), dt(dt), total_step(total_step)
{
}

void Simu::evol()
{
}

void Simu::evol_Bz(int i, int j)
{
    double Bz = field.getBz(i, j);
    double Ey_left = field.getEy(i - 1, j), Ey_right = field.getEy(i, j);
    double Ex_up = field.getEx(i, j), Ex_down = field.getEx(i, j - 1);
    field.setBz(i, j, Bz - dt / dx * (Ey_right - Ey_left) + dt / dy * (Ex_up - Ex_down));
}

void Simu::evol_Ex(int i, int j)
{
    double Ex = field.getEx(i, j);
    double Bz_up = field.getBz(i, j + 1), Bz_down = field.getBz(i, j);
    field.setEx(i, j, Ex + dt / dy * (Bz_up - Bz_down));
}

void Simu::evol_Ey(int i, int j)
{
    double Ey = field.getEy(i, j);
    double Bz_left = field.getBz(i, j), Bz_right = field.getBz(i + 1, j);
    field.setEy(i, j, Ey - dt / dx * (Bz_left - Bz_right));
}

double Simu::get_inject(int i, int j, int t)
{
    // TODO: need to change
    double omega = 1;
    return std::sin(omega * t);
}
