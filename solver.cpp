#include "solver.h"
#include "field.h"
#include <cmath>

Simu::Simu(int nx, int ny, double dt, int total_step) : field(nx, ny), dt(dt), total_step(total_step)
{
}

void Simu::evol()
{
    int i, j;
    double Bz_right_old[ny+1], Bz_left_old[ny+1], Bz_up_old[nx+1], Bz_down_old[nx+1];
    double Bz_old_bl[ny+1], Bz_old_br[ny+1], Bz_old_bu[nx+1], Bz_old_bd[nx+1];

    for ( i = 0; i <= nx; i++)
    {
        Bz_down_old[i] = field.getBz(i, ny-1);
        Bz_up_old[i] = field.getBz(i, 1);
        Bz_old_bd[i] = field.getBz(i, 0);
        Bz_old_bu[i] = field.getBz(i, ny);
    }

    for (j = 0; j <= ny; j++)
    {
        Bz_left_old[j] = field.getBz(nx-1, j);
        Bz_right_old[j] = field.getBz(1, j);
        Bz_old_bl[j] = field.getBz(0, j);
        Bz_old_br[j] = field.getBz(nx, j);
    }
    
    for (i = 1; i < nx; i++)
    {
        for (j = 1; j < ny; j++)
        {
            evol_Bz(i, j);
        }
    }
    for (i = 0; i <= nx; i++)
    {
        for (j = 0; j <= ny; j++)
        {
            evol_Ex(i, j);
            evol_Ey(i, j);
        }
    }
    for (j = 1; j < ny; j++)
    {
        evol_boundary_left(j, Bz_old_bl[j], Bz_right_old[j]);
        evol_boundary_right(j, Bz_old_br[j], Bz_left_old[j]);
    }
    for (i = 0; i < nx; i++)
    {
        evol_boundary_up(i, Bz_old_bu[i], Bz_down_old[i]);
        evol_boundary_down(i, Bz_old_bd[i], Bz_up_old[i]);
    }
    evol_boundary_angle_00(Bz_old_bl[0], Bz_right_old[1]);
    evol_boundary_angle_01(Bz_old_br[0], Bz_left_old[1]);
    evol_boundary_angle_10(Bz_old_bu[0], Bz_down_old[1]);
    evol_boundary_angle_11(Bz_old_bu[nx], Bz_down_old[nx-1]);
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

void evol_boundary_left(int j, double Bz_old, double Bz_right_old){
    double Bz = field.getBz(0, j);
    double Bz_right = field.getBz(1, j);
    double Ex_up = field.getEy(0, j), Ex_down = field.getEy(0, j-1);
    double Ex_right_up = field.getEy(1, j), Ex_right_down = field.getEy(1, j-1);
    field.setBz(0, j, Bz_right_old + ((dt-dx)/(dt+dx))*(Bz_right-Bz_old) + (dx*dt/(2*dy*(dt+dx)))*(Ex_up+Ex_right_up-Ex_down-Ex_right_down));
}

void evol_boundary_right(int j, double Bz_old, double Bz_left_old){
    double Bz = field.getBz(nx, j);
    double Bz_left = field.getBz(nx-1, j);
    double Ex_up = field.getEy(nx, j), Ex_down = field.getEy(nx, j-1);
    double Ex_left_up = field.getEy(nx-1, j), Ex_left_down = field.getEy(nx-1, j-1);
    field.setBz(nx, j, Bz_left_old + ((dt-dx)/(dt+dx))*(Bz_left-Bz_old) + (dx*dt/(2*dy*(dt+dx)))*(Ex_up+Ex_left_up-Ex_down-Ex_left_down));
}

void evol_boundary_up(int i, double Bz_old, double Bz_down_old){
    double Bz = field.getBz(i, ny);
    double Bz_down = field.getBz(i, ny-1);
    double Ex_left = field.getEy(i-1, ny), Ex_right = field.getEy(i, ny);
    double Ex_down_left = field.getEy(i-1, ny-1), Ex_down_right = field.getEy(i, ny-1);
    field.setBz(i, ny, Bz_down_old + ((dt-dy)/(dt+dy))*(Bz_down-Bz_old) + (dy*dt/(2*dx*(dt+dy)))*(Ex_left+Ex_down_left-Ex_right-Ex_down_right));
}

void evol_boundary_down(int i, double Bz_old, double Bz_up_old){
    double Bz = field.getBz(i, 0);
    double Bz_up = field.getBz(i, 1);
    double Ex_left = field.getEy(i-1, 0), Ex_right = field.getEy(i, 0);
    double Ex_up_left = field.getEy(i-1, 1), Ex_up_right = field.getEy(i, 1);
    field.setBz(i, 0, Bz_up_old + ((dt-dy)/(dt+dy))*(Bz_up-Bz_old) + (dy*dt/(2*dx*(dt+dy)))*(Ex_left+Ex_up_left-Ex_right-Ex_up_right));
}

void evol_boundary_angle_00(double Bz_old, double Bz_right_up_old){
    double Bz = field.getBz(0, 0);
    double Bz_right_up = field.getBz(1, 1);
    field.setBz(0, 0, Bz_right_up_old + ((dt-sqrt(2.0)*dx)/(dt+sqrt(2.0)*dx))*(Bz_right_up-Bz_old));
}


void evol_boundary_angle_01(double Bz_old, double Bz_left_up_old){
    double Bz = field.getBz(nx, 0);
    double Bz_left_up = field.getBz(nx-1, 1);
    field.setBz(nx, 0, Bz_left_up_old + ((dt-sqrt(2.0)*dx)/(dt+sqrt(2.0)*dx))*(Bz_left_up-Bz_old));
}

void evol_boundary_angle_10(double Bz_old, double Bz_right_down_old){
    double Bz = field.getBz(0, ny);
    double Bz_right_down = field.getBz(1, ny-1);
    field.setBz(0, ny, Bz_right_down_old + ((dt-sqrt(2.0)*dx)/(dt+sqrt(2.0)*dx))*(Bz_right_down-Bz_old));
}

void evol_boundary_angle_11(double Bz_old, double Bz_left_down_old){
    double Bz = field.getBz(nx, ny);
    double Bz_left_down = field.getBz(nx-1, ny-1);
    field.setBz(nx, ny, Bz_left_down_old + ((dt-sqrt(2.0)*dx)/(dt+sqrt(2.0)*dx))*(Bz_left_down-Bz_old));
}

double Simu::get_inject(int i, int j, int t)
{
    // TODO: need to change
    double omega = 1;
    return std::sin(omega * t);
}
