#include "solver.h"
#include "field.h"
#include <cmath>

#define SQRT2 1.414214
#define PI 3.1415926

Simu::Simu(int nx, int ny, double dt, int total_step) : nx(nx), ny(ny), field(nx, ny), dt(dt), total_step(total_step)
{
}

double Simu::get_inject(int i, int j, int t)
{
    // TODO: need to change
    double omega = 2 * PI;
    return std::sin(omega * t);
}

void Simu::evol()
{
    int x0, y0, x1, y1;
    x0 = y0 = 10;
    x1 = nx - 10;
    y1 = ny - 10; // 表示总场边界位置的两个对角点

    std::vector<double> Ey_inject_left(y1 - y0 + 1, 0.0),
        Ey_inject_right(y1 - y0 + 1, 0.0), // 顺序均从下往上，分别用于更新左、右总场边界的Bz值
        Ex_inject_up(x1 - x0 + 1, 0.0),
        Ex_inject_down(x1 - x0 + 1, 0.0), // 顺序均从左往右，分别用于更新上、下总场边界的Bz值
        Bz_inject_left(y1 - y0 + 1, 0.0),
        Bz_inject_right(y1 - y0 + 1, 0.0), // 顺序均从下往上，分别用于更新右、左总场边界的Ex值
        Bz_inject_up(x1 - x0 + 1, 0.0),
        Bz_inject_down(x1 - x0 + 1, 0.0); // 顺序均从左往右，分别用于更新下、上总场边界的Ey值

    for (int j = y0; j < y1 + 1; j++)
    {
       // Bz_inject_left[j - y0] = get_inject(x0, j, field.getT());
    }

    std::vector<double> Bz_right_old(ny + 1), Bz_left_old(ny + 1), Bz_up_old(nx + 1), Bz_down_old(nx + 1);

    // init Bz_old
    for (int i = 0; i < nx + 1; i++)
    {
        Bz_down_old[i] = field.getBz(i, ny - 1);
        Bz_up_old[i] = field.getBz(i, 1);
    }
    for (int j = 0; j < ny + 1; j++)
    {
        Bz_left_old[j] = field.getBz(nx - 1, j);
        Bz_right_old[j] = field.getBz(1, j);
    }

    // evol Ex
    for (int i = 0; i < nx + 1; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            evol_Ex(i, j);
        }
    }
    for (int j = y0; j < y1 + 1; j++)
    {
        // evol_TSBC_Ex_left(j, Bz_inject_right[j - y0]);
        field.setEx(x0 - 1, j, field.getEx(x0 - 1, j) - (dt / dy) * Bz_inject_right[j - y0]);
    }
    for (int j = y0; j < y1 + 1; j++)
    {
        // evol_TSBC_Ex_right(j, Bz_inject_left[j - y0]);
        field.setEx(x1, j, field.getEx(x1, j) + (dt / dy) * Bz_inject_left[j - y0]);
    }

    // evol Ey
    for (int j = 0; j < ny + 1; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            evol_Ey(i, j);
        }
    }
    for (int i = x0; i < x1 + 1; i++)
    {
        // evol_TSBC_Ey_down(i, Bz_inject_up[i - x0]);
        field.setEy(i, y0 - 1, field.getEy(i, y0 - 1) + (dt / dx) * Bz_inject_up[i - x0]);
    }
    for (int i = x0; i < x1 + 1; i++)
    {
        // evol_TSBC_Ey_up(i, Bz_inject_down[i - x0]);
        field.setEy(i, y1, field.getEy(i, y1) - (dt / dx) * Bz_inject_down[i - x0]);
    }

    // evol Bz
    for (int i = 1; i < nx; i++)
    {
        for (int j = 1; j < ny; j++)
        {
            evol_Bz(i, j);
        }
    }
    // evol_TSBC_Bz_00(Ex_inject_down[0], Ey_inject_left[0]);
    field.setBz(x0, y0, field.getBz(x0, y0) + (dt / dy) * Ex_inject_down[0] - (dt / dx) * Ey_inject_left[0]);
    // evol_TSBC_Bz_01(Ex_inject_down[x1 - x0 + 1], Ey_inject_right[0]);
    field.setBz(x1, y0, field.getBz(x1, y0) + (dt / dy) * Ex_inject_down[x1 - x0 + 1] + (dt / dx) * Ey_inject_right[0]);
    // evol_TSBC_Bz_10(Ex_inject_up[0], Ey_inject_left[y1 - y0 + 1]);
    field.setBz(x0, y1, field.getBz(x0, y1) - (dt / dy) * Ex_inject_up[0] - (dt / dx) * Ey_inject_left[y1 - y0 + 1]);
    // evol_TSBC_Bz_11(Ex_inject_up[x1 - x0 + 1], Ey_inject_right[y1 - y0 + 1]);
    field.setBz(x1, y1,
                field.getBz(x1, y1) - (dt / dy) * Ex_inject_up[x1 - x0 + 1] + (dt / dx) * Ey_inject_right[y1 - y0 + 1]);

    for (int i = x0 + 1; i < x1; i++)
    { // evol_TSBC_Bz_down(i, Ex_inject_down[i - x0]);
        field.setBz(i, y0, field.getBz(i, y0) + (dt / dy) * Ex_inject_down[i - x0]);
        // evol_TSBC_Bz_up(i, Ex_inject_up[i - x0]);
        field.setBz(i, y1, field.getBz(i, y1) - (dt / dy) * Ex_inject_up[i - x0]);
    }
    for (int j = y0 + 1; j < y1; j++)
    {
        // evol_TSBC_Bz_left(j, Ey_inject_left[j - y0]);
        field.setBz(x0, j, field.getBz(x0, j) - (dt / dx) * Ey_inject_left[j - y0]);
        // evol_TSBC_Bz_right(j, Ey_inject_right[j - y0]);
        field.setBz(x1, j, field.getBz(x1, j) + (dt / dx) * Ey_inject_right[j - y0]);
    }

    for (int j = 1; j < ny; j++)
    {
        evol_boundary_left(j, Bz_right_old[j]);
        evol_boundary_right(j, Bz_left_old[j]);
    }
    for (int i = 1; i < nx; i++)
    {
        evol_boundary_up(i, Bz_down_old[i]);
        evol_boundary_down(i, Bz_up_old[i]);
    }
    evol_boundary_angle_00(Bz_right_old[1]);
    evol_boundary_angle_01(Bz_left_old[1]);
    evol_boundary_angle_10(Bz_down_old[1]);
    evol_boundary_angle_11(Bz_down_old[nx - 1]);
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

void Simu::evol_boundary_left(int j, double Bz_right_old)
{
    double Bz = field.getBz(0, j);
    double Bz_right = field.getBz(1, j);
    double Ex_up = field.getEy(0, j), Ex_down = field.getEy(0, j - 1);
    double Ex_right_up = field.getEy(1, j), Ex_right_down = field.getEy(1, j - 1);
    field.setBz(0, j,
                Bz_right_old + ((dt - dx) / (dt + dx)) * (Bz_right - Bz) +
                    (dx * dt / (2 * dy * (dt + dx))) * (Ex_up + Ex_right_up - Ex_down - Ex_right_down));
}

void Simu::evol_boundary_right(int j, double Bz_left_old)
{
    double Bz = field.getBz(nx, j);
    double Bz_left = field.getBz(nx - 1, j);
    double Ex_up = field.getEy(nx, j), Ex_down = field.getEy(nx, j - 1);
    double Ex_left_up = field.getEy(nx - 1, j), Ex_left_down = field.getEy(nx - 1, j - 1);
    field.setBz(nx, j,
                Bz_left_old + ((dt - dx) / (dt + dx)) * (Bz_left - Bz) +
                    (dx * dt / (2 * dy * (dt + dx))) * (Ex_up + Ex_left_up - Ex_down - Ex_left_down));
}

void Simu::evol_boundary_up(int i, double Bz_down_old)
{
    double Bz = field.getBz(i, ny);
    double Bz_down = field.getBz(i, ny - 1);
    double Ex_left = field.getEy(i - 1, ny), Ex_right = field.getEy(i, ny);
    double Ex_down_left = field.getEy(i - 1, ny - 1), Ex_down_right = field.getEy(i, ny - 1);
    field.setBz(i, ny,
                Bz_down_old + ((dt - dy) / (dt + dy)) * (Bz_down - Bz) +
                    (dy * dt / (2 * dx * (dt + dy))) * (Ex_left + Ex_down_left - Ex_right - Ex_down_right));
}

void Simu::evol_boundary_down(int i, double Bz_up_old)
{
    double Bz = field.getBz(i, 0);
    double Bz_up = field.getBz(i, 1);
    double Ex_left = field.getEy(i - 1, 0), Ex_right = field.getEy(i, 0);
    double Ex_up_left = field.getEy(i - 1, 1), Ex_up_right = field.getEy(i, 1);
    field.setBz(i, 0,
                Bz_up_old + ((dt - dy) / (dt + dy)) * (Bz_up - Bz) +
                    (dy * dt / (2 * dx * (dt + dy))) * (Ex_left + Ex_up_left - Ex_right - Ex_up_right));
}

void Simu::evol_boundary_angle_00(double Bz_right_up_old)
{
    double Bz = field.getBz(0, 0);
    double Bz_right_up = field.getBz(1, 1);
    field.setBz(0, 0, Bz_right_up_old + ((dt - SQRT2 * dx) / (dt + SQRT2 * dx)) * (Bz_right_up - Bz));
}

void Simu::evol_boundary_angle_01(double Bz_left_up_old)
{
    double Bz = field.getBz(nx, 0);
    double Bz_left_up = field.getBz(nx - 1, 1);
    field.setBz(nx, 0, Bz_left_up_old + ((dt - SQRT2 * dx) / (dt + SQRT2 * dx)) * (Bz_left_up - Bz));
}

void Simu::evol_boundary_angle_10(double Bz_right_down_old)
{
    double Bz = field.getBz(0, ny);
    double Bz_right_down = field.getBz(1, ny - 1);
    field.setBz(0, ny, Bz_right_down_old + ((dt - SQRT2 * dx) / (dt + SQRT2 * dx)) * (Bz_right_down - Bz));
}

void Simu::evol_boundary_angle_11(double Bz_left_down_old)
{
    double Bz = field.getBz(nx, ny);
    double Bz_left_down = field.getBz(nx - 1, ny - 1);
    field.setBz(nx, ny, Bz_left_down_old + ((dt - SQRT2 * dx) / (dt + SQRT2 * dx)) * (Bz_left_down - Bz));
}

// void Simu::evol_TSBC_Ex_left(int j, double Bz_inject_right){
//     field.setEx(x0 - 1, j, field.getEx(x0-1, j) - (dt/dy) * Bz_inject_right);
// }

// void Simu::evol_TSBC_Ex_right(int j, double Bz_inject_left){
//     field.setEx(x1, j, field.getEx(x1, j) + (dt/dy) * Bz_inject_left);
// }

// void Simu::evol_TSBC_Ey_down(int i, double Bz_inject_up){
//     field.setEy(i, y0-1, field.getEy(i, y0-1) + (dt/dx) * Bz_inject_up);
// }

// void Simu::evol_TSBC_Ey_up(int i, double Bz_inject_down){
//     field.setEy(i, y1, field.getEy(i, y1) - (dt/dx) * Bz_inject_down);
// }

// void Simu::evol_TSBC_Bz_left(int j, double Ey_inject_left){
//     field.setBz(x0, j, field.getBz(x0, j) - (dt/dx) * Ey_inject_left);
// }

// void Simu::evol_TSBC_Bz_right(int j, double Ey_inject_right){
//     field.setBz(x1, j, field.getBz(x1, j) + (dt/dx) * Ey_inject_right);
// }

// void Simu::evol_TSBC_Bz_up(int i, double Ex_inject_up){
//     field.setBz(i, y1, field.getBz(i, y1) - (dt/dy) * Ex_inject_up);
// }

// void Simu::evol_TSBC_Bz_down(int i, double Ex_inject_down){
//     field.setBz(i, y0, field.getBz(i, y0) + (dt/dy) * Ex_inject_down);
// }

// void Simu::evol_TSBC_Bz_00(double Ex_inject_down, double Ey_inject_left){
//     field.setBz(x0, y0, field.getBz(x0, y0) + (dt/dy) * Ex_inject_down - (dt/dx) * Ey_inject_left);
// }

// void Simu::evol_TSBC_Bz_01(double Ex_inject_down, double Ey_inject_right){
//     field.setBz(x1, y0, field.getBz(x1, y0) + (dt/dy) * Ex_inject_down + (dt/dx) * Ey_inject_right);
// }

// void Simu::evol_TSBC_Bz_10(double Ex_inject_up, double Ey_inject_left){
//     field.setBz(x0, y1, field.getBz(x0, y1) - (dt/dy) * Ex_inject_up - (dt/dx) * Ey_inject_left);
// }

// void Simu::evol_TSBC_Bz_11(double Ex_inject_up, double Ey_inject_right){
//     field.setBz(x1, y1, field.getBz(x1, y1) - (dt/dy) * Ex_inject_up + (dt/dx) * Ey_inject_right);
// }
