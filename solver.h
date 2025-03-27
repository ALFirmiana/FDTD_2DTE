#ifndef SOLVER_H
#define SOLVER_H

#include "field.h"
#include <cmath>

class Simu
{
  public:
    Simu(int nx, int ny, double dt, int total_step);
    void evol();
    void evol_wrong();
    int const total_step;
    Field field;

  private:
    void evol_Bz(int i, int j);
    void evol_Ex(int i, int j);
    void evol_Ey(int i, int j);

    void evol_boundary_left(int j, double Bz_right_old);
    void evol_boundary_right(int j, double Bz_left_old);
    void evol_boundary_up(int i, double Bz_down_old);
    void evol_boundary_down(int i, double Bz_up_old);
    void evol_boundary_angle_00(double Bz_right_up_old);
    void evol_boundary_angle_01(double Bz_left_up_old);
    void evol_boundary_angle_10(double Bz_right_down_old);
    void evol_boundary_angle_11(double Bz_left_down_old);

    void evol_inject();
    void evol_inject_wrong();

    double get_inject(double i, double j, double t);
    double const dt;
    double const dx = 1.0, dy = 1.0; // dx=dy=1 for Default, set this for maintenance
    int const nx, ny;
};

#endif
