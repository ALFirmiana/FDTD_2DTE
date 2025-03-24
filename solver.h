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

    void evol_TSBC_Ex_left(int j, double Bz_inject_right);
    void evol_TSBC_Ex_right(int j, double Bz_inject_left);
    void evol_TSBC_Ey_up(int i, double Bz_inject_down);
    void evol_TSBC_Ey_down(int i, double Bz_inject_up);
    void evol_TSBC_Bz_left(int j, double Ey_inject_left);
    void evol_TSBC_Bz_right(int j, double Ey_inject_right);
    void evol_TSBC_Bz_up(int i, double Ex_inject_up);
    void evol_TSBC_Bz_down(int i, double Ex_inject_down);
    void evol_TSBC_Bz_00(double Ex_inject_down, double Ey_inject_left);  //左下边界
    void evol_TSBC_Bz_01(double Ex_inject_down, double Ey_inject_right); //右下边界
    void evol_TSBC_Bz_10(double Ex_inject_up, double Ey_inject_left);    //左上边界
    void evol_TSBC_Bz_11(double Ex_inject_up, double Ey_inject_right);  //右上边界


    void evol_inject_wrong();

    double get_inject(int i, int j, int t);
    double const dt;
    int const dx = 1, dy = 1; // dx=dy=1 for Default, set this for maintenance
    int const nx,ny;
};

#endif
