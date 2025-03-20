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
    void evol_boundary(int i, int j);
    void evol_boundary_angle(int i, int j);
    void evol_inject_TSBC(int i, int j);
    void evol_inject_wrong();
    double get_inject(int i, int j, int t);
    double const dt;
    int const dx = 1, dy = 1; // dx=dy=1 for Default, set this for maintenance
};

#endif
