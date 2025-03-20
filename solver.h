#ifndef SOLVER_H
#define SOLVER_H

#include "field.h"

class Simu
{
  public:
    Simu(int nx, int ny, double dt, int total_step);
    void evol();
    void evol_wrong();

  private:
    void evol_Bz(int i, int j);
    void evol_Ex(int i, int j);
    void evol_Ey(int i, int j);
    void evol_boundary();
    void evol_inject();
    void evel_inject_wrong();
    Field field;
    double const dt;
    int const dx = 1, dy = 1; // dx=dy=1 for Default, set this for maintenance
    int const total_step;
};

#endif
