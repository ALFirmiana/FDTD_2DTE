#ifndef SOLVER_H
#define SOLVER_H

#include "field.h"

class Simu
{
  public:
    Simu(int nx, int ny, double dx, double dy, double dt);
    void evol();

  private:
    void evol_inner(int i, int j);
    void evol_boundary();
    void evol_inject();
    Field field;
    double const dx, dy, dt;
};

#endif
