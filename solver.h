#ifndef SOLVER_H
#define SOLVER_H

#include "field.h"

class Simu
{
  public:
    Simu(int nx, int ny, double dt);
    void evol();
    void evol_wrong();

  private:
    void evol_inner(int i, int j);
    void evol_boundary();
    void evol_inject();
    void evel_inject_wrong();
    Field field;
    double const dx, dy, dt;
};

#endif
