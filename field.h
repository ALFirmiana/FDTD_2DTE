#ifndef FIELD_H
#define FIELD_H

#include <vector>

class Field
{
  public:
    Field(int nx, int ny);

  private:
    int nx, ny;
    std::vector<double> Ex;
    std::vector<double> Ey;
    std::vector<double> Bz;
};

#endif // !FIELD_H
