#ifndef FIELD_H
#define FIELD_H

#include <vector>

class Field
{
  public:
    Field(int nx, int ny);
    double getEx(int i, int j) const;
    double getEy(int i, int j) const;
    double getBz(int i, int j) const;
    void setEx(int i, int j, double new_Ex);
    void setEy(int i, int j, double new_Ey);
    void setBz(int i, int j, double new_Bz);

  private:
    const int nx, ny;
    std::vector<double> Ex, Ey, Bz;
    void checkBounds(int i, int j) const;
};

#endif // !FIELD_H
