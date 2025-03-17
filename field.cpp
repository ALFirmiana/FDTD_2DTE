#include "field.h"

#include <stdexcept>
#include <vector>

Field::Field(int nx, int ny) : nx(nx), ny(ny), Ex(nx * ny, 0.0), Ey(nx * ny, 0.0), Bz(nx * ny, 0.0)
{
}

double Field::getEx(int i, int j) const
{
    return Ex[i * ny + j];
}

void Field::checkBounds(int i, int j) const
{
    if (i < 0 || i >= nx || j < 0 || j >= ny)
    {
        throw std::out_of_range("Index out of bounds");
    }
}
