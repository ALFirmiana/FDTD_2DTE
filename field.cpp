#include "field.h"

#include <H5Cpp.h>
#include <stdexcept>
#include <vector>

Field::Field(int nx, int ny) : nx(nx), ny(ny), Ex(nx * ny, 0.0), Ey(nx * ny, 0.0), Bz(nx * ny, 0.0)
{
}

double Field::getEx(int i, int j) const
{
    checkBounds(i, j);
    return Ex[i * ny + j];
}

double Field::getEy(int i, int j) const
{
    checkBounds(i, j);
    return Ey[i * ny + j];
}

double Field::getBz(int i, int j) const
{
    checkBounds(i, j);
    return Bz[i * ny + j];
}

void Field::setEx(int i, int j, double new_Ex)
{
    checkBounds(i, j);
    Ex[i * ny + j] = new_Ex;
}

void Field::setEy(int i, int j, double new_Ex)
{
    checkBounds(i, j);
    Ey[i * ny + j] = new_Ex;
}

void Field::setBz(int i, int j, double new_Ex)
{
    checkBounds(i, j);
    Bz[i * ny + j] = new_Ex;
}

void Field::checkBounds(int i, int j) const
{
    if (i < 0 || i >= nx || j < 0 || j >= ny)
    {
        throw std::out_of_range("Index out of bounds");
    }
}

// TODO: generate by deepseek, need to check
void Field::writeToHDF5(const std::string &filename, double dx, double dy, double dt) const
{
    // 创建 HDF5 文件
    H5::H5File file(filename, H5F_ACC_TRUNC);

    // 创建根组
    H5::Group root = file.openGroup("/");

    // 写入元数据
    H5::Attribute dx_attr = root.createAttribute("dx", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
    dx_attr.write(H5::PredType::NATIVE_DOUBLE, &dx);

    H5::Attribute dy_attr = root.createAttribute("dy", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
    dy_attr.write(H5::PredType::NATIVE_DOUBLE, &dy);

    H5::Attribute dt_attr = root.createAttribute("dt", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
    dt_attr.write(H5::PredType::NATIVE_DOUBLE, &dt);

    // 创建网格数据集
    hsize_t dims[2] = {static_cast<hsize_t>(nx), static_cast<hsize_t>(ny)};
    H5::DataSpace dataspace(2, dims);

    // 写入电场数据
    H5::DataSet Ex_dataset = file.createDataSet("/Ex", H5::PredType::NATIVE_DOUBLE, dataspace);
    Ex_dataset.write(Ex.data(), H5::PredType::NATIVE_DOUBLE);

    // 写入磁场数据
    H5::DataSet Bz_dataset = file.createDataSet("/Bz", H5::PredType::NATIVE_DOUBLE, dataspace);
    Bz_dataset.write(Bz.data(), H5::PredType::NATIVE_DOUBLE);

    // 关闭文件
    file.close();
}
