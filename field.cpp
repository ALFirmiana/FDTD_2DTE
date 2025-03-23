#include "field.h"

#include <H5Cpp.h>
#include <vector>

Field::Field(int nx, int ny)
    : nx(nx), ny(ny), t(0), Ex((nx + 1) * ny, 0.0), Ey(nx * (ny + 1), 0.0), Bz((nx + 1) * (ny + 1), 0.0)
{
}

// Bz设置在格点上，i表示x坐标，j表示y坐标；Bz每行有(nx+1)个元素
// Ex设置在“y轴”上的半格点上，(i,j)表示物理上的(i,j+1/2)；Ex每行有(nx+1)个元素。
// Ey设置在“x轴”上的半格点上，(i,j)表示物理上的(i+1/2,j)；Ey每行有nx个元素。
double Field::getEx(int i, int j) const
{
    return Ex[j * (nx + 1) + i];
}

double Field::getEy(int i, int j) const
{
    return Ey[j * nx + i];
}

double Field::getBz(int i, int j) const
{
    return Bz[j * (nx + 1) + i];
}

int Field::getT() const
{
    return t;
}

void Field::setEx(int i, int j, double new_Ex)
{
    Ex[j * (nx + 1) + i] = new_Ex;
}

void Field::setEy(int i, int j, double new_Ey)
{
    Ey[j * nx + i] = new_Ey;
}

void Field::setBz(int i, int j, double new_Bz)
{
    Bz[j * (nx + 1) + i] = new_Bz;
}

void Field::push()
{
    t++;
}

// TODO: generate by deepseek, need to check
// void Field::writeToHDF5(const std::string &filename, double dx, double dy, double dt) const
// {
//     // 创建 HDF5 文件
//     H5::H5File file(filename, H5F_ACC_TRUNC);
//
//     // 创建根组
//     H5::Group root = file.openGroup("/");
//
//     // 写入元数据
//     H5::Attribute dx_attr = root.createAttribute("dx", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
//     dx_attr.write(H5::PredType::NATIVE_DOUBLE, &dx);
//
//     H5::Attribute dy_attr = root.createAttribute("dy", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
//     dy_attr.write(H5::PredType::NATIVE_DOUBLE, &dy);
//
//     H5::Attribute dt_attr = root.createAttribute("dt", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
//     dt_attr.write(H5::PredType::NATIVE_DOUBLE, &dt);
//
//     // 创建网格数据集
//     hsize_t dims[2] = {static_cast<hsize_t>(nx), static_cast<hsize_t>(ny)};
//     H5::DataSpace dataspace(2, dims);
//
//     // 写入电场数据
//     H5::DataSet Ex_dataset = file.createDataSet("/Ex", H5::PredType::NATIVE_DOUBLE, dataspace);
//     Ex_dataset.write(Ex.data(), H5::PredType::NATIVE_DOUBLE);
//
//     // 写入磁场数据
//     H5::DataSet Bz_dataset = file.createDataSet("/Bz", H5::PredType::NATIVE_DOUBLE, dataspace);
//     Bz_dataset.write(Bz.data(), H5::PredType::NATIVE_DOUBLE);
//
//     // 关闭文件
//     file.close();
// }
