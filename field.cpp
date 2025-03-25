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

void Field::writeToHDF5(const std::string &filename) const
{
    // Check if the file exists
    bool file_exists = false;
    try
    {
        H5::H5File file(filename, H5F_ACC_RDONLY);
        file_exists = true;
        file.close();
    }
    catch (H5::FileIException &)
    {
        file_exists = false;
    }

    // Open the file in append mode if it exists, otherwise create it
    H5::H5File file;
    if (file_exists)
    {
        file.openFile(filename, H5F_ACC_RDWR);
    }
    else
    {
        file = H5::H5File(filename, H5F_ACC_TRUNC);
    }

    // Create dataset names using the current time step
    std::string Ex_name = "/Ex_" + std::to_string(t);
    std::string Ey_name = "/Ey_" + std::to_string(t);
    std::string Bz_name = "/Bz_" + std::to_string(t);

    // Define the dataspace for the arrays
    hsize_t dims_Ex[2] = {static_cast<hsize_t>(ny), static_cast<hsize_t>(nx + 1)};
    hsize_t dims_Ey[2] = {static_cast<hsize_t>(ny + 1), static_cast<hsize_t>(nx)};
    hsize_t dims_Bz[2] = {static_cast<hsize_t>(ny + 1), static_cast<hsize_t>(nx + 1)};

    H5::DataSpace dataspace_Ex(2, dims_Ex);
    H5::DataSpace dataspace_Ey(2, dims_Ey);
    H5::DataSpace dataspace_Bz(2, dims_Bz);

    // Create datasets and write data
    H5::DataSet dataset_Ex = file.createDataSet(Ex_name, H5::PredType::NATIVE_DOUBLE, dataspace_Ex);
    dataset_Ex.write(Ex.data(), H5::PredType::NATIVE_DOUBLE);

    H5::DataSet dataset_Ey = file.createDataSet(Ey_name, H5::PredType::NATIVE_DOUBLE, dataspace_Ey);
    dataset_Ey.write(Ey.data(), H5::PredType::NATIVE_DOUBLE);

    H5::DataSet dataset_Bz = file.createDataSet(Bz_name, H5::PredType::NATIVE_DOUBLE, dataspace_Bz);
    dataset_Bz.write(Bz.data(), H5::PredType::NATIVE_DOUBLE);

    // Close resources
    dataset_Ex.close();
    dataset_Ey.close();
    dataset_Bz.close();
    file.close();
}
