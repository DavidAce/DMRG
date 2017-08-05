//
// Created by david on 8/1/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_HDF5_H
#define FINITE_DMRG_EIGEN_CLASS_HDF5_H

#include <hdf5.h>
#include <hdf5_hl.h>
#include <iostream>
#include <n_tensor_extra.h>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem::v1;

using namespace Textra;



/*!
 \brief Writes and reads data to a binary hdf5-file.

 # HDF5 Class


*/

class class_hdf5 {
private:
    hid_t       file_id;
    fs::path    executable_path    = fs::current_path();
    fs::path    output_path;
    fs::path    output_dirname;
    fs::path    file_name;
    bool create_dir;
    void set_output_path();

public:

    class_hdf5(fs::path file_name_ = fs::path("dmrg_output.h5"), fs::path output_dirname_ = fs::path("output"), bool create_dir_ = true):
            file_name(file_name_),
            output_dirname(output_dirname_),
            create_dir(create_dir_)
    {
        set_output_path();
        file_id = H5Fcreate (file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

    void open_file(){
        file_id = H5Fopen (file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }


    template <typename TensorType>
    void write_to_hdf5(const TensorType &tensor, const std::string dataset_name){
        int rank = tensor.rank();
        hsize_t     dims[rank];
        hid_t scalar = get_ScalarType(tensor);
        std::copy(tensor.dimensions().begin(), tensor.dimensions().end(), dims);
        H5LTmake_dataset(file_id,dataset_name.c_str(), rank, dims,scalar,tensor.data());
    }

    ~class_hdf5(){
        H5Fclose (file_id);
    }

private:

    template<typename ObjectType>
    auto get_ScalarType(ObjectType &object){
        if(std::is_same<typename ObjectType::Scalar, double>::value){
            return H5Tcopy(H5T_NATIVE_DOUBLE);
        }
        else
        if(std::is_same<typename ObjectType::Scalar, int>::value ){
            return H5Tcopy(H5T_NATIVE_INT);
        }else{
            std::cout << "Error: unsupported Scalar type for HDF5." <<std::endl;
            exit(1);
        }
    }



};


// Specializations
template<>
void class_hdf5::write_to_hdf5<MatrixType>(const MatrixType &matrix, const std::string dataset_name);

#endif //FINITE_DMRG_EIGEN_CLASS_HDF5_H
