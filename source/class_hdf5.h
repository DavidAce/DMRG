//
// Created by david on 8/1/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_HDF5_H
#define FINITE_DMRG_EIGEN_CLASS_HDF5_H

#include <hdf5.h>
#include <hdf5_hl.h>
#include <iostream>
#include <n_tensor_extra.h>
using namespace Textra;



/*!
 \brief Writes and reads data to a binary hdf5-file.

 # HDF5 Class

   
*/

class class_hdf5 {
public:
    hid_t       file_id;
    std::string file_name = "dmrg_output.h5";
    class_hdf5(std::string file_name_ =  "dmrg_output.h5"):file_name(file_name_){
        file_id = H5Fcreate (file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }


    ~class_hdf5(){
        H5Fclose (file_id);
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
        std::cout << " Tensor : Rank " << rank << " dims: " << dims[0] << dims[1] << dims[2] << dims[3] << "\n" << tensor << std::endl;
        H5LTmake_dataset(file_id,dataset_name.c_str(), rank, dims,scalar,tensor.data());
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
