//
// Created by david on 8/1/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_HDF5_H
#define FINITE_DMRG_EIGEN_CLASS_HDF5_H

//#include <hdf5.h>
//#include <hdf5_hl.h>
#include <iostream>
#include <n_tensor_extra.h>
#include <experimental/filesystem>
//#include <H5CommonFG.h>
//#include <H5File.h>
//#include <H5Cpp.h>
#include <hdf5.h>
namespace fs = std::experimental::filesystem::v1;

using namespace Textra;



/*!
 \brief Writes and reads data to a binary hdf5-file.

 # HDF5 Class


*/

class class_hdf5 {
private:
    hid_t       file;          //
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
        file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC,  H5P_DEFAULT, H5P_DEFAULT);
    }

    void open_file(){
        file = H5Fopen (file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }

    void create_group(const std::string group_relative_name){
        hid_t   lcpl    = H5Pcreate(H5P_LINK_CREATE);   //Create missing intermediate group if they don't exist
        herr_t  status  = H5Pset_create_intermediate_group(lcpl, 1);
        hid_t   group   = H5Gcreate(file,group_relative_name.c_str(), lcpl,H5P_DEFAULT,H5P_DEFAULT);
        H5Pclose(lcpl);
        H5Gclose(group);

    }

    template <typename DataType>
    void write_to_hdf5(const DataType &data, const std::string dataset_relative_name){
        hid_t dataspace = get_DataSpace(data);
        hid_t datatype  = get_DataType(data);

        hid_t lcpl      = H5Pcreate(H5P_LINK_CREATE);   //Create missing intermediate group if they don't exist
        herr_t status1  = H5Pset_create_intermediate_group(lcpl, 1);
        hid_t dataset   = H5Dcreate(file,
                                    dataset_relative_name.c_str(),
                                    datatype,
                                    dataspace,
                                    lcpl,
                                    H5P_DEFAULT,
                                    H5P_DEFAULT);
        herr_t status2 = H5Dwrite(dataset,
                                 datatype,
                                 H5S_ALL,
                                 H5S_ALL,
                                 H5P_DEFAULT,
                                 data.data());
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Dclose(dataset);
        H5Pclose(lcpl);

    }

    ~class_hdf5(){
        H5Fclose(file);
    }

private:

    template<typename ObjectType>
    auto get_DataType(ObjectType &object){
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

    template<typename TensorType>
    hid_t get_DataSpace(const TensorType &tensor){
        int rank = tensor.rank();
        hsize_t     dims[rank];
        std::copy(tensor.dimensions().begin(), tensor.dimensions().end(), dims);
        return H5Screate_simple(rank, dims, nullptr);
    }

};


// Specializations
template<>
hid_t class_hdf5::get_DataSpace<MatrixType>(const MatrixType &matrix);
#endif //FINITE_DMRG_EIGEN_CLASS_HDF5_H
