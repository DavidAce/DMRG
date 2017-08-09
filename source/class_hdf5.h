//
// Created by david on 8/1/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_HDF5_H
#define FINITE_DMRG_EIGEN_CLASS_HDF5_H

//#include <hdf5.h>
//#include <hdf5_hl.h>
#include <iostream>
#include <n_tensor_extra.h>
#include <funcs.h>

#include <experimental/filesystem>
#include <experimental/type_traits>

//#include <H5CommonFG.h>
//#include <H5File.h>
#include <H5Cpp.h>
//#include <hdf5.h>
#include <type_traits>
#include <funcs.h>
namespace fs = std::experimental::filesystem::v1;
using namespace Textra;



namespace TypeCheck{
    template <typename T> using Data_t          = decltype(std::declval<T>().data());
    template <typename T> using Size_t          = decltype(std::declval<T>().size());
    template <typename T> using Dims_t          = decltype(std::declval<T>().dimensions());
    template <typename T> using Setz_t          = decltype(std::declval<T>().setZero());
    template <typename T> using Matr_t          = decltype(std::declval<T>().matrix());
    template <typename T> using Push_t          = decltype(std::declval<T>().push_back());

    template <typename T> using has_member_data         = std::experimental::is_detected<Data_t, T>;
    template <typename T> using has_member_size         = std::experimental::is_detected<Size_t, T>;
    template <typename T> using has_member_dimensions   = std::experimental::is_detected<Dims_t , T>;
    template <typename T> using is_eigen                = std::experimental::is_detected<Setz_t , T>;
    template <typename T> using is_tensor               = std::experimental::is_detected<Dims_t , T>;
    template <typename T> using is_matrix               = std::experimental::is_detected<Matr_t, T>;
    template <typename T> using is_stdvector            = std::experimental::is_detected<Push_t, T>;
}
namespace tc = TypeCheck;



/*!
 \brief Writes and reads data to a binary hdf5-file.

 # HDF5 Class


*/

class class_hdf5 {
private:
    hid_t       file;          //
    hid_t       group, lcpl;
    fs::path    file_name;
    fs::path    executable_path    = fs::current_path();
    fs::path    output_path;
    fs::path    output_dirname;
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
        Tensor3 test;
        MatrixType test2;
//        test.setR
    }
    ~class_hdf5(){
        H5Fclose(file);
    }


    void open_file(){
        file = H5Fopen (file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        Tensor3 test;
    }

    void create_group_and_open(const std::string group_relative_name){
        lcpl            = H5Pcreate(H5P_LINK_CREATE);   //Create missing intermediate group if they don't exist
        herr_t status   = H5Pset_create_intermediate_group(lcpl, 1);
        group           = H5Gcreate(file,group_relative_name.c_str(), lcpl,H5P_DEFAULT,H5P_DEFAULT);
        H5Eclose_stack(status);
    }


    void close_group(const std::string group_relative_name){
        H5Pclose(lcpl);
        H5Gclose(group);
    }

    template <typename DataType>
    void create_group_w_attribute(const std::string group_relative_name, const std::string attribute_name, const DataType &data){
        lcpl                 = H5Pcreate(H5P_LINK_CREATE);   //Create missing intermediate group if they don't exist
        herr_t status1       = H5Pset_create_intermediate_group(lcpl, 1);
        group                = H5Gcreate(file,group_relative_name.c_str(), lcpl,H5P_DEFAULT,H5P_DEFAULT);
        hid_t datatype       = get_DataType<DataType>();
        hid_t dataspace      = get_DataSpace(data);
        hid_t attribute      = H5Acreate(group, attribute_name.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT );
        herr_t status2       = H5Awrite(attribute, datatype, &data);

        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Aclose(attribute);
        H5Pclose(lcpl);
        H5Gclose(group);
        H5Eclose_stack(status1);
        H5Eclose_stack(status2);

    }


//    template <typename DataType, typename = typename std::enable_if<std::is_base_of<Eigen::TensorBase<double>,DataType>::value>::type>
    template <typename DataType, typename std::enable_if<tc::has_member_data<DataType>::value>::type* = nullptr>
    void write_to_file(const DataType &data, const std::string dataset_relative_name){
//        std::cout << tc::has_member_dimensions<DataType>::value << std::endl; // false
        hid_t dataspace = get_DataSpace(data);
        hid_t datatype  = get_DataType<DataType>();
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
        H5Eclose_stack(status1);
        H5Eclose_stack(status2);
    }


    //Overloads will wrap the input data inside an Eigen-Type.
    template <typename T, typename std::enable_if<std::is_arithmetic<T>::value>::type>
    void write_to_file(const T &data, const std::string dataset_relative_name){
        Eigen::TensorMap<Eigen::Tensor<T,1>> wrapper(*data,1);
        write_to_file(wrapper, dataset_relative_name);
    };
    //Overloads will wrap the input data inside an Eigen-Type.
    template<typename T>
    void write_to_file(const std::vector<T> &data, const std::string dataset_relative_name){
        Eigen::TensorMap<Eigen::Tensor<T,1>> wrapper(data.data(),data.size());
        write_to_file(wrapper, dataset_relative_name);
    }

private:

    template<typename DataType, typename std::enable_if<std::is_arithmetic<DataType>::value>::type* = nullptr>
    hid_t get_DataType(){
        if(std::is_same<DataType, int>::value)          {return  H5Tcopy(H5T_NATIVE_INT);}
        if(std::is_same<DataType, long>::value)         {return  H5Tcopy(H5T_NATIVE_LONG);}
        if(std::is_same<DataType, double>::value)       {return  H5Tcopy(H5T_NATIVE_DOUBLE);}
        if(std::is_same<DataType, float>::value)       {return  H5Tcopy(H5T_NATIVE_FLOAT);}
        std::cout << "get_DataType can't match the type provided" << std::endl;
        exit(1);
    }

    template<typename DataType, typename  std::enable_if<tc::is_eigen<DataType>::value>::type* = nullptr>
    hid_t get_DataType(){return  get_DataType<typename DataType::Scalar>();}

    template<typename DataType, typename  std::enable_if<tc::is_stdvector<DataType>::value>::type* = nullptr>
    hid_t get_DataType(){return  get_DataType<typename DataType::value_type>();}


    template<typename DataType, typename  std::enable_if<tc::is_tensor<DataType>::value>::type* = nullptr>
    hid_t get_DataSpace(const DataType &data){
        int rank = data.rank();
        hsize_t dims[rank];
        std::copy(data.dimensions().begin(), data.dimensions().end(), dims);
        return H5Screate_simple(rank, dims, nullptr);
    }

    template<typename DataType, typename  std::enable_if<std::is_arithmetic<DataType>::value>::type* = nullptr>
    hid_t get_DataSpace(const DataType &data) {
        int rank = 1;
        hsize_t dims[rank] = {1};
        return H5Screate_simple(rank, dims, nullptr);;
    }

    template<typename DataType, typename  std::enable_if<tc::is_matrix<DataType>::value>::type* =nullptr>
    hid_t get_DataSpace(const DataType &data) {
        int rank = 2;
        hsize_t dims[rank] = {(hsize_t)data.rows(),(hsize_t) data.cols()};
        return H5Screate_simple(rank, dims, nullptr);
    }

    template<typename DataType, typename  std::enable_if<tc::is_stdvector<DataType>::value>::type* =nullptr>
    hid_t get_DataSpace(const DataType &data) {
        int rank = 1;
        hsize_t dims[rank] = {1};
        return H5Screate_simple(rank, dims, nullptr);
    }

};

#endif //FINITE_DMRG_EIGEN_CLASS_HDF5_H
