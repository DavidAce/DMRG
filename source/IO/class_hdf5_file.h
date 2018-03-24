//
// Created by david on 8/1/17.
//

#ifndef CLASS_HDF5_H
#define CLASS_HDF5_H

#include <hdf5.h>
#include <hdf5_hl.h>
#include <iostream>
#include <string>
#include <iomanip>

#include <experimental/filesystem>
#include <experimental/type_traits>

#include <type_traits>
#include <typeinfo>
#include "class_custom_cout.h"
#include <general/nmspc_type_check.h>
#include <directory.h>
#include <gitversion.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <iterator>


namespace fs = std::experimental::filesystem::v1;
namespace tc = TypeCheck;

/*!
 \brief Writes and reads data to a binary hdf5-file.
*/

class class_hdf5_file {
private:
    herr_t      retval;
    fs::path    output_filename;
    fs::path    output_folder;
    fs::path    output_file_path;
    bool        create_dir;
    bool        overwrite;
    bool        mpi_io = false;

    void set_output_file_path();
    class_custom_cout ccout;


    //Mpi related constants
    int mpi_size, mpi_rank;
    hid_t plist_facc;
    hid_t plist_xfer;
    hid_t plist_lncr;


public:
    hid_t       file;

    explicit class_hdf5_file(const fs::path output_filename_, const fs::path output_dirname_ , bool create_dir_ = true, bool overwrite_ = false, bool mpi_io_=false);
    explicit class_hdf5_file(bool mpi_io_=false);

    ~class_hdf5_file(){
        H5Pclose(plist_facc);
        H5Pclose(plist_xfer);
        H5Pclose(plist_lncr);
        H5Fclose(file);
        std::cout << "Data written to file: " << output_file_path << std::endl;
    }
    template<typename DataType>
    void extend_dataset(const DataType &data, const std::string & dataset_relative_name);

    void create_group_link(const std::string &group_relative_name);

    template <typename AttrType>
    void create_group_link(const std::string &group_relative_name, const std::string &attribute_name, const AttrType attr);


    template <typename DataType>
    void write_dataset(const DataType &data, const std::string &dataset_relative_name);


    template <typename DataType>
    void read_dataset(DataType &data, const std::string &dataset_relative_name);

    template <typename AttrType>
    void write_attribute_to_dataset(const std::string &dataset_relative_name, const AttrType attribute,
                                    const std::string &attribute_name);

    template <typename AttrType>
    void write_attribute_to_group(const std::string &group_relative_name, const AttrType attribute,
                                  const std::string &attribute_name);



private:
    void initialize();
    void initialize_mpi();

    void extend_dataset(const std::string & dataset_relative_name, const int dim, const int extent );

    template <typename DataType>
    void create_dataset_link(const std::string &dataset_relative_name, hsize_t *chunk_dims);

    template <typename DataType>
    void create_dataset_link(const DataType &data, const std::string &dataset_relative_name);


    void select_hyperslab(hid_t &filespace, hid_t &memspace);


// MPI FUNCTIONS -- COMMENT OUT
//    template<typename DataType>
//    void extend_dataset_mpi(const DataType &data, const std::string & dataset_relative_name);

//    template <typename DataType>
//    void write_dataset_mpi(const DataType &data, const std::string &dataset_relative_name, bool each_mpi_thread = true);

//    template <typename DataType>
//    void read_dataset_mpi(DataType &data, const std::string &dataset_relative_name, bool each_mpi_thread = true);


//    template <typename AttrType>
//    void write_attribute_to_dataset_mpi(const std::string &dataset_relative_name, const AttrType attribute,
//    const std::string &attribute_name, bool each_mpi_thread = true);

//    void extend_dataset_mpi(const std::string & dataset_relative_name, const int dim, const int extent);
//    template <typename DataType>
//    void create_dataset_link_mpi(const DataType &data, const std::string &dataset_relative_name);



    template<typename DataType>
    hid_t get_DataType(){
        if constexpr (std::is_same<DataType, int>::value)          {return  H5Tcopy(H5T_NATIVE_INT);}
        if constexpr (std::is_same<DataType, long>::value)         {return  H5Tcopy(H5T_NATIVE_LONG);}
        if constexpr (std::is_same<DataType, double>::value)       {return  H5Tcopy(H5T_NATIVE_DOUBLE);}
        if constexpr (std::is_same<DataType, float>::value)        {return  H5Tcopy(H5T_NATIVE_FLOAT);}
        if constexpr (std::is_same<DataType, char>::value)         {return  H5Tcopy(H5T_C_S1);}
        if constexpr (std::is_same<DataType, std::string>::value)  {return  H5Tcopy(H5T_C_S1);}
        if constexpr (tc::has_member_scalar <DataType>::value)   {return  get_DataType<typename DataType::Scalar>();}
        if constexpr (tc::has_member_value_type<DataType>::value){return  get_DataType<typename DataType::value_type>();}
        std::cerr << "get_DataType could not match the type provided" << std::endl;
        exit(1);
    }



    template<typename DataType>
    int get_Size(const DataType &data){

        if constexpr (tc::is_tensor<DataType>::value){
            return (int) data.size();
        }

        if constexpr (std::is_arithmetic<DataType>::value){
            return 1;
        }
        if constexpr (tc::is_matrix <DataType>::value) {
            return data.size();
        }
        if constexpr(tc::is_instance_of<std::vector,DataType>::value){
            return data.size();
        }

        if constexpr(std::is_same<std::string, DataType>::value){
            return data.size();
        }
        std::cerr << "get_Rank can't match the type provided: " << typeid(data).name() << '\n';
        exit(1);
    }

    template<typename DataType>
    constexpr int get_Rank() const{

//        if constexpr (tc::is_tensor<DataType>::value){
//            return (int) data.rank();
//        }

        if constexpr (std::is_arithmetic<DataType>::value){
            return 1;
        }
        if constexpr (tc::is_matrix <DataType>::value) {
            return 2;
        }
        if constexpr(tc::is_instance_of<std::vector,DataType>::value){
            return 1;
        }

        if constexpr(std::is_same<std::string, DataType>::value){
            return 1;
        }
        std::cerr << "get_Rank can't match the type provided: " << typeid(DataType).name() << '\n';
        exit(1);
    }


    template<typename DataType>
    hid_t get_DataSpace_unlimited() const {

//        if constexpr (tc::is_tensor<DataType>::value){
//            int rank = data.rank();
//            hsize_t dims[rank];
////            hsize_t max_dims[rank] = {H5S_UNLIMITED};
//            std::copy(data.dimensions().begin(), data.dimensions().end(), dims);
//            return H5Screate_simple(rank, dims, nullptr);
//        }
        if constexpr (std::is_arithmetic<DataType>::value){
            const int rank = 1;
            hsize_t dims[rank] = {0};
            hsize_t max_dims[rank] = {H5S_UNLIMITED};
            return H5Screate_simple(rank, dims, max_dims);
        }
        if constexpr (tc::is_matrix <DataType>::value) {
            const int rank = 2;
            hsize_t dims[rank]      = {0, 0};
            hsize_t max_dims[rank]  = { H5S_UNLIMITED, H5S_UNLIMITED };
            return H5Screate_simple(rank, dims, max_dims);
        }
        if constexpr(tc::is_instance_of<std::vector,DataType>::value){
            const int rank = 1;
            hsize_t dims[rank] = {0};
            hsize_t max_dims[rank] = {H5S_UNLIMITED};
            return H5Screate_simple(rank, dims, max_dims);
        }

        if constexpr(std::is_same<std::string, DataType>::value){
            const int rank = 1;
            hsize_t dims[rank] = {0};
            hsize_t max_dims[rank] = {H5S_UNLIMITED};
            return H5Screate_simple(rank,dims, max_dims);
        }
        std::cerr << "get_DataSpace_unlimited can't match the type provided: " << typeid(DataType).name() << '\n';
        exit(1);
    }


    template<typename DataType>
    hid_t get_MemSpace(const DataType &data){

//        if constexpr (tc::is_tensor<DataType>::value){
//            int rank = get_Rank<DataType>();
//            hsize_t dims[rank];
//            std::copy(data.dimensions().begin(), data.dimensions().end(), dims);
//            return H5Screate_simple(rank, dims, nullptr);
//        }
        if constexpr (std::is_arithmetic<DataType>::value){
            constexpr int rank = 1;
            hsize_t dims[rank] = {1};
            return H5Screate_simple(rank, dims, nullptr);
        }
        if constexpr (tc::is_matrix <DataType>::value) {
            constexpr int rank = 2;
            hsize_t dims[rank] = { (hsize_t) data.rows(), (hsize_t) data.cols() };
            return H5Screate_simple(rank, dims, nullptr);
        }
        if constexpr(tc::is_instance_of<std::vector,DataType>::value){
            constexpr int rank = 1;
            hsize_t dims[rank] = {data.size()};
            return H5Screate_simple(rank, dims, NULL);
        }

        if constexpr(std::is_same<std::string, DataType>::value){
            constexpr int rank = 1;
            hsize_t dims[rank] = {1};
            return H5Screate_simple(rank,dims, nullptr);
        }
        std::cerr << "get_MemSpace can't match the type provided: " << typeid(DataType).name() << '\n';
        exit(1);
    }


    template <typename DataType>
    void set_ChunkDims(const DataType &data, hsize_t* dims){
//        if constexpr (tc::is_tensor<DataType>::value){
//            std::copy(data.dimensions().begin(), data.dimensions().end(), dims);
//            return;
//        }
        if constexpr (std::is_arithmetic<DataType>::value){
            dims[0]={100};
            return;
        }
        if constexpr (tc::is_matrix <DataType>::value) {
            dims[0] = (hsize_t) data.rows();
            dims[1] = (hsize_t) data.cols();
            return;
        }
        if constexpr(tc::is_instance_of<std::vector,DataType>::value){
            dims[0]={data.size()};
            return;
        }

        if constexpr(std::is_same<std::string, DataType>::value){
            dims[0]={data.size()};
            return;
        }
        std::cerr << "get_ChunkDims can't match the type provided: " << typeid(DataType).name() << '\n';
        exit(1);
    }


};




template <typename AttrType>
void class_hdf5_file::create_group_link(const std::string &group_relative_name, const std::string &attribute_name,
                                        const AttrType attr){
    if (!H5Lexists(file, group_relative_name.c_str(), H5P_DEFAULT)) {
        hid_t group          = H5Gcreate(file,group_relative_name.c_str(), plist_lncr,H5P_DEFAULT,H5P_DEFAULT);
        hid_t datatype       = get_DataType<AttrType>();
        hid_t dataspace      = get_MemSpace(attr);
        hid_t attribute      = H5Acreate(group, attribute_name.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT );
        retval               = H5Awrite(attribute, datatype, &attr);
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Aclose(attribute);
        H5Gclose(group);
    }
}



template <typename DataType>
void class_hdf5_file::create_dataset_link(const std::string &dataset_relative_name, hsize_t *chunk_dims){
    if (!H5Lexists(file, dataset_relative_name.c_str(), H5P_DEFAULT)){
        hid_t dataspace = get_DataSpace_unlimited<DataType>();
        hid_t datatype  = get_DataType<DataType>();
        hid_t dset_cpl  = H5Pcreate(H5P_DATASET_CREATE);
        auto ndims      = H5Sget_simple_extent_ndims(dataspace);
        H5Pset_layout(dset_cpl, H5D_CHUNKED);
        H5Pset_chunk(dset_cpl, ndims, chunk_dims);
        hid_t dataset = H5Dcreate(file,
                                  dataset_relative_name.c_str(),
                                  datatype,
                                  dataspace,
                                  plist_lncr,
                                  dset_cpl,
                                  H5P_DEFAULT);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Pclose(dset_cpl);
    }
}





template <typename DataType>
void class_hdf5_file::create_dataset_link(const DataType &data, const std::string &dataset_relative_name){
    if (!H5Lexists(file, dataset_relative_name.c_str(), H5P_DEFAULT)){
        hid_t dataspace = get_DataSpace_unlimited<DataType>();
        hid_t datatype  = get_DataType<DataType>();
        hid_t dset_cpl  = H5Pcreate(H5P_DATASET_CREATE);
        auto ndims      = H5Sget_simple_extent_ndims(dataspace);
        hsize_t chunk_dims[ndims];
        set_ChunkDims(data,chunk_dims);
        H5Pset_layout(dset_cpl, H5D_CHUNKED);
        H5Pset_chunk(dset_cpl, ndims, chunk_dims);
        hid_t dataset = H5Dcreate(file,
                                  dataset_relative_name.c_str(),
                                  datatype,
                                  dataspace,
                                  plist_lncr,
                                  dset_cpl,
                                  H5P_DEFAULT);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Pclose(dset_cpl);
    }
}






template<typename DataType>
void class_hdf5_file::extend_dataset(const DataType &data, const std::string & dataset_relative_name){
    if constexpr (tc::is_matrix<DataType>::value){
        extend_dataset(dataset_relative_name, 0, data.rows());
        hid_t dataset   = H5Dopen(file, dataset_relative_name.c_str(), H5P_DEFAULT);
        hid_t filespace = H5Dget_space(dataset);
        int ndims = H5Sget_simple_extent_ndims(filespace);
        hsize_t dims[ndims];
        H5Sget_simple_extent_dims(filespace,dims,NULL);
        H5Dclose(dataset);
        H5Sclose(filespace);
        if (dims[1] < (hsize_t) data.cols()){
            extend_dataset(dataset_relative_name, 1, data.cols());
        }

    }
    else{
        extend_dataset(dataset_relative_name, 0, get_Size(data));
    }
}




template <typename DataType>
void class_hdf5_file::write_dataset(const DataType &data, const std::string &dataset_relative_name){
    const int ndims = get_Rank<DataType>();
    hsize_t chunk_dims[ndims];
    set_ChunkDims(data,chunk_dims);
    create_dataset_link<DataType>(dataset_relative_name, chunk_dims);
    extend_dataset(data,dataset_relative_name);
    hid_t datatype  = get_DataType<DataType>();
    hid_t dataset   = H5Dopen(file,dataset_relative_name.c_str(), H5P_DEFAULT);
    hid_t filespace = H5Dget_space(dataset);
    hid_t memspace  = get_MemSpace(data);
    select_hyperslab(filespace,memspace);
    if constexpr(tc::has_member_data<DataType>::value){
        retval = H5Dwrite(dataset, datatype, memspace, filespace, H5P_DEFAULT, data.data());
    }
    if constexpr(std::is_arithmetic<DataType>::value){
        retval = H5Dwrite(dataset, datatype, memspace, filespace, H5P_DEFAULT, &data);
    }
    H5Tclose(datatype);
    H5Dclose(dataset);
}


template <typename DataType>
void class_hdf5_file::read_dataset(DataType &data, const std::string &dataset_relative_name){
    if (H5Lexists(file, dataset_relative_name.c_str(), H5P_DEFAULT)) {
        hid_t dataset = H5Dopen(file, dataset_relative_name.c_str(), H5P_DEFAULT);
        hid_t dataspace = H5Dget_space(dataset);
        hid_t datatype = H5Dget_type(dataset);
        int ndims = H5Sget_simple_extent_ndims(dataspace);
        hsize_t dims[ndims];
        H5Sget_simple_extent_dims(dataspace, dims, NULL);
        if constexpr(tc::is_matrix<DataType>::value) {
            data.resize(dims[0], dims[1]);
        }
        if constexpr(tc::is_instance_of<std::vector, DataType>::value) {
            assert(ndims == 1 and "Vector cannot take 2D datasets");
            data.resize(dims[0]);
        }
        H5LTread_dataset(file, dataset_relative_name.c_str(), datatype, data.data());
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Tclose(datatype);
    }else{
        std::cerr << "Attempted to read dataset that doesn't exist." << std::endl;
    }
}




template <typename AttrType>
void class_hdf5_file::write_attribute_to_dataset(const std::string &dataset_relative_name, const AttrType attribute,
                                                 const std::string &attribute_name){
    if (H5Lexists(file, dataset_relative_name.c_str(), H5P_DEFAULT)) {
        hid_t datatype = get_DataType<AttrType>();
        hid_t dataspace = get_MemSpace(attribute);
        hid_t dataset = H5Dopen(file, dataset_relative_name.c_str(), H5P_DEFAULT);
        hid_t attribute_id = H5Acreate(dataset, attribute_name.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        retval = H5Awrite(attribute_id, datatype, &attribute);
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Dclose(dataset);
        H5Aclose(attribute_id);

    }
}



template <typename AttrType>
void class_hdf5_file::write_attribute_to_group(const std::string &group_relative_name, const AttrType attribute,
                                               const std::string &attribute_name){
    if (H5Lexists(file, group_relative_name.c_str(), H5P_DEFAULT)) {

        hid_t datatype = get_DataType<AttrType>();
        hid_t dataspace = get_MemSpace(attribute);
        hid_t group = H5Gopen(file, group_relative_name.c_str(), H5P_DEFAULT);
        hid_t attribute_id = H5Acreate(group, attribute_name.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        retval = H5Awrite(attribute_id, datatype, &attribute);
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Dclose(group);
        H5Aclose(attribute_id);
    }
}






//
//
//template <typename DataType>
//void class_hdf5_file::create_dataset_link_mpi(const DataType &data, const std::string &dataset_relative_name){
//
//    for (int i = 0; i < mpi_size; i++) {
//        std::string dataset_relative_name_n = dataset_relative_name + std::to_string(i);
//        if (!H5Lexists(file, dataset_relative_name_n.c_str(), H5P_DEFAULT)) {
//            auto rank = get_Rank<DataType>();
//            hsize_t chunk_dims[rank];
//            set_ChunkDims(data,chunk_dims);
//            create_dataset_link<DataType>(dataset_relative_name_n, chunk_dims);
//        }
//    }
//}
//
//
//template<typename DataType>
//void class_hdf5_file::extend_dataset_mpi(const DataType &data, const std::string & dataset_relative_name){
//    if constexpr (tc::is_matrix<DataType>::value){
//        extend_dataset_mpi(dataset_relative_name, 0, data.rows());
//        std::string dataset_relative_name_n = dataset_relative_name + std::to_string(0);
//        hid_t dataset   = H5Dopen(file, dataset_relative_name_n.c_str(), H5P_DEFAULT);
//        hid_t filespace = H5Dget_space(dataset);
//        int ndims = H5Sget_simple_extent_ndims(filespace);
//        hsize_t dims[ndims];
//        H5Sget_simple_extent_dims(filespace,dims,NULL);
//        H5Dclose(dataset);
//        H5Sclose(filespace);
//        if (dims[1] < (hsize_t)data.cols()){
//            extend_dataset_mpi(dataset_relative_name, 1, data.cols());
//        }
//    }
//    else{
//        extend_dataset_mpi(dataset_relative_name, 0, get_Size(data));
//    }
//}
//
//
//
//
//template <typename DataType>
//void class_hdf5_file::write_dataset_mpi(const DataType &data, const std::string &dataset_relative_name, bool each_mpi_thread){
//    if (mpi_io){
//        if(!each_mpi_thread){
//            write_dataset(data,dataset_relative_name);
//            return;
//        }
//
//        create_dataset_link_mpi(data, dataset_relative_name);
//        extend_dataset_mpi(data,dataset_relative_name);
//        std::string dataset_relative_name_n = dataset_relative_name + std::to_string(mpi_rank);
//        hid_t dataset   = H5Dopen(file, dataset_relative_name_n.c_str(), H5P_DEFAULT);
//        hid_t datatype  = get_DataType<DataType>();
//        hid_t filespace = H5Dget_space(dataset);
//        hid_t memspace  = get_MemSpace(data);
//        select_hyperslab(filespace,memspace);
//
//        if constexpr(std::is_arithmetic<DataType>::value){
//            retval = H5Dwrite(dataset, datatype, memspace, filespace, plist_xfer, &data);
//        }
//
//        if constexpr(tc::has_member_data<DataType>::value){
//            retval = H5Dwrite(dataset, datatype, memspace, filespace, plist_xfer, data.template data());
//        }
//        if (retval < 0){
//            std::cerr << "Could not write to file!!!" <<std::endl;
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//
//        H5Tclose(datatype);
//        H5Dclose(dataset);
//        H5Sclose(filespace);
//        H5Sclose(memspace);
//
//    }else{
//        std::cerr << "mpi_io is set to false, please initialize hdf5 object with mpi_io=true" << std::endl;
//        exit(1);
//    }
//}
//
//
//
//template <typename DataType>
//void class_hdf5_file::read_dataset_mpi(DataType &data, const std::string &dataset_relative_name, bool each_mpi_thread){
//    if(mpi_io) {
//        if(!each_mpi_thread){
//            read_dataset(data,dataset_relative_name);
//            return;
//        }
//        for (int i = 0; i < mpi_size; i++){
//            if(mpi_rank == i) {
//                std::string dataset_relative_name_n = dataset_relative_name + std::to_string(i);
//                read_dataset(data, dataset_relative_name_n);
//            }
//        }
//
//    }else{
//        std::cerr << "mpi_io is set to false, please initialize hdf5 object with mpi_io=true" << std::endl;
//        exit(1);
//    }
//}
//
//
//template <typename AttrType>
//void class_hdf5_file::write_attribute_to_dataset_mpi(const std::string &dataset_relative_name, const AttrType attribute,
//                                                     const std::string &attribute_name, bool each_mpi_thread){
//    if(mpi_io) {
//        if(!each_mpi_thread){
//            write_attribute_to_dataset(dataset_relative_name, attribute,attribute_name);
//            return;
//        }
//        for (int i = 0; i < mpi_size; i++){
//            std::string dataset_relative_name_n = dataset_relative_name + std::to_string(i);
//            write_attribute_to_dataset(dataset_relative_name_n, attribute,attribute_name);
//
//        }
//    }else{
//        std::cerr << "mpi_io is set to false, please initialize hdf5 object with mpi_io=true" << std::endl;
//        exit(1);
//    }
//}

#endif
