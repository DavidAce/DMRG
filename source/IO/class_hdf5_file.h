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
#include <general/nmspc_type_check.h>


namespace fs = std::experimental::filesystem;
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

    void set_output_file_path();


    //Mpi related constants
    hid_t plist_facc;
    hid_t plist_xfer;
    hid_t plist_lncr;


    struct H5T_COMPLEX_STRUCT {
        double real;   /*real part*/
        double imag;   /*imaginary part*/
    };
    hid_t H5T_COMPLEX_DOUBLE;

    class DatasetProperties {
    public:
        hid_t                   datatype;
        hsize_t                 size;
        int                     ndims;
        std::vector<hsize_t>    dims;
        std::string             dset_name;

        ~DatasetProperties(){
            H5Tclose(datatype);
        }
    };

    class AttributeProperties {
    public:
        hid_t     datatype;
        hsize_t   size;
        int       ndims;
        std::vector<hsize_t> dims;
        std::string attr_name;
        std::string link_name;
        ~AttributeProperties(){
            H5Tclose(datatype);
        }
    };


public:
    hid_t       file;

    explicit class_hdf5_file(const std::string output_filename_, const std::string output_dirname_ , bool create_dir_ = true, bool overwrite_ = false);

    ~class_hdf5_file(){
        H5Pclose(plist_facc);
        H5Pclose(plist_xfer);
        H5Pclose(plist_lncr);
        H5Tclose(H5T_COMPLEX_DOUBLE);
        H5Fclose(file);
        std::cout << "Data written to file: " << output_file_path << std::endl;
    }

    template <typename DataType>
    void write_dataset(const DataType &data, const std::string &dataset_relative_name);

    template <typename DataType>
    void read_dataset(DataType &data, const std::string &dataset_relative_name);


    void create_group_link(const std::string &group_relative_name);


    template <typename AttrType>
    void write_attribute_to_file(const AttrType &attribute, const std::string attribute_name);

    template <typename AttrType>
    void write_attribute_to_dataset(const std::string &dataset_relative_name, const AttrType attribute,
                                    const std::string &attribute_name);


    template <typename AttrType>
    void write_attribute_to_group(const std::string &group_relative_name, const AttrType attribute,
                                  const std::string &attribute_name);


    template<typename DataType>
    std::vector<H5T_COMPLEX_STRUCT> convert_complex_data(const DataType &data);


private:
    void initialize();

    void extend_dataset(const std::string & dataset_relative_name, const int dim, const int extent );

    template<typename DataType>
    void extend_dataset(const DataType &data, const std::string & dataset_relative_name);

    void set_extent_dataset(const DatasetProperties &props);

    bool check_link_exists_recursively(std::string path);

    void create_dataset_link(const DatasetProperties &props);

    void select_hyperslab(const hid_t &filespace, const hid_t &memspace);

    template <typename DataType>
    void write_dataset(const DataType &data, const DatasetProperties &props);

    template <typename AttrType>
    void write_attribute_to_dataset(const AttrType &attribute, const AttributeProperties &aprops);

    template <typename AttrType>
    void write_attribute_to_group(const AttrType attribute, const AttributeProperties &aprops);

    template<typename DataType>
    constexpr hid_t get_DataType() const {
        if constexpr (std::is_same<DataType, int>::value)                 {return  H5Tcopy(H5T_NATIVE_INT);}
        if constexpr (std::is_same<DataType, long>::value)                {return  H5Tcopy(H5T_NATIVE_LONG);}
        if constexpr (std::is_same<DataType, unsigned int>::value)        {return  H5Tcopy(H5T_NATIVE_UINT);}
        if constexpr (std::is_same<DataType, unsigned long>::value)       {return  H5Tcopy(H5T_NATIVE_ULONG);}
        if constexpr (std::is_same<DataType, double>::value)              {return  H5Tcopy(H5T_NATIVE_DOUBLE);}
        if constexpr (std::is_same<DataType, float>::value)               {return  H5Tcopy(H5T_NATIVE_FLOAT);}
        if constexpr (std::is_same<DataType, std::complex<double>>::value){return  H5Tcopy(H5T_COMPLEX_DOUBLE);}
        if constexpr (std::is_same<DataType, H5T_COMPLEX_STRUCT>::value)  {return  H5Tcopy(H5T_COMPLEX_DOUBLE);}
        if constexpr (std::is_same<DataType, char>::value)                {return  H5Tcopy(H5T_C_S1);}
        if constexpr (std::is_same<DataType, std::string>::value)         {return  H5Tcopy(H5T_C_S1);}
        if constexpr (tc::is_eigen_tensor<DataType>::value)               {return  get_DataType<typename DataType::Scalar>();}
        if constexpr (tc::is_eigen_matrix<DataType>::value)               {return  get_DataType<typename DataType::Scalar>();}
        if constexpr (tc::is_eigen_array<DataType>::value)                {return  get_DataType<typename DataType::Scalar>();}
        if constexpr (tc::is_vector<DataType>::value)                     {return  get_DataType<typename DataType::value_type>();}
        if constexpr (tc::has_member_scalar <DataType>::value)            {return  get_DataType<typename DataType::Scalar>();}
        if constexpr (tc::has_member_value_type <DataType>::value)        {return  get_DataType<typename DataType::value_type>();}

        std::cerr << "get_DataType could not match the type provided" << std::endl;
        exit(1);
    }



    template<typename DataType>
    hsize_t get_Size(const DataType &data)const{
        if constexpr (tc::has_member_size<DataType>::value){return data.size();}
        if constexpr (std::is_arithmetic<DataType>::value) { return 1;}
        std::cerr << "get_Size can't match the type provided: " << typeid(data).name() << '\n';
        exit(1);
    }
    template<typename DataType>
    constexpr int get_Rank() const {
        std::cout << "Type is: " << tc::type_name<DataType>() << std::endl;
        if      constexpr (tc::is_eigen_tensor<DataType>::value){return (int) DataType::NumIndices;}
        else if constexpr (tc::is_eigen_matrix<DataType>::value){return 2; }
        else if constexpr (tc::is_eigen_array<DataType>::value){return 2; }
        else if constexpr (std::is_arithmetic<DataType>::value){return 1; }
        else if constexpr(tc::is_vector<DataType>::value){return 1;}
        else if constexpr(std::is_same<std::string, DataType>::value){return 1;}
        else {
            std::cerr << "get_Rank can't match the type provided: " << tc::type_name<DataType>() << '\n';
            tc::print_type_and_exit_compile_time<DataType>();
            exit(1);
        }
    }


    hid_t get_DataSpace_unlimited(int rank) const{
        std::vector<hsize_t> dims(rank);
        std::vector<hsize_t> max_dims(rank);
        std::fill_n(dims.begin(), rank, 0);
        std::fill_n(max_dims.begin(), rank, H5S_UNLIMITED);
        return H5Screate_simple(rank, dims.data(), max_dims.data());
    }

    template<typename DataType>
    hid_t get_MemSpace(const DataType &data) const {
        auto rank = get_Rank<DataType>();
        auto dims = get_Dimensions<DataType>(data);
        return H5Screate_simple(rank, dims.data(), nullptr);
    }


    template <typename DataType>
    std::vector<hsize_t> get_Dimensions(const DataType &data) const {
        int rank = get_Rank<DataType>();
        std::vector<hsize_t> dims(rank);
        if constexpr (tc::is_eigen_tensor<DataType>::value){
            std::copy(data.dimensions().begin(), data.dimensions().end(), dims.begin());
            return dims;
        }
        else if constexpr (std::is_arithmetic<DataType>::value){
            dims[0]={1};
            return dims;
        }
        else if constexpr (tc::is_eigen_matrix <DataType>::value) {
            dims[0] = (hsize_t) data.rows();
            dims[1] = (hsize_t) data.cols();
            return dims;
        }
        else if constexpr(tc::is_vector<DataType>::value){
            dims[0]={data.size()};
            return dims;
        }
        else if constexpr(std::is_same<std::string, DataType>::value){
            dims[0]={data.size()};
            return dims;
        }
        else{
            tc::print_type_and_exit_compile_time<DataType>();
            std::cerr << "get_Dimensions can't match the type provided: " << typeid(DataType).name() << '\n';
            exit(1);
        }

    }


};





template<typename DataType>
void class_hdf5_file::extend_dataset(const DataType &data, const std::string & dataset_relative_name){
    if constexpr (tc::is_eigen_matrix_or_array<DataType>()){
        extend_dataset(dataset_relative_name, 0, data.rows());
        hid_t dataset   = H5Dopen(file, dataset_relative_name.c_str(), H5P_DEFAULT);
        hid_t filespace = H5Dget_space(dataset);
        int ndims = H5Sget_simple_extent_ndims(filespace);
        std::vector<hsize_t> dims(ndims);
        H5Sget_simple_extent_dims(filespace,dims.data(),NULL);
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


template<typename DataType>
std::vector<class_hdf5_file::H5T_COMPLEX_STRUCT> class_hdf5_file::convert_complex_data(const DataType &data){

    if constexpr(tc::is_eigen_tensor<DataType>::value or tc::is_eigen_matrix_or_array<DataType>()){
        if constexpr(std::is_same<typename DataType::Scalar, std::complex<double>>::value) {
            std::vector<H5T_COMPLEX_STRUCT> new_data(data.size());
            for (int i = 0; i < data.size(); i++) {
                new_data[i].real = data(i).real();
                new_data[i].imag = data(i).imag();
            }
            return new_data;
        }
    }
    if constexpr(tc::is_vector<DataType>::value) {
        if constexpr(std::is_same<typename DataType::value_type, std::complex<double>>::value) {
            std::vector<H5T_COMPLEX_STRUCT> new_data(data.size());
            for (int i = 0; i < data.size(); i++) {
                new_data[i].real = data[i].real();
                new_data[i].imag = data[i].imag();
            }
            return new_data;
        }
    }

    if constexpr(std::is_arithmetic<DataType>::value and std::is_same<DataType, std::complex<double>>::value){
        std::vector<H5T_COMPLEX_STRUCT> new_data(1);
        new_data[0].real = data.real();
        new_data[0].imag = data.imag();
        return new_data;
    }
    //This should never happen, but is here so that we compile successfully.
    assert(NAN == NAN and "Big error! Tried to convert non-complex data_struct to complex");
    return     std::vector<H5T_COMPLEX_STRUCT>();
}


template <typename DataType>
void class_hdf5_file::write_dataset(const DataType &data, const DatasetProperties &props){
    create_dataset_link(props);
    set_extent_dataset(props);
    hid_t dataset   = H5Dopen(file,props.dset_name.c_str(), H5P_DEFAULT);
    hid_t filespace = H5Dget_space(dataset);
    hid_t memspace  = H5Screate_simple(props.ndims, props.dims.data(), nullptr);
    select_hyperslab(filespace,memspace);
    if constexpr(tc::has_member_data<DataType>::value){
        retval = H5Dwrite(dataset, props.datatype, memspace, filespace, H5P_DEFAULT, data.data());
    }
    if constexpr(std::is_arithmetic<DataType>::value){
        retval = H5Dwrite(dataset, props.datatype, memspace, filespace, H5P_DEFAULT, &data);
    }
    H5Dclose(dataset);
}

template <typename DataType>
void class_hdf5_file::write_dataset(const DataType &data, const std::string &dataset_relative_name){
    DatasetProperties props;
    props.datatype   = get_DataType<DataType>();
    props.size       = get_Size<DataType>(data);
    props.dset_name  = dataset_relative_name;
    props.ndims      = get_Rank<DataType>();
    props.dims       = get_Dimensions<DataType>(data);


    if (H5Tequal(get_DataType<DataType>(), H5T_COMPLEX_DOUBLE) and not std::is_same<std::vector<H5T_COMPLEX_STRUCT>, DataType>::value) {
        //If complex, convert data_struct to complex struct H5T_COMPLEX_DOUBLE
        auto new_data = convert_complex_data(data);
        write_dataset(new_data, props);
    }else if(H5Tequal(props.datatype, H5T_C_S1)){
        //If string, pad datatype
        retval = H5Tset_size(props.datatype, props.size);
        retval = H5Tset_strpad(props.datatype,H5T_STR_NULLTERM);
        write_dataset(data, props);
    }
    else{
        write_dataset(data, props);
    }
}


template <typename DataType>
void class_hdf5_file::read_dataset(DataType &data, const std::string &dataset_relative_name){
    if (H5Lexists(file, dataset_relative_name.c_str(), H5P_DEFAULT)) {
        hid_t dataset = H5Dopen(file, dataset_relative_name.c_str(), H5P_DEFAULT);
        hid_t dataspace = H5Dget_space(dataset);
        hid_t datatype = H5Dget_type(dataset);
        int ndims = H5Sget_simple_extent_ndims(dataspace);
        std::vector<hsize_t> dims(ndims);
        H5Sget_simple_extent_dims(dataspace, dims.data(), NULL);
        if constexpr(tc::is_eigen_matrix_or_array<DataType>()) {
            data.resize(dims[0], dims[1]);
        }
        if constexpr(tc::is_vector<DataType>::value) {
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
void class_hdf5_file::write_attribute_to_file(const AttrType &attribute, const std::string attribute_name){
    hid_t datatype          = get_DataType<AttrType>();
    hid_t dataspace         = get_MemSpace(attribute);
    auto size               = get_Size(attribute);
    hid_t attribute_id      = H5Acreate(file, attribute_name.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT );
    if constexpr (tc::has_member_data<AttrType>::value){
        retval                  = H5Awrite(attribute_id, datatype, attribute.data());
    }
    if constexpr (tc::has_member_c_str<AttrType>::value){
        retval                  = H5Tset_size(datatype, size);
        retval                  = H5Tset_strpad(datatype,H5T_STR_NULLTERM);
        retval                  = H5Awrite(attribute_id, datatype, attribute.c_str());
    }
    else{
        retval                  = H5Awrite(attribute_id, datatype, &attribute);
    }
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Aclose(attribute_id);
}



template <typename AttrType>
void class_hdf5_file::write_attribute_to_dataset(const AttrType &attribute, const AttributeProperties &aprops){
    if (check_link_exists_recursively(aprops.link_name)) {
        hid_t dataspace = get_MemSpace(attribute);
        hid_t dataset   = H5Dopen(file, aprops.link_name.c_str(), H5P_DEFAULT);
        hid_t attribute_id = H5Acreate(dataset, aprops.attr_name.c_str(), aprops.datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        retval = H5Awrite(attribute_id, aprops.datatype, &attribute);
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Aclose(attribute_id);
    }else{
        std::cerr << "Link '" << aprops.link_name << "' does not exist, yet attribute is being written." << std::endl;
        exit(1);
    }
}
template <typename AttrType>
void class_hdf5_file::write_attribute_to_dataset(const std::string &dataset_relative_name, const AttrType attribute,
                                                 const std::string &attribute_name){
    AttributeProperties aprops;
    aprops.datatype  = get_DataType<AttrType>();
    aprops.size      = get_Size(attribute);
    aprops.attr_name = attribute_name;
    aprops.link_name = dataset_relative_name;
    aprops.ndims     = get_Rank<AttrType>();
    aprops.dims      = get_Dimensions(attribute);

    write_attribute_to_dataset(attribute,aprops);


}

template <typename AttrType>
void class_hdf5_file::write_attribute_to_group(const AttrType attribute,
                                               const AttributeProperties &aprops){
    if (check_link_exists_recursively(aprops.link_name)) {
        hid_t dataspace = get_MemSpace(attribute);
        hid_t group = H5Gopen(file, aprops.link_name.c_str(), H5P_DEFAULT);
        hid_t attribute_id = H5Acreate(group, aprops.link_name.c_str(), aprops.datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        retval = H5Awrite(attribute_id, aprops.datatype, &attribute);
        H5Sclose(dataspace);
        H5Dclose(group);
        H5Aclose(attribute_id);
    }
}
template <typename AttrType>
void class_hdf5_file::write_attribute_to_group(const std::string &group_relative_name,
                                               const AttrType attribute,
                                               const std::string &attribute_name){

    AttributeProperties aprops;
    aprops.datatype  = get_DataType<AttrType>();
    aprops.size      = get_Size(attribute);
    aprops.attr_name = attribute_name;
    aprops.link_name = group_relative_name;
    aprops.ndims     = get_Rank<AttrType>();
    aprops.dims      = get_Dimensions(attribute);
    write_attribute_to_dataset(attribute,aprops);
}


#endif
