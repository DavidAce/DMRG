//
// Created by david on 8/1/17.
//

#ifndef DMRG_CLASS_HDF5_H
#define DMRG_CLASS_HDF5_H

//#include <hdf5.h>
//#include <hdf5_hl.h>
#include <iostream>
#include <string>
#include <iomanip>
#include "general/n_tensor_extra.h"
#include "general/n_math.h"

#include <experimental/filesystem>
#include <experimental/type_traits>

#include <H5Cpp.h>
#include <type_traits>
#include <typeinfo>
#include "general/n_math.h"
#include <gitversion.h>

namespace fs = std::experimental::filesystem::v1;



namespace TypeCheck{
    template <typename T> using Data_t          = decltype(std::declval<T>().data());
    template <typename T> using Size_t          = decltype(std::declval<T>().size());
    template <typename T> using Dims_t          = decltype(std::declval<T>().dimensions());
    template <typename T> using Rank_t          = decltype(std::declval<T>().rank());
    template <typename T> using Setz_t          = decltype(std::declval<T>().setZero());
    template <typename T> using Matr_t          = decltype(std::declval<T>().matrix());
    template <typename T> using Scal_t          = typename T::Scalar;
    template <typename T> using Valt_t          = typename T::value_type;

    template <typename T> using has_member_data         = std::experimental::is_detected<Data_t, T>;
    template <typename T> using has_member_size         = std::experimental::is_detected<Size_t, T>;
    template <typename T> using has_member_dimensions   = std::experimental::is_detected<Dims_t , T>;
    template <typename T> using has_member_scalar       = std::experimental::is_detected<Scal_t , T>;
    template <typename T> using has_member_value_type   = std::experimental::is_detected<Valt_t , T>;
    template <typename T> using is_eigen                = std::experimental::is_detected<Setz_t , T>;
    template <typename T> using is_matrix               = std::experimental::is_detected<Matr_t , T>;
    template <typename T> using is_tensor               = std::experimental::is_detected<Rank_t , T>;


    //This does not work for "non-type" class template parameters.
    //In fact it doesn't seem to work very well at all...
    template < template <typename...> class Template, typename T >
    struct is_instance_of : std::false_type {};

    template < template <typename...> class Template, typename... Args >
    struct is_instance_of< Template, Template<Args...> > : std::true_type {};

}
namespace tc = TypeCheck;



/*!
 \brief Writes and reads data to a binary hdf5-file.

 # HDF5 Class


*/

class class_hdf5 {
private:
    hid_t       file;
    herr_t      retval;
    fs::path    file_name          = fs::path("output.h5");
    fs::path    executable_path    = fs::current_path();
    fs::path    output_path;
    fs::path    output_dirname     = fs::path("output");
    bool create_dir;
    void set_output_path();

public:

    explicit class_hdf5(const fs::path &file_name_, const fs::path &output_dirname_ , bool create_dir_ = true);

    ~class_hdf5(){H5Fclose(file);}


    void open_file();
    void create_group(const std::string &group_relative_name);

    template <typename AttrType>
    void create_group(const std::string &group_relative_name, const std::string &attribute_name, const AttrType &attr);

    template <typename DataType, typename AttrType>
    void write_to_file(const DataType &data, const std::string &dataset_relative_name, const AttrType &attribute, const std::string &attribute_name ) {
        write_to_file(data,dataset_relative_name);
        write_attribute_to_dataset(dataset_relative_name, attribute, attribute_name);
    }

    template <typename DataType>
    void write_to_file(const DataType &data, const std::string &dataset_relative_name);

    template <typename AttrType>
    void write_attribute_to_dataset(const std::string &dataset_relative_name,  const AttrType &attribute ,const std::string &attribute_name);



private:

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
        std::cout << "get_DataType can't match the type provided" << '\n';
        exit(1);
    }


    template<typename DataType>
    hid_t get_DataSpace(const DataType &data){

        if constexpr (tc::is_tensor<DataType>::value){
                int rank = data.rank();
                hsize_t dims[rank];
                std::copy(data.dimensions().begin(), data.dimensions().end(), dims);
                return H5Screate_simple(rank, dims, nullptr);
        }
        if constexpr (std::is_arithmetic<DataType>::value){
            const int rank = 1;
            hsize_t dims[rank] = {1};
            return H5Screate_simple(rank, dims, nullptr);
        }
        if constexpr (tc::is_matrix <DataType>::value) {
                const int rank = 2;
                hsize_t dims[rank] = { (hsize_t) data.rows(), (hsize_t) data.cols() };
                return H5Screate_simple(rank, dims, nullptr);
        }
        if constexpr(tc::is_instance_of<std::vector,DataType>::value){
            const int rankk = 1;
            hsize_t dimss[rankk] = {data.size()};
            return H5Screate_simple(rankk, dimss, nullptr);
        }

        if constexpr(std::is_same<std::string, DataType>::value){
            const int rankk = 1;
            hsize_t dimss[rankk] = {1};
            return H5Screate_simple(rankk, dimss, nullptr);
        }
        std::cout << "get_DataSpace can't match the type provided: " << typeid(data).name() << '\n';
        exit(1);
    }



};










#endif
