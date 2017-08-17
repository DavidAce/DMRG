//
// Created by david on 8/1/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_HDF5_H
#define FINITE_DMRG_EIGEN_CLASS_HDF5_H

//#include <hdf5.h>
//#include <hdf5_hl.h>
#include <iostream>
#include <n_tensor_extra.h>
#include <n_math.h>

#include <experimental/filesystem>
#include <experimental/type_traits>

//#include <H5CommonFG.h>
//#include <H5File.h>
#include <H5Cpp.h>
//#include <hdf5.h>
#include <type_traits>
#include <typeinfo>
#include <n_math.h>
namespace fs = std::experimental::filesystem::v1;
using namespace Textra;



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
    fs::path    file_name;
    fs::path    executable_path    = fs::current_path();
    fs::path    output_path;
    fs::path    output_dirname;
    bool create_dir;
    void set_output_path();

public:

    explicit class_hdf5(fs::path file_name_ = fs::path("dmrg_output.h5"), fs::path output_dirname_ = fs::path("output"), bool create_dir_ = true):
            file_name(file_name_),
            output_dirname(output_dirname_),
            create_dir(create_dir_)
    {
        set_output_path();
        file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC,  H5P_DEFAULT, H5P_DEFAULT);

        //Put git revision in file attribute
        std::string gitrevision = std::to_string(GIT_REVISION);
        hid_t datatype          = get_DataType<std::string>();
        hid_t dataspace         = get_DataSpace(gitrevision);
        retval                  = H5Tset_size(datatype, gitrevision.size());
        retval                  = H5Tset_strpad(datatype,H5T_STR_NULLTERM);
        hid_t attribute         = H5Acreate(file, "GIT REVISION", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT );
        retval                  = H5Awrite(attribute, datatype, gitrevision.c_str());
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Aclose(attribute);
    }

    ~class_hdf5(){
        H5Fclose(file);
    }


    void open_file(){
        file = H5Fopen (file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }

    void create_group(const std::string group_relative_name){
        hid_t lcpl            = H5Pcreate(H5P_LINK_CREATE);   //Create missing intermediate group if they don't exist
        retval                = H5Pset_create_intermediate_group(lcpl, 1);
        hid_t group           = H5Gcreate(file,group_relative_name.c_str(), lcpl,H5P_DEFAULT,H5P_DEFAULT);
        H5Pclose(lcpl);
        H5Gclose(group);
    }


    template <typename AttrType>
    void create_group(const std::string group_relative_name, const std::string attribute_name, const AttrType &attr){
        hid_t lcpl           = H5Pcreate(H5P_LINK_CREATE);   //Create missing intermediate group if they don't exist
        retval               = H5Pset_create_intermediate_group(lcpl, 1);
        hid_t group          = H5Gcreate(file,group_relative_name.c_str(), lcpl,H5P_DEFAULT,H5P_DEFAULT);
        hid_t datatype       = get_DataType<AttrType>();
        hid_t dataspace      = get_DataSpace(attr);
        hid_t attribute      = H5Acreate(group, attribute_name.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT );
        retval               = H5Awrite(attribute, datatype, &attr);
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Aclose(attribute);
        H5Pclose(lcpl);
        H5Gclose(group);
    }


    template <typename DataType>
    void write_to_file(const DataType &data, const std::string dataset_relative_name){
        hid_t dataspace = get_DataSpace(data);
        hid_t datatype  = get_DataType<DataType>();
        hid_t lcpl      = H5Pcreate(H5P_LINK_CREATE);   //Create missing intermediate group if they don't exist
        retval          = H5Pset_create_intermediate_group(lcpl, 1);
        hid_t dataset   = H5Dcreate(file,
                                    dataset_relative_name.c_str(),
                                    datatype,
                                    dataspace,
                                    lcpl,
                                    H5P_DEFAULT,
                                    H5P_DEFAULT);

        if constexpr(tc::has_member_data<DataType>::value){
            retval = H5Dwrite(dataset,
                              datatype,
                              H5S_ALL,
                              H5S_ALL,
                              H5P_DEFAULT,
                              data.data());
        }
        if constexpr(std::is_arithmetic<DataType>::value){
                retval = H5Dwrite(dataset,
                                  datatype,
                                  H5S_ALL,
                                  H5S_ALL,
                                  H5P_DEFAULT,
                                  &data);
        }
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Dclose(dataset);
        H5Pclose(lcpl);
    }

    template <typename DataType>
    void write_attribute_to_dataset(const std::string dataset_relative_name,  const DataType &data ,const std::string attribute_name){
        hid_t datatype       = get_DataType<DataType>();
        hid_t dataspace      = get_DataSpace(data);
        hid_t dataset        = H5Dopen(file, dataset_relative_name.c_str(),H5P_DEFAULT);
        hid_t attribute      = H5Acreate(dataset, attribute_name.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT );
        retval               = H5Awrite(attribute, datatype, &data);
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Dclose(dataset);
        H5Aclose(attribute);
    }

    template <typename DataType, typename AttrType>
    void write_to_file(const DataType &data, const std::string dataset_relative_name, const AttrType attribute, const std::string attribute_name ) {
        write_to_file(data,dataset_relative_name);
        write_attribute_to_dataset(dataset_relative_name, attribute, attribute_name);
    }

private:

    template<typename DataType>
    hid_t get_DataType(){
        if(std::is_same<DataType, int>::value)          {return  H5Tcopy(H5T_NATIVE_INT);}
        if(std::is_same<DataType, long>::value)         {return  H5Tcopy(H5T_NATIVE_LONG);}
        if(std::is_same<DataType, double>::value)       {return  H5Tcopy(H5T_NATIVE_DOUBLE);}
        if(std::is_same<DataType, float>::value)        {return  H5Tcopy(H5T_NATIVE_FLOAT);}
        if(std::is_same<DataType, char>::value)         {return  H5Tcopy(H5T_C_S1);}
        if(std::is_same<DataType, std::string>::value)  {return  H5Tcopy(H5T_C_S1);}
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
            int rank = 1;
            hsize_t dims[rank] = {1};
            return H5Screate_simple(rank, dims, nullptr);
        }
        if constexpr (tc::is_matrix <DataType>::value) {
                int rank = 2;
                hsize_t dims[rank] = { (hsize_t) data.rows(), (hsize_t) data.cols() };
                return H5Screate_simple(rank, dims, nullptr);
        }
        if constexpr(tc::is_instance_of<std::vector,DataType>::value){
            int rankk = 1;
            hsize_t dimss[rankk] = {data.size()};
            return H5Screate_simple(rankk, dimss, nullptr);
        }

        if constexpr(std::is_same<std::string, DataType>::value){
            int rankk = 1;
            hsize_t dimss[rankk] = {1};
            return H5Screate_simple(rankk, dimss, nullptr);
        }
        std::cout << "get_DataSpace can't match the type provided" << '\n';
        exit(1);
    }

};

#endif //FINITE_DMRG_EIGEN_CLASS_HDF5_H
