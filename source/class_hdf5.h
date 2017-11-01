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
#include <n_tensor_extra.h>
#include <n_math.h>

#include <experimental/filesystem>
#include <experimental/type_traits>

#include <H5Cpp.h>
#include <type_traits>
#include <typeinfo>
#include <n_math.h>
#include <gitversion.h>

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
    fs::path    file_name          = fs::path("output.h5");
    fs::path    executable_path    = fs::current_path();
    fs::path    output_path;
    fs::path    output_dirname     = fs::path("output");
    bool create_dir;

public:

    explicit class_hdf5(const fs::path &file_name_, const fs::path &output_dirname_ , bool create_dir_ = true):
            file_name(file_name_),
            output_dirname(output_dirname_),
            create_dir(create_dir_)
    {
        set_output_path();
        file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC,  H5P_DEFAULT, H5P_DEFAULT);

        //Put git revision in file attribute
        std::string gitversion = "Branch: " + GIT::BRANCH + " | Commit hash: " + GIT::COMMIT_HASH + " | Revision: " + GIT::REVISION;
        hid_t datatype          = get_DataType<std::string>();
        hid_t dataspace         = get_DataSpace(gitversion);
        retval                  = H5Tset_size(datatype, gitversion.size());
        retval                  = H5Tset_strpad(datatype,H5T_STR_NULLTERM);
        hid_t attribute         = H5Acreate(file, "GIT REVISION", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT );
        retval                  = H5Awrite(attribute, datatype, gitversion.c_str());
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

    void create_group(const std::string &group_relative_name){
        hid_t lcpl            = H5Pcreate(H5P_LINK_CREATE);   //Create missing intermediate group if they don't exist
        retval                = H5Pset_create_intermediate_group(lcpl, 1);
        hid_t group           = H5Gcreate(file,group_relative_name.c_str(), lcpl,H5P_DEFAULT,H5P_DEFAULT);
        H5Pclose(lcpl);
        H5Gclose(group);
    }


    template <typename AttrType>
    void create_group(const std::string &group_relative_name, const std::string &attribute_name, const AttrType &attr){
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
    void write_to_file(const DataType &data, const std::string &dataset_relative_name){
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

    template <typename AttrType>
    void write_attribute_to_dataset(const std::string &dataset_relative_name,  const AttrType &attribute ,const std::string &attribute_name){
        hid_t datatype       = get_DataType<AttrType>();
        hid_t dataspace      = get_DataSpace(attribute);
        hid_t dataset        = H5Dopen(file, dataset_relative_name.c_str(),H5P_DEFAULT);
        hid_t attribute_id   = H5Acreate(dataset, attribute_name.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT );
        retval               = H5Awrite(attribute_id, datatype, &attribute);
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Dclose(dataset);
        H5Aclose(attribute_id);
    }

    template <typename DataType, typename AttrType>
    void write_to_file(const DataType &data, const std::string &dataset_relative_name, const AttrType &attribute, const std::string &attribute_name ) {
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

//        if constexpr(tc::is_specialization_of<std::vector,DataType>::value){
//                int rankk = 1;
//                hsize_t dimss[rankk] = {data.size()};
//                return H5Screate_simple(rankk, dimss, nullptr);
//        }


        if constexpr(std::is_same<std::string, DataType>::value){
            const int rankk = 1;
            hsize_t dimss[rankk] = {1};
            return H5Screate_simple(rankk, dimss, nullptr);
        }
//        if constexpr(true) {
//            std::cout << "get_DataSpace can't match the type provided: " <<  decltype(data)::dummy << '\n';
//        }
        std::cout << "get_DataSpace can't match the type provided: " << typeid(data).name() << '\n';
//        std::cout << "is_specialization_of: " << tc::is_specialization_of<std::vector,DataType>::value << '\n';
        exit(1);
    }

    void set_output_path() {
        output_path = fs::system_complete(output_dirname);
        if (create_dir) {
            //Create directory
            fs::create_directories(output_path);
            output_path = fs::canonical(output_dirname);
            file_name = output_path / file_name;
        } else {
            //Try to find the directory
            output_path = executable_path / output_dirname.stem();
            while (true) {
                std::cout << "Searching for directory: " << output_path << '\n';
                if (fs::exists(output_path)) {
                    output_path = fs::canonical(output_path);
                    std::cout << "Found " << output_dirname << " directory: " << output_path << '\n';
                    break;
                } else if (output_path.parent_path().has_parent_path()) {
                    output_path = output_path.parent_path().parent_path() / output_dirname.stem();
                } else {
                    std::cout << "ERROR: " << output_dirname << " folder does not exist and create_dir == false. \n";
                    std::cout << "Either create an " << output_dirname << "  directory manually or pass create_dir = true.\n";
                    exit(1);
                }
            }
        }
    }

};







template<typename DataType, typename AttrType = int>
class class_hdf5_dataset_buffer : public std::vector<DataType>{
private:
    class_hdf5 &hdf5;
    std::string group_name      = "default_group";
    std::string dataset_name    = "default_data";
    std::string attribute_name  = "default_attr";
    AttrType attribute          = 0;
    bool attribute_set          = false;
    int iteration               = 0;
    int pad_to_width            = 3;
    std::string dataset_relative_name;

    std::string left_pad(const std::string &temp){
        std::ostringstream ostr;
        ostr << std::setfill('0') << std::setw(pad_to_width) << temp;
        return ostr.str();
    }

public:
    explicit class_hdf5_dataset_buffer(class_hdf5 &hdf5_):hdf5(hdf5_){
    }

    class_hdf5_dataset_buffer(class_hdf5 &hdf5_,
                              const std::string &group_name_,
                              const int &iteration_,
                              const std::string &dataset_name_
                              )
            :hdf5(hdf5_)
    {
        hdf5_set_all(group_name_, iteration_,dataset_name_);
    }

    class_hdf5_dataset_buffer(class_hdf5 &hdf5_,
                              const std::string &group_name_,
                              const int &iteration_,
                              const std::string &dataset_name_,
                              const AttrType &attribute_,
                              const std::string &attribute_name_
    )
            :hdf5(hdf5_),
             group_name(group_name_),
             iteration(iteration_),
             dataset_name(dataset_name_),
             attribute(attribute_),
             attribute_name(attribute_name_),

             attribute_set(true)
    {

    }

    ~class_hdf5_dataset_buffer(){
        hdf5_write();
    }

    void hdf5_set_group_name    (const std::string &group_name_)   {
        group_name     = group_name_;
    }
    void hdf5_set_dataset_name  (const std::string &dataset_name_) {
        dataset_name   = dataset_name_;
    }
    void hdf5_set_iteration     (const int &iteration_)             {
        iteration      = iteration_;
    }

    void hdf5_set_width     (const int &pad_to_width_)             {
        pad_to_width      = pad_to_width_;
    }

    void hdf5_set_all(const std::string &group_name_, const int iteration_, const std::string &dataset_name_){
        group_name     = group_name_;
        iteration      = iteration_;
        dataset_name   = dataset_name_;
    }

    void hdf5_set_attribute(const AttrType &attribute_, const std::string &attribute_name_){
        attribute_name  = attribute_name_;
        attribute       = attribute_;
        attribute_set   = true;
    }

    void hdf5_refresh_relative_name(){
        dataset_relative_name = "/" + group_name + "_" + left_pad(std::to_string(iteration)) + "/" + dataset_name ;
    }


    void hdf5_write() {
        if (!this->empty()) {
            hdf5_refresh_relative_name();
            if (attribute_set) {
                hdf5.write_to_file(static_cast<std::vector<DataType>>(*this), dataset_relative_name, attribute,
                                   attribute_name);
            } else {
                hdf5.write_to_file(static_cast<std::vector<DataType>>(*this), dataset_relative_name);
            }
            this->clear();
        }
    }
};





#endif
