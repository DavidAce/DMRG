//
// Created by david on 2017-12-02.
//

#include "class_hdf5.h"

class_hdf5::class_hdf5(const fs::path &file_name_, const fs::path &output_dirname_ , bool create_dir_):
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


void class_hdf5::set_output_path() {
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


void class_hdf5::open_file(){
    file = H5Fopen (file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
}

void class_hdf5::create_group(const std::string &group_relative_name){
    hid_t lcpl            = H5Pcreate(H5P_LINK_CREATE);   //Create missing intermediate group if they don't exist
    retval                = H5Pset_create_intermediate_group(lcpl, 1);
    hid_t group           = H5Gcreate(file,group_relative_name.c_str(), lcpl,H5P_DEFAULT,H5P_DEFAULT);
    H5Pclose(lcpl);
    H5Gclose(group);
}


template <typename AttrType>
void class_hdf5::create_group(const std::string &group_relative_name, const std::string &attribute_name, const AttrType &attr){
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

//Explicit instantiations
template void class_hdf5::create_group(const std::string &, const std::string &, const int &);
//template void class_hdf5::create_group(const std::string &, const std::string &, const long &);
//template void class_hdf5::create_group(const std::string &, const std::string &, const double &);




template <typename DataType>
void class_hdf5::write_to_file(const DataType &data, const std::string &dataset_relative_name){
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

//Explicit instantiations
template void class_hdf5::write_to_file(const int &, const std::string &);
template void class_hdf5::write_to_file(const long &, const std::string &);
template void class_hdf5::write_to_file(const double &, const std::string &);
template void class_hdf5::write_to_file(const float &, const std::string &);
template void class_hdf5::write_to_file(const char &, const std::string &);
template void class_hdf5::write_to_file(const std::string &, const std::string &);
template void class_hdf5::write_to_file(const std::vector<int> &, const std::string &);
template void class_hdf5::write_to_file(const std::vector<long> &, const std::string &);
template void class_hdf5::write_to_file(const std::vector<double> &, const std::string &);
template void class_hdf5::write_to_file(const std::vector<std::string> &, const std::string &);



template <typename AttrType>
void class_hdf5::write_attribute_to_dataset(const std::string &dataset_relative_name,  const AttrType &attribute ,const std::string &attribute_name){
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

//Explicit instantiations
template void class_hdf5::write_attribute_to_dataset(const std::string &,  const int & ,const std::string &);
template void class_hdf5::write_attribute_to_dataset(const std::string &,  const long & ,const std::string &);
