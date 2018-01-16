//
// Created by david on 2017-12-02.
//

#include "class_hdf5.h"
#include <fstream>
#include <istream>
#include <directory.h>

class_hdf5::class_hdf5(const fs::path &output_filename_, const fs::path &output_folder_ , bool create_dir_, bool overwrite_){
    output_filename = output_filename_;
    output_folder   = output_folder_;
    create_dir      = create_dir_;
    overwrite       = overwrite_;
    initialize();
}

class_hdf5::class_hdf5(){
    output_filename = settings::hdf5::output_filename;
    output_folder   = settings::hdf5::output_folder;
    create_dir      = settings::hdf5::create_dir_if_not_found;
    overwrite       = settings::hdf5::overwrite_file_if_found;
    initialize();
}

void class_hdf5::initialize(){
    set_output_file_path();
    file = H5Fcreate(output_file_path.c_str(), H5F_ACC_TRUNC,  H5P_DEFAULT, H5P_DEFAULT);
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


void class_hdf5::set_output_file_path() {
    fs::path project_folder  = fs::canonical(directory::PROJECT_DIR);
    if(output_folder.has_filename() and output_folder.has_extension()){output_folder.remove_filename();}
    ccout(1) << "output_folder     = " << output_folder << std::endl;
    ccout(1) << "output_filename   = " << output_filename << std::endl;

    //Convert the output_folder to a full path. If relative, make it relative to project root folder!
    fs::path output_folder_abs;
    if (output_folder.is_relative()) {
        output_folder_abs = fs::system_complete(project_folder / output_folder);
    } else if (output_folder.is_absolute()) {
        output_folder_abs = fs::system_complete(output_folder);
    }else{
        std::cerr << "Error: Provided path for hdf5 output folder is neither relative nor absolute: [output_folder]: " << output_folder << std::endl;
        exit(1);
    }

    //Create directory if create_dir == true, always relative to project root folder
    if(create_dir) {
        if (fs::create_directories(output_folder_abs)){
            ccout(1) << "Created directory: " <<  output_folder_abs << std::endl;
        }else{
            ccout(1) << "Directory already exists: " <<  output_folder_abs << std::endl;
        }
    }

    //Now the output_folder_abs may or may not exist
    //Now the output_file may or may not exist

    //If the folder exists use it.
    output_file_path = fs::system_complete(output_folder_abs / output_filename);
    if (fs::exists(output_folder_abs)) {
        if(fs::exists(output_file_path)) {
            ccout(1) << "File exists already: " << output_file_path << std::endl;
            if (settings::hdf5::overwrite_file_if_found) {
                ccout(1) << "Overwrite mode is TRUE." << std::endl;
                return;
            }
            else{
                int i = 1;
                while (fs::exists(output_file_path)){
                    output_file_path.replace_filename(output_filename.stem().string() + "-" + std::to_string(i++) + output_filename.extension().string() );
                }
                ccout(1) << "Renamed output file: " << output_filename << "  -->  " << output_file_path.filename() << std::endl;
                return;
            }
        }
        ccout(1) << "Creating new file: " << output_file_path << std::endl;
        return;
    }
    else{
        ccout(1) << "Output folder does not exist in the given path: " << output_folder_abs << std::endl;
    }

    //Now the output_folder_abs definitely does not exist in the given path
    //As a last resort, try finding the output folder somewhere inside the project root folder, excluding .git/ and libs/ and docs/
    ccout(1) << "Searching recursively for folder " <<  output_folder_abs.stem() << " in project root directory..." << std::endl;
    for(auto& p: fs::recursive_directory_iterator(project_folder)) {
        if (p.path().has_filename() and p.path().has_extension()){continue;}
        fs::path trimmed_path =  p.path().string().substr(p.path().string().find(project_folder) + project_folder.string().size());
        if (trimmed_path.string().find(std::string(".git")) != std::string::npos){continue;}
        if (trimmed_path.string().find(std::string("libs")) != std::string::npos){continue;}
        if (trimmed_path.string().find(std::string("docs")) != std::string::npos){continue;}
        if (trimmed_path.string().find(std::string("cmake")) != std::string::npos){continue;}
        if (trimmed_path.string().find(std::string("build")) != std::string::npos){continue;}
        if (p.path().stem() == output_folder_abs.stem())  {
            ccout(1) << "Found output folder: " << p.path() << std::endl;
            output_folder_abs = p.path();
            output_file_path = fs::system_complete(output_folder_abs / output_filename);
            break;
        }
    }


    if(fs::exists(output_file_path)) {
        ccout(1) << "File exists already: " << output_file_path << std::endl;
        if (overwrite) {
            ccout(1) << "Overwrite mode is TRUE." << std::endl;
            return;
        } else {
            int i = 1;
            while (fs::exists(output_file_path)){
                output_file_path.replace_filename(output_filename.stem().string() + "-" + std::to_string(i++) + output_filename.extension().string() );
            }
            ccout(1) << "Renamed file: " << output_filename << "  -->  " << output_file_path.filename() << std::endl;
            return;
        }

    }
}


void class_hdf5::open_file(){
    file = H5Fopen (output_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
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
