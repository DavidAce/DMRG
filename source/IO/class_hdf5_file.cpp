//
// Created by david on 2018-02-07.
//
//#include <directory.h>
#include <gitversion.h>
#include <IO/class_hdf5_file.h>

class_hdf5_file::class_hdf5_file(const std::string output_filename_, const std::string output_folder_, bool create_dir_, bool overwrite_){
    output_filename = output_filename_;
    output_folder   = output_folder_;
    create_dir      = create_dir_;
    overwrite       = overwrite_;
  /*
  * Check if gzip compression is available and can be used for both
  * compression and decompression.  Normally we do not perform error
  * checking in these examples for the sake of clarity, but in this
  * case we will make an exception because this filter is an
  * optional part of the hdf5 library.
  */
//    herr_t          status;
//    htri_t          avail;
//    H5Z_filter_t    filter_type;
//    unsigned int    flags, filter_info;
//    avail = H5Zfilter_avail(H5Z_FILTER_DEFLATE);
//    if (!avail) {
//        std::cerr << "gzip filter not available.\n" << std::endl;
//        exit(1);
//    }
//    status = H5Zget_filter_info (H5Z_FILTER_DEFLATE, &filter_info);
//    if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ||
//         !(filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED) ) {
//        std::cerr <<  "gzip filter not available for encoding and decoding." << std::endl;
//        exit(1);
//    }

    initialize();

    H5T_COMPLEX_DOUBLE = H5Tcreate (H5T_COMPOUND, sizeof(H5T_COMPLEX_STRUCT));
    H5Tinsert (H5T_COMPLEX_DOUBLE, "real", HOFFSET(H5T_COMPLEX_STRUCT,real), H5T_NATIVE_DOUBLE);
    H5Tinsert (H5T_COMPLEX_DOUBLE, "imag", HOFFSET(H5T_COMPLEX_STRUCT,imag), H5T_NATIVE_DOUBLE);
}

void class_hdf5_file::initialize(){
    set_output_file_path();
    plist_facc = H5Pcreate(H5P_FILE_ACCESS);
    plist_lncr = H5Pcreate(H5P_LINK_CREATE);   //Create missing intermediate group if they don't exist
    plist_xfer = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_create_intermediate_group(plist_lncr, 1);

    file = H5Fcreate(output_file_path.c_str(), H5F_ACC_TRUNC,  H5P_DEFAULT, H5P_DEFAULT);
    //Put git revision in file attribute
    std::string gitversion = "Branch: " + GIT::BRANCH + " | Commit hash: " + GIT::COMMIT_HASH + " | Revision: " + GIT::REVISION;
    write_attribute_to_file(gitversion, "GIT REVISION");
}

void class_hdf5_file::set_output_file_path() {
//    if(output_folder.has_filename() and output_folder.has_extension()){output_folder.remove_filename();}
    fs::path current_folder    = fs::current_path();
    fs::path output_folder_abs = fs::system_complete(output_folder);
    std::cout << "output_folder     = " << output_folder << std::endl;
    std::cout << "output_folder_abs = " << output_folder_abs << std::endl;
    std::cout << "output_filename   = " << output_filename << std::endl;

      //Create directory if create_dir == true, always relative to executable
    if(create_dir) {
        if (fs::create_directories(output_folder_abs)){
            output_folder_abs = fs::canonical(output_folder_abs);
            std::cout << "Created directory: " <<  output_folder_abs << std::endl;
        }else{
            output_folder_abs = fs::canonical(output_folder_abs);
            std::cout << "Directory already exists: " <<  output_folder_abs << std::endl;
        }
    }

    //Now the output_folder_abs may or may not exist
    //Now the output_file may or may not exist

    //If the folder exists use it.
    output_file_path = fs::system_complete(output_folder_abs / output_filename);
    if (fs::exists(output_folder_abs)) {
        if(fs::exists(output_file_path)) {
            output_file_path = fs::canonical(output_file_path);
            std::cout << "File already exists: " << output_file_path << std::endl;
            if (overwrite) {
                std::cout << "Overwrite mode is TRUE." << std::endl;
                return;
            }
            else{
                int i = 1;
                while (fs::exists(output_file_path)){
                    output_file_path.replace_filename(output_filename.stem().string() + "-" + std::to_string(i++) + output_filename.extension().string() );
                }
                std::cout << "Renamed output file: " << output_filename << "  -->  " << output_file_path.filename() << std::endl;
                return;
            }
        }
        std::cout << "Creating new file: " << output_file_path << std::endl;
        return;
    }
    else{
        std::cout << "Output folder does not exist in the given path: " << output_folder_abs << std::endl;
    }

    //Now the output_folder_abs definitely does not exist in the given path
    //As a last resort, try finding the output folder somewhere inside the project root folder, excluding .git/ and libs/ and docs/
//    std::cout << "Searching recursively for folder " <<  output_folder_abs.stem() << " in project root directory..." << std::endl;
//    for(auto& p: fs::recursive_directory_iterator(current_folder)) {
//        if (p.path().has_filename() and p.path().has_extension()){continue;}
//        fs::path trimmed_path =  p.path().string().substr(p.path().string().find(current_folder) + current_folder.string().size());
//        if (trimmed_path.string().find(std::string(".git")) != std::string::npos){continue;}
//        if (trimmed_path.string().find(std::string("libs")) != std::string::npos){continue;}
//        if (trimmed_path.string().find(std::string("docs")) != std::string::npos){continue;}
//        if (trimmed_path.string().find(std::string("cmake")) != std::string::npos){continue;}
//        if (trimmed_path.string().find(std::string("build")) != std::string::npos){continue;}
//        if (trimmed_path.string().find(std::string("backup")) != std::string::npos){continue;}
//        if (trimmed_path.string().find(std::string("source")) != std::string::npos){continue;}
//        if (p.path().stem() == output_folder_abs.stem())  {
//            std::cout << "Found output folder: " << p.path() << std::endl;
//            output_folder_abs = p.path();
//            output_file_path = fs::system_complete(output_folder_abs / output_filename);
//            break;
//        }
//    }


    if(fs::exists(output_file_path)) {
        output_file_path = fs::canonical(output_file_path);
        std::cout << "File exists already: " << output_file_path << std::endl;
        if (overwrite) {
            std::cout << "Overwrite mode is TRUE." << std::endl;
            return;
        } else {
            int i = 1;
            while (fs::exists(output_file_path)){
                output_file_path.replace_filename(output_filename.stem().string() + "-" + std::to_string(i++) + output_filename.extension().string() );
            }
            std::cout << "Renamed file: " << output_filename << "  -->  " << output_file_path.filename() << std::endl;
            return;
        }

    }
}

bool class_hdf5_file::check_link_exists_recursively(std::string path){
    std::stringstream path_stream(path);
    std::vector<std::string> split_path;
    std::string part;
    while(std::getline(path_stream, part, '/')){
        split_path.push_back(part);
    }
    bool exists = true;
    part.clear();
    for(unsigned long i = 0; i < split_path.size(); i++){
        part += "/" + split_path[i];
        exists = exists and H5Lexists(file, part.c_str(), H5P_DEFAULT);
        if (not exists){break;}
    }
    return exists;
}

void class_hdf5_file::create_dataset_link(const DatasetProperties &props){
    if (not check_link_exists_recursively(props.dset_name)){
        hid_t dataspace = get_DataSpace_unlimited(props.ndims);
        hid_t dset_cpl  = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(dset_cpl, H5D_CHUNKED);
        H5Pset_chunk(dset_cpl, props.ndims, props.chunk_size.data());
//        H5Pset_deflate (dset_cpl ,props.compression_level);
        hid_t dataset = H5Dcreate(file,
                                  props.dset_name.c_str(),
                                  props.datatype,
                                  dataspace,
                                  plist_lncr,
                                  dset_cpl,
                                  H5P_DEFAULT);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Pclose(dset_cpl);
    }
}

void class_hdf5_file::create_group_link(const std::string &group_relative_name) {
    //Check if group exists already
    if (not check_link_exists_recursively(group_relative_name)) {
        hid_t group = H5Gcreate(file, group_relative_name.c_str(), plist_lncr, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(group);
    }
}

void class_hdf5_file::select_hyperslab(const hid_t &filespace, const hid_t &memspace){
    const int ndims = H5Sget_simple_extent_ndims(filespace);
    std::vector<hsize_t> mem_dims(ndims);
    std::vector<hsize_t> file_dims(ndims);
    std::vector<hsize_t> start(ndims);

    H5Sget_simple_extent_dims(memspace , mem_dims.data(), nullptr);
    H5Sget_simple_extent_dims(filespace, file_dims.data(), nullptr);
    for(int i = 0; i < ndims; i++){
        start[i] = file_dims[i] - mem_dims[i];
    }
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start.data(), nullptr, mem_dims.data(), nullptr);
}

void class_hdf5_file::set_extent_dataset(const DatasetProperties &props){
    if (check_link_exists_recursively(props.dset_name)) {
        hid_t dataset = H5Dopen(file, props.dset_name.c_str(), H5P_DEFAULT);
        H5Dset_extent(dataset, props.dims.data());
        H5Dclose(dataset);
    }else{
        std::cerr << "Link does not exist, yet the extent is being set." << std::endl;
        exit(1);
    }
}

void class_hdf5_file::extend_dataset(const std::string & dataset_relative_name, const int dim, const int extent){
    if (H5Lexists(file, dataset_relative_name.c_str(), H5P_DEFAULT)) {
        hid_t dataset = H5Dopen(file, dataset_relative_name.c_str(), H5P_DEFAULT);
        // Retrieve the current size of the dataspace (act as if you don't know it's size and want to append)
        hid_t filespace = H5Dget_space(dataset);
        const int ndims = H5Sget_simple_extent_ndims(filespace);
        std::vector<hsize_t> old_dims(ndims);
        std::vector<hsize_t> new_dims(ndims);
        H5Sget_simple_extent_dims(filespace, old_dims.data(), nullptr);
        new_dims = old_dims;
        new_dims[dim] += extent;
        H5Dset_extent(dataset, new_dims.data());
        H5Dclose(dataset);
        H5Sclose(filespace);
    }

}