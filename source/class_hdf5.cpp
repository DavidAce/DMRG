//
// Created by david on 8/1/17.
//

#include "class_hdf5.h"



// Specializations
template<>
void class_hdf5::write_to_hdf5<MatrixType>(const MatrixType &matrix, const std::string dataset_name){
    int rank = 2;
    hsize_t     dims[rank] = {(hsize_t)matrix.rows(),(hsize_t) matrix.cols()};
    hid_t scalar = get_ScalarType(matrix);
    H5LTmake_dataset(file_id,dataset_name.c_str(), rank, dims,scalar,matrix.data());
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
            std::cout << "Searching for directory: " << output_path << std::endl;
            if (fs::exists(output_path)) {
                output_path = fs::canonical(output_path);
                std::cout << "Found " << output_dirname << " directory: " << output_path << std::endl;
                break;
            } else if (output_path.parent_path().has_parent_path()) {
                output_path = output_path.parent_path().parent_path() / output_dirname.stem();
            } else {
                std::cout << "ERROR: " << output_dirname << " folder does not exist and create_dir == false"
                          << std::endl;
                std::cout << "Either create an " << output_dirname << "  directory manually or pass create_dir = true."
                          << std::endl;
                exit(1);
            }
        }
    }
}
