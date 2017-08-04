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