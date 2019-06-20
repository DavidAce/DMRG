//
// Created by david on 2018-07-04.
//

#include "class_hamiltonian_base.h"
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_random_numbers.h>
#include <iomanip>

using namespace qm;
using Scalar = std::complex<double>;

const Eigen::Tensor<Scalar,4> & class_hamiltonian_base::MPO() const{
    if (all_mpo_parameters_have_been_set){
        return mpo_internal;
    }else{
        throw std::runtime_error("All MPO parameters haven't been set yet.");
    }
}

bool class_hamiltonian_base::isReal()const{
    return Textra::isReal(MPO(),"MPO");
}

class_hamiltonian_base::class_hamiltonian_base(std::string logName){
    log = Logger::setLogger(logName);
}

void class_hamiltonian_base::set_position(size_t position_){
    position = position_;
}



size_t class_hamiltonian_base::get_position() const{
    return position;
}

Eigen::MatrixXcd class_hamiltonian_base::MPO_matrix_view(){
    auto rows = MPO().dimension(0)*MPO().dimension(2);
    auto cols = MPO().dimension(1)*MPO().dimension(3);
    Eigen::Tensor<Scalar,4> MPO_temp = MPO().shuffle(Textra::array4{0,2,1,3});
    return Textra::Tensor_to_Matrix<Scalar>(MPO_temp, rows ,cols);
}