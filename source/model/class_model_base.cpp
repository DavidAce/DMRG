//
// Created by david on 2018-07-04.
//

#include "class_model_base.h"
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
#include <math/nmspc_random.h>
#include <iomanip>

using namespace qm;
using Scalar = std::complex<double>;

const Eigen::Tensor<Scalar,4> & class_model_base::MPO() const{
    if (all_mpo_parameters_have_been_set){
        return mpo_internal;
    }else{
        throw std::runtime_error("All MPO parameters haven't been set yet.");
    }
}

bool class_model_base::isReal()const{
    return Textra::isReal(MPO(),"MPO");
}

bool class_model_base::isReduced() const{
    return e_reduced != 0.0;
}

double class_model_base::get_reduced_energy() const {
    return e_reduced;
}

void   class_model_base::set_reduced_energy(double site_energy){
    e_reduced    = site_energy;
    mpo_internal = MPO_reduced_view();
}



class_model_base::class_model_base(size_t position_, std::string logName){
    position = position_;
    log = Logger::setLogger(logName);
}

void class_model_base::set_position(size_t position_){
    position = position_;
}



size_t class_model_base::get_position() const{
    if (position){return position.value();}
    else{throw std::runtime_error("Position of MPO has not been set");}
}

Eigen::MatrixXcd class_model_base::MPO_matrix_view(){
    auto rows = MPO().dimension(0)*MPO().dimension(2);
    auto cols = MPO().dimension(1)*MPO().dimension(3);
    Eigen::Tensor<Scalar,2> MPO_temp = MPO().shuffle(Textra::array4{0,2,1,3}).reshape(Textra::array2{rows,cols});
    return Eigen::Map<Eigen::MatrixXcd>(MPO_temp.data(), rows,cols);
//    return Textra::Tensor_to_Matrix_Map<Scalar>(MPO_temp, rows ,cols);
}