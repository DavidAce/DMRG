//
// Created by david on 7/17/17.
//

#include <class_TwoSiteHamiltonian.h>
#include <n_model.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <n_settings.h>


class_TwoSiteHamiltonian::class_TwoSiteHamiltonian(){
//    asMatrix        = Model::Hamiltonian(sites);
    asMatrix       = Model::Hamiltonian2(sites);
    asTensor4       = Eigen::TensorMap<Tensor4> (asMatrix.data(),2,2,2,2);
    W               = Model::W(sites);
    local_dimension = Model::local_dimension;
    picture         = "**";
//    asTimeEvolution = matrix_to_tensor<4>((-settings::itebd_delta_t * asMatrix).exp().eval(), array4{2,2,2,2});
//    asTimeEvolution = matrix_to_tensor<4>((-settings::itebd_delta_t * asMatrix).exp().eval(), array4{2,2,2,2});
    asTimeEvolution = Model::TimeEvolution_4th_order(sites, settings::itebd::delta_t);

};

