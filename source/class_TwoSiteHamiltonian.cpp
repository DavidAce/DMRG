//
// Created by david on 7/17/17.
//

#include <class_TwoSiteHamiltonian.h>
#include <iomanip>
#include <n_model.h>



class_TwoSiteHamiltonian::class_TwoSiteHamiltonian(){
    asMatrix        = Model::Hamiltonian(sites);
    asTensor4       = Eigen::TensorMap<Tensor4> (asMatrix.data(),2,2,2,2);
    W               = Model::W(sites);
    local_dimension = Model::local_dimension;
    pic         = "**";
};

