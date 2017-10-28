//
// Created by david on 7/17/17.
//

#include <class_TwoSiteHamiltonian.h>
#include <n_model.h>


class_TwoSiteHamiltonian::class_TwoSiteHamiltonian(){

    asMatrix       = Model::H();
    asTensor4      = Eigen::TensorMap<Tensor4d> (asMatrix.data(),2,2,2,2);
    W              = Model::W();
    local_dimension= Model::local_dimension;
    picture        = "**";
};

