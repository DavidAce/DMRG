//
// Created by david on 7/17/17.
//

#include <class_TwoSiteHamiltonian.h>
#include <iomanip>
#include <n_model.h>
#include <unsupported/Eigen/MatrixFunctions>



class_TwoSiteHamiltonian::class_TwoSiteHamiltonian(){
    asMatrix        = Model::Hamiltonian(sites);
    asTensor4       = Eigen::TensorMap<Tensor4> (asMatrix.data(),2,2,2,2);
    W               = Model::W(sites);
    local_dimension = Model::local_dimension;
    picture         = "**";
    asTimeEvolution = matrix_to_tensor4((-delta * asMatrix).exp().eval(), array4{2,2,2,2});
};

