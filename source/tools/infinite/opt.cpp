//
// Created by david on 2019-07-06.
//

#include <tools/nmspc_tools.h>
#include <state/class_infinite_state.h>
#include <state/class_environment.h>
#include <state/class_mps_2site.h>
#include <state/class_vidal_site.h>
#include <model/class_model_base.h>
#include <math/class_eigsolver.h>
#include <math/class_svd_wrapper.h>
#include <math/arpack_extra/matrix_product_hamiltonian.h>


Eigen::Tensor<tools::infinite::Scalar,4>
        tools::infinite::opt::find_ground_state(const class_infinite_state & state, std::string ritzstring){
    auto theta = state.get_theta();
    std::array<long,4> shape_theta4 = theta.dimensions();
    std::array<long,4> shape_mpo4   = state.HA->MPO().dimensions();

    tools::common::profile::t_eig.tic();
    int nev = 1;
    using namespace settings::precision;
    using namespace eigutils::eigSetting;
    Ritz ritz = stringToRitz(ritzstring);

    DenseHamiltonianProduct<Scalar>  matrix (state.Lblock->block.data(), state.Rblock->block.data(), state.HA->MPO().data(), state.HB->MPO().data(), shape_theta4, shape_mpo4);
    class_eigsolver solver;
    solver.eigs_dense(matrix,nev,eigMaxNcv,NAN,Form::SYMMETRIC,ritz,Side::R,true,true);
    auto eigvec  = Eigen::TensorMap<const Eigen::Tensor<Scalar,1>>  (solver.solution.get_eigvecs<Type::CPLX, Form::SYMMETRIC>().data(),solver.solution.meta.rows);

    tools::common::profile::t_eig.toc();
    tools::common::profile::t_eig.print_delta();

    return eigvec.reshape(theta.dimensions());

}


//============================================================================//
// Do unitary evolution on an MPS
//============================================================================//


Eigen::Tensor<tools::infinite::Scalar, 4> tools::infinite::opt::time_evolve_theta(const class_infinite_state & state, const Eigen::Tensor<Scalar, 4> &U)
/*!
@verbatim
  1--[ Θ ]--3
     |   |
     0   2
                   1--[ Θ ]--3
     0   1   --->     |   |
     |   |            0   2
     [ U ]
     |   |
     2   3
@endverbatim
*/
{
    return U.contract(state.get_theta(), Textra::idx({0,1},{0,2}))
            .shuffle(Textra::array4{0,2,1,3});
}





void tools::infinite::opt::truncate_theta(Eigen::Tensor<Scalar,4> &theta, class_infinite_state & state, long chi_, double SVDThreshold){
    class_SVD SVD;
    SVD.setThreshold(SVDThreshold);
    auto[U, S, V] = SVD.schmidt(theta,chi_);
    state.MPS->truncation_error = SVD.get_truncation_error();
    state.MPS->LC  = S;
    Eigen::Tensor<Scalar,3> L_U = Textra::asDiagonalInversed(state.MPS->MPS_A->get_L()).contract(U,Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
    Eigen::Tensor<Scalar,3> V_L = V.contract(Textra::asDiagonalInversed(state.MPS->MPS_B->get_L()), Textra::idx({2},{0}));
    state.MPS->MPS_A->set_G(L_U);
    state.MPS->MPS_B->set_G(V_L);
}



