//
// Created by david on 2019-06-24.
//

#include <mps_tools/finite/opt.h>
#include <mps_state/class_finite_chain_state.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <general/arpack_extra/matrix_product_hamiltonian.h>
#include <general/class_eigsolver.h>

std::tuple<Eigen::Tensor<std::complex<double>,4>, double> mpstools::finite::opt::internals::ground_state_optimization(const class_finite_chain_state & state, const class_simulation_state & sim_state, std::string ritz){
    mpstools::log->trace("Starting ground state optimization");
    using Scalar = std::complex<double>;
    using namespace internals;
    using namespace settings::precision;
    using namespace eigutils::eigSetting;

    Ritz ritz_enum = Ritz::SR;
    if(ritz == "SR"){ritz_enum = Ritz::SR;}
    if(ritz == "LR"){ritz_enum = Ritz::LR;}

    auto theta = state.get_theta();
    std::array<long,4> shape_theta4  = theta.dimensions();
    std::array<long,4> shape_mpo4   = state.MPO_L.back()->MPO().dimensions();

    t_eig->tic();
    int nev = 1;
    DenseHamiltonianProduct<Scalar>  matrix (
            state.ENV_L.back().block.data(),
            state.ENV_R.front().block.data(),
            state.MPO_L.back()->MPO().data(),
            state.MPO_R.front()->MPO().data(),
            shape_theta4,
            shape_mpo4);



    class_eigsolver solver;
    solver.eigs_dense(matrix,nev,eigMaxNcv,NAN,Form::SYMMETRIC,ritz_enum,Side::R,true,true);

    auto eigvals           = Eigen::TensorMap<const Eigen::Tensor<double,1>>  (solver.solution.get_eigvals<Form::SYMMETRIC>().data() ,solver.solution.meta.cols);
    auto eigvecs           = Eigen::TensorMap<const Eigen::Tensor<Scalar,1>>  (solver.solution.get_eigvecs<Type::CPLX, Form::SYMMETRIC>().data(),solver.solution.meta.rows);

    t_eig->toc();
//    t_eig->print_delta();

    return std::make_tuple(eigvecs.reshape(theta.dimensions()), std::real(eigvals(0))/state.get_length());
}
