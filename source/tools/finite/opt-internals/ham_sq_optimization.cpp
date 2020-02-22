//
// Created by david on 2019-06-24.
//

#include <tools/finite/opt.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <state/class_state_finite.h>
#include <simulation/nmspc_settings.h>
#include <math/arpack_extra/matrix_product_hamiltonian_sq.h>
#include <math/class_eigsolver.h>

Eigen::Tensor<class_state_finite::Scalar,3> tools::finite::opt::internal::ham_sq_optimization(const class_state_finite & state, OptType optType, OptMode optMode, OptSpace optSpace,std::string ritzstring){
    tools::log->trace("Starting hamiltonian squared optimization");
    using Scalar = std::complex<double>;
    using namespace internal;
    using namespace settings::precision;
    using namespace eigutils::eigSetting;
//    optType = OptType::CPLX;
    switch(optType){
        case OptType::REAL: {
            Eigen::Tensor<double,3> theta   = state.get_multitheta().real();
            Eigen::Tensor<double,4> mpo   = state.get_multimpo().real();
            Eigen::Tensor<double,4> env2L = state.get_ENV2L(state.active_sites.front()).block.real();
            Eigen::Tensor<double,4> env2R = state.get_ENV2R(state.active_sites.back()).block.real();
            auto shape_theta  = state.active_dimensions();
            auto shape_mpo    = mpo.dimensions();

            tools::common::profile::t_eig->tic();
            DenseHamiltonianSqProduct<double>  matrix (
                env2L.data(),
                env2R.data(),
                mpo.data(),
                shape_theta,
                shape_mpo,
                settings::threading::num_threads);
            class_eigsolver solver;
            solver.solverConf.eigMaxIter = 5000;
            solver.solverConf.eigThreshold = 1e-6;
            int nev = 1;
            Ritz ritz = stringToRitz("SM");
            solver.eigs_dense(matrix, nev, eig_max_ncv, NAN, Form::SYMMETRIC, ritz, Side::R, true, false,theta.data());
            [[maybe_unused]] auto eigvals = Eigen::Map<const Eigen::VectorXd> (solver.solution.get_eigvals<Form::SYMMETRIC>().data() ,solver.solution.meta.cols);
            [[maybe_unused]] auto eigvecs = Eigen::Map<const Eigen::VectorXd> (solver.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solver.solution.meta.rows);

            Eigen::VectorXcd theta_new (eigvecs.size());
            for(size_t idx = 0; idx < eigvecs.size();idx++)
                theta_new(idx) = eigvecs(idx);
            tools::common::profile::t_eig->toc();
            return  Textra::MatrixTensorMap(theta_new, state.active_dimensions());
        }
        case OptType::CPLX: {
            auto & env2L = state.get_ENV2L(state.active_sites.front());
            auto & env2R = state.get_ENV2R(state.active_sites.back());
            auto mpo   = state.get_multimpo();
            auto shape_theta  = state.active_dimensions();
            auto shape_mpo    = mpo.dimensions();
            tools::common::profile::t_eig->tic();
            DenseHamiltonianSqProduct<Scalar>  matrix (
                env2L.block.data(),
                env2R.block.data(),
                mpo.data(),
                shape_theta,
                shape_mpo,
                settings::threading::num_threads);
            class_eigsolver solver;
            int nev = 1;
            Ritz ritz = stringToRitz("SM");
            solver.eigs_dense(matrix, nev, eig_max_ncv, NAN, Form::SYMMETRIC, ritz, Side::R, true, true);
            [[maybe_unused]] auto eigvals = Eigen::Map<const Eigen::VectorXd>  (solver.solution.get_eigvals<Form::SYMMETRIC>().data() ,solver.solution.meta.cols);
            [[maybe_unused]] auto eigvecs = Eigen::Map<const Eigen::VectorXcd> (solver.solution.get_eigvecs<Type::CPLX, Form::SYMMETRIC>().data(),solver.solution.meta.rows);
            tools::common::profile::t_eig->toc();
            return  Textra::MatrixTensorMap(eigvecs, state.active_dimensions());
        }
    }
}
