#include <config/nmspc_settings.h>
#include <math/eig/arpack_solver/matrix_product_dense.h>
#include <math/eig.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/opt.h>
#include <vector>

std::vector<int> tools::finite::opt::internal::generate_size_list(int shape) {
    std::vector<int> nev_list = {8};
    if(shape <= 512) nev_list = {32,128};
    if(512 < shape and shape <= 1024) nev_list = {64,128,256};
    if(1024 < shape and shape <= 2048) nev_list = {8,64,256};
    if(2048 < shape and shape <= 3072)  nev_list = {8,128,256};
    if(3072 < shape and shape <= 4096)  nev_list = {8,128};
    return nev_list;
}

template<typename Scalar>
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
    tools::finite::opt::internal::subspace::find_subspace_part(const MatrixType<Scalar> &H_local, const TensorType<cplx, 3> &multisite_tensor,
                                                               double energy_target, double subspace_error_threshold, OptMode optMode, OptSpace optSpace) {
    tools::log->trace("Finding subspace -- partial");
    // You need to copy the data into StlMatrixProduct, because the PartialPivLU will overwrite the data in H_local otherwise.
    MatrixProductDense<Scalar> hamiltonian(H_local.data(), H_local.rows(), false, eig::Form::SYMM, eig::Side::R);
    hamiltonian.set_shift(energy_target);
    hamiltonian.FactorOP();
    *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_lu"] += *hamiltonian.t_factorOP;
    double time_lu  = hamiltonian.t_factorOP->get_last_interval();
    double time_ham = tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_ham"]->get_last_interval();

    eig::solver solver;
    solver.config.eigThreshold = settings::precision::eig_threshold;
    std::string                        reason = "exhausted";
    Eigen::VectorXd                    eigvals;
    Eigen::MatrixXcd                   eigvecs;
    Eigen::Map<const Eigen::VectorXcd> multisite_vector(multisite_tensor.data(), multisite_tensor.size());
    for(auto nev : generate_size_list(static_cast<int>(multisite_vector.size()))) {
        tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_eig"]->tic();
        solver.config.clear();
        solver.eigs(hamiltonian,nev, -1,eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, energy_target,  eig::Shinv::ON, eig::Vecs::ON, eig::Dephase::OFF);
        eigvals = eig::view::get_eigvals<eig::real>(solver.result);
        eigvecs = eig::view::get_eigvecs<eig::real>(solver.result, eig::Side::R);
        tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_eig"]->toc();

        Eigen::VectorXd overlaps       = (multisite_vector.adjoint() * eigvecs).cwiseAbs().real();
        double          max_overlap    = overlaps.maxCoeff();
        double          min_overlap    = overlaps.minCoeff();
        double          sq_sum_overlap = overlaps.cwiseAbs2().sum();
        double          subspace_error = 1.0 - sq_sum_overlap;
        reports::eigs_add_entry(nev, max_overlap, min_overlap, std::log10(subspace_error),
                                tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_eig"]->get_last_interval(),time_ham, time_lu);
        time_lu  = 0;
        time_ham = 0;
        if(max_overlap > 1.0 + 1e-6) throw std::runtime_error(fmt::format("max_overlap larger than one: {:.16f}",max_overlap));
        if(sq_sum_overlap > 1.0 + 1e-6) throw std::runtime_error(fmt::format("eps larger than one: {:.16f}", sq_sum_overlap));
        if(min_overlap < 0.0) throw std::runtime_error(fmt::format("min_overlap smaller than zero: {:.16f}", min_overlap));
        if(subspace_error < subspace_error_threshold) {
            reason = fmt::format("subspace error is low enough: {:.3e} < threshold {:.3e}", subspace_error, subspace_error_threshold);
            break;
        }
        if(optSpace == OptSpace::SUBSPACE_AND_DIRECT and subspace_error < 1e-3) {
            reason = fmt::format("subspace error sufficient for SUBSPACE_AND_DIRECT mode: {:.3e} < threshold {:.3e}", subspace_error, 1e-3);
            break;
        }
        if(optSpace == OptSpace::SUBSPACE_ONLY and optMode == OptMode::OVERLAP and max_overlap >= 1.0 / std::sqrt(2.0)) {
            reason = fmt::format("Overlap sufficient for SUBSPACE_ONLY OVERLAP mode:  {:.16f} >= threshold {:.16f}", max_overlap , 1.0 / std::sqrt(2.0) );
            break;
        }
    }
    tools::log->debug("Finished partial eigensolver -- reason: {}", reason);
    return std::make_tuple(eigvecs, eigvals);
}


template std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
tools::finite::opt::internal::subspace::find_subspace_part(const MatrixType<real> &H_local, const TensorType<cplx,3> &multisite_tensor, double energy_target, double subspace_error_threshold, OptMode optMode, OptSpace optSpace);

template std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
tools::finite::opt::internal::subspace::find_subspace_part(const MatrixType<cplx> &H_local, const TensorType<cplx,3> &multisite_tensor, double energy_target, double subspace_error_threshold, OptMode optMode, OptSpace optSpace);

