#include "opt-internal.h"
#include "report.h"
#include <config/settings.h>
#include <math/eig.h>
#include <math/eig/matvec/matvec_dense.h>
#include <tensors/state/StateFinite.h>
#include <tensors/TensorsFinite.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>
#include <vector>

namespace tools::finite::opt::internal {

    std::vector<int> generate_size_list(int shape) {
        std::vector<int> nev_list = {8};
        if(shape <= 512) nev_list = {32, 128};
        if(512 < shape and shape <= 1024) nev_list = {64, 128, 256};
        if(1024 < shape and shape <= 2048) nev_list = {8, 64, 256};
        if(2048 < shape and shape <= 3072) nev_list = {8, 64, 128};
        if(3072 < shape and shape <= 4096) nev_list = {8, 64};
        return nev_list;
    }

    template<typename Scalar>
    std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_iter(const TensorsFinite &tensors, double energy_target,
                                                                               double subspace_error_threshold, OptMode optMode, OptSpace optSpace) {
        tools::log->trace("Finding subspace -- iterative");
        auto  t_iter = tid::tic_scope("iter");
        auto &t_lu   = tid::get("lu_decomp");

        // Generate the Hamiltonian matrix
        double energy_shift = 0.0;
        if(settings::precision::use_reduced_mpo_energy and settings::precision::use_shifted_mpo_energy)
            energy_shift = tools::finite::measure::energy_minus_energy_reduced(tensors);

        MatrixType<Scalar> H_local       = tools::finite::opt::internal::get_multisite_hamiltonian_matrix<Scalar>(*tensors.model, *tensors.edges, energy_shift);
        double             time_ham      = tid::get("ham").get_time();
        const auto        &multisite_mps = tensors.state->get_multisite_mps();

        // Define the dense matrix for the eigenvalue solver
        MatVecDense<Scalar> hamiltonian(H_local.data(), H_local.rows(), false, eig::Form::SYMM, eig::Side::R);
#pragma message "Do we really need to set shift a second time?"
        hamiltonian.set_shift(energy_target);
        eig::solver solver;
        solver.config.tol                         = settings::precision::eig_tolerance;
        std::string                        reason = "exhausted";
        Eigen::VectorXd                    eigvals;
        Eigen::MatrixXcd                   eigvecs;
        Eigen::Map<const Eigen::VectorXcd> multisite_vector(multisite_mps.data(), multisite_mps.size());
        for(auto nev : generate_size_list(static_cast<int>(multisite_vector.size()))) {
            hamiltonian.FactorOP();
            t_lu += *hamiltonian.t_factorOP;

            solver.config.clear();
            solver.eigs(hamiltonian, nev, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, energy_target, eig::Shinv::ON, eig::Vecs::ON, eig::Dephase::OFF);
            eigvals = eig::view::get_eigvals<eig::real>(solver.result);
            eigvecs = eig::view::get_eigvecs<eig::cplx>(solver.result, eig::Side::R);

            Eigen::VectorXd overlaps       = (multisite_vector.adjoint() * eigvecs).cwiseAbs().real();
            double          max_overlap    = overlaps.maxCoeff();
            double          min_overlap    = overlaps.minCoeff();
            double          sq_sum_overlap = overlaps.cwiseAbs2().sum();
            double          subspace_error = 1.0 - sq_sum_overlap;
            reports::eigs_add_entry(nev, max_overlap, min_overlap, std::log10(subspace_error), solver.result.meta.time_total, time_ham,
                                    t_lu.get_last_interval(), static_cast<size_t>(solver.result.meta.counter));
            time_ham = 0;
            if(max_overlap > 1.0 + 1e-6) throw std::runtime_error(fmt::format("max_overlap larger than one: {:.16f}", max_overlap));
            if(sq_sum_overlap > 1.0 + 1e-6) throw std::runtime_error(fmt::format("eps larger than one: {:.16f}", sq_sum_overlap));
            if(min_overlap < 0.0) throw std::runtime_error(fmt::format("min_overlap smaller than zero: {:.16f}", min_overlap));
            if(subspace_error < subspace_error_threshold) {
                reason = fmt::format("subspace error is low enough: {:.3e} < threshold {:.3e}", subspace_error, subspace_error_threshold);
                break;
            }
            if(optSpace == OptSpace::SUBSPACE_AND_DIRECT and subspace_error < 1e-3) {
                reason = fmt::format("subspace error is low enough for SUBSPACE_AND_DIRECT mode: {:.3e} < threshold {:.3e}", subspace_error, 1e-3);
                break;
            }
            if(optSpace == OptSpace::SUBSPACE_ONLY and optMode == OptMode::OVERLAP and max_overlap >= 1.0 / std::sqrt(2.0)) {
                reason = fmt::format("Overlap is sufficient for SUBSPACE_ONLY OVERLAP mode:  {:.16f} >= threshold {:.16f}", max_overlap, 1.0 / std::sqrt(2.0));
                break;
            }
        }
        tools::log->debug("Finished iterative eigensolver -- reason: {}", reason);
        return {eigvecs, eigvals};
    }

    template std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_iter<cplx>(const TensorsFinite &tensors, double energy_target,
                                                                                              double subspace_error_threshold, OptMode optMode,
                                                                                              OptSpace optSpace);

    template std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_iter<real>(const TensorsFinite &tensors, double energy_target,
                                                                                              double subspace_error_threshold, OptMode optMode,
                                                                                              OptSpace optSpace);

}
