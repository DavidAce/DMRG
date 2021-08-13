#include "../opt_meta.h"
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
    std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_iter(const TensorsFinite &tensors, double energy_target, double target_subspace_error,
                                                                              const OptMeta &conf) {
        tools::log->trace("Finding subspace -- iterative");
        auto  t_iter = tid::tic_scope("iter");
        auto &t_lu   = tid::get("lu_decomp");

        // Generate the Hamiltonian matrix
        MatrixType<Scalar> H_local       = tools::finite::opt::internal::get_multisite_hamiltonian_matrix<Scalar>(*tensors.model, *tensors.edges);
        double             time_ham      = tid::get("ham").get_last_interval();
        const auto        &multisite_mps = tensors.state->get_multisite_mps();

        // Define the dense matrix for the eigenvalue solver
        MatVecDense<Scalar>      hamiltonian(H_local.data(), H_local.rows(), false, eig::Form::SYMM, eig::Side::R);
        Eigen::Tensor<Scalar, 3> init;
        if constexpr(std::is_same_v<Scalar, real>) init = tensors.state->get_multisite_mps().real();
        if constexpr(std::is_same_v<Scalar, cplx>) init = tensors.state->get_multisite_mps();

        eig::settings config;
        config.tol             = settings::precision::eig_tolerance;
        config.shift_invert    = eig::Shinv::ON;
        config.sigma           = energy_target;
        config.compute_eigvecs = eig::Vecs::ON;
        config.ritz            = eig::Ritz::LM;
        config.initial_guess.push_back({init.data(), 0});
        std::string                        reason = "exhausted";
        Eigen::VectorXd                    eigvals;
        Eigen::MatrixXcd                   eigvecs;
        Eigen::Map<const Eigen::VectorXcd> multisite_vector(multisite_mps.data(), multisite_mps.size());
        for(auto nev : generate_size_list(static_cast<int>(multisite_vector.size()))) {
            eig::solver solver;
            solver.config        = config;
            solver.config.maxNev = nev;
            // Set the new initial guess if we are doing one more round
            if(eigvecs.cols() != 0) {
                solver.config.initial_guess.clear();
                for(long n = 0; n < eigvecs.cols(); n++) { config.initial_guess.push_back({eigvecs.col(n).data(), n}); }
            }


            solver.eigs(hamiltonian);
            t_lu += *hamiltonian.t_factorOP;

            eigvals = eig::view::get_eigvals<eig::real>(solver.result);
            eigvecs = eig::view::get_eigvecs<eig::cplx>(solver.result, eig::Side::R);

            Eigen::VectorXd overlaps       = (multisite_vector.adjoint() * eigvecs).cwiseAbs().real();
            double          max_overlap    = overlaps.maxCoeff();
            double          min_overlap    = overlaps.minCoeff();
            double          sq_sum_overlap = overlaps.cwiseAbs2().sum();
            double          subspace_error = 1.0 - sq_sum_overlap;
            reports::eigs_add_entry(nev, max_overlap, min_overlap, std::log10(subspace_error), solver.result.meta.time_total, time_ham,
                                    t_lu.get_last_interval(), static_cast<size_t>(solver.result.meta.matvecs));
            time_ham = 0;
            if(max_overlap > 1.0 + 1e-6) throw std::runtime_error(fmt::format("max_overlap larger than one: {:.16f}", max_overlap));
            if(sq_sum_overlap > 1.0 + 1e-6) throw std::runtime_error(fmt::format("eps larger than one: {:.16f}", sq_sum_overlap));
            if(min_overlap < 0.0) throw std::runtime_error(fmt::format("min_overlap smaller than zero: {:.16f}", min_overlap));
            if(subspace_error < target_subspace_error) {
                reason = fmt::format("subspace error is low enough: {:.3e} < threshold {:.3e}", subspace_error, target_subspace_error);
                break;
            }
            if(conf.optSpace == OptSpace::SUBSPACE and conf.optMode == OptMode::OVERLAP and max_overlap >= 1.0 / std::sqrt(2.0)) {
                reason = fmt::format("Overlap is sufficient for SUBSPACE OVERLAP mode:  {:.16f} >= threshold {:.16f}", max_overlap, 1.0 / std::sqrt(2.0));
                break;
            }
        }
        tools::log->debug("Finished iterative eigensolver -- reason: {}", reason);
        return {eigvecs, eigvals};
    }

    template std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_iter<cplx>(const TensorsFinite &tensors, double energy_target,
                                                                                             double target_subspace_error, const OptMeta &conf);

    template std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_iter<real>(const TensorsFinite &tensors, double energy_target,
                                                                                             double target_subspace_error, const OptMeta &conf);

    template<typename Scalar>
    std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_iter2(const TensorsFinite &tensors, double energy_target, double target_subspace_error,
                                                                               const OptMeta &conf) {
        tools::log->trace("Finding subspace -- iterative");
        auto  t_iter = tid::tic_scope("iter");
        auto &t_lu   = tid::get("lu_decomp");

        // Generate the Hamiltonian matrix
        MatrixType<Scalar> H_local       = tools::finite::opt::internal::get_multisite_hamiltonian_matrix<Scalar>(*tensors.model, *tensors.edges);
        double             time_ham      = tid::get("ham").get_last_interval();
        const auto        &multisite_mps = tensors.state->get_multisite_mps();

        // Define the dense matrix for the eigenvalue solver
        MatVecDense<Scalar>      hamiltonian(H_local.data(), H_local.rows(), false, eig::Form::SYMM, eig::Side::R);
        Eigen::Tensor<Scalar, 3> init;
        if constexpr(std::is_same_v<Scalar, real>) init = tensors.state->get_multisite_mps().real();
        if constexpr(std::is_same_v<Scalar, cplx>) init = tensors.state->get_multisite_mps();

        eig::settings config;
        config.lib           = eig::Lib::PRIMME;
        config.primme_method = "PRIMME_DYNAMIC";
        //        config.tol = settings::precision::eig_tolerance;
        config.tol               = 1e-4;
        config.ritz              = eig::Ritz::primme_closest_abs;
        config.sigma             = energy_target;
        config.compute_eigvecs   = eig::Vecs::ON;
        config.shift_invert      = eig::Shinv::OFF;
        config.maxIter           = 80000;
        config.primme_projection = "primme_proj_refined";
        config.primme_locking    = 1;
        config.initial_guess.push_back({init.data(), 0});
        std::string                        reason = "exhausted";
        Eigen::VectorXd                    eigvals;
        Eigen::MatrixXcd                   eigvecs;
        Eigen::Map<const Eigen::VectorXcd> multisite_vector(multisite_mps.data(), multisite_mps.size());
        for(auto nev : generate_size_list(static_cast<int>(multisite_vector.size()))) {
            eig::solver solver;
            solver.config        = config;
            solver.config.maxNev = nev;
            solver.config.maxNcv = std::clamp(nev * 4, nev, hamiltonian.rows());
            // Set the new initial guess if we are doing one more round
            if(eigvecs.cols() != 0) {
                solver.config.initial_guess.clear();
                for(long n = 0; n < eigvecs.cols(); n++) { config.initial_guess.push_back({eigvecs.col(n).data(), n}); }
            }

            solver.eigs(hamiltonian);
            t_lu += *hamiltonian.t_factorOP;

            eigvals = eig::view::get_eigvals<eig::real>(solver.result);
            eigvecs = eig::view::get_eigvecs<eig::cplx>(solver.result, eig::Side::R);

            Eigen::VectorXd overlaps       = (multisite_vector.adjoint() * eigvecs).cwiseAbs().real();
            double          max_overlap    = overlaps.maxCoeff();
            double          min_overlap    = overlaps.minCoeff();
            double          sq_sum_overlap = overlaps.cwiseAbs2().sum();
            double          subspace_error = 1.0 - sq_sum_overlap;
            reports::eigs_add_entry(nev, max_overlap, min_overlap, std::log10(subspace_error), solver.result.meta.time_total, time_ham,
                                    t_lu.get_last_interval(), static_cast<size_t>(solver.result.meta.matvecs));
            time_ham = 0;
            if(max_overlap > 1.0 + 1e-6) throw std::runtime_error(fmt::format("max_overlap larger than one: {:.16f}", max_overlap));
            if(sq_sum_overlap > 1.0 + 1e-6) throw std::runtime_error(fmt::format("eps larger than one: {:.16f}", sq_sum_overlap));
            if(min_overlap < 0.0) throw std::runtime_error(fmt::format("min_overlap smaller than zero: {:.16f}", min_overlap));
            if(subspace_error < target_subspace_error) {
                reason = fmt::format("subspace error is low enough: {:.3e} < threshold {:.3e}", subspace_error, target_subspace_error);
                break;
            }
            if(conf.optSpace == OptSpace::SUBSPACE and conf.optMode == OptMode::OVERLAP and max_overlap >= 1.0 / std::sqrt(2.0)) {
                reason = fmt::format("Overlap is sufficient for SUBSPACE OVERLAP mode:  {:.16f} >= threshold {:.16f}", max_overlap, 1.0 / std::sqrt(2.0));
                break;
            }
        }
        tools::log->debug("Finished iterative eigensolver -- reason: {}", reason);
        return {eigvecs, eigvals};
    }

    template std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_iter2<cplx>(const TensorsFinite &tensors, double energy_target,
                                                                                              double target_subspace_error, const OptMeta &conf);

    template std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_iter2<real>(const TensorsFinite &tensors, double energy_target,
                                                                                              double target_subspace_error, const OptMeta &conf);

}
