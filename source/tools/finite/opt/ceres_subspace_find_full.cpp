#include "opt-internal.h"
#include "report.h"
#include <config/enums.h>
#include <config/settings.h>
#include <math/eig.h>
#include <tensors/state/StateFinite.h>
#include <tensors/TensorsFinite.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>

namespace tools::finite::opt::internal {
    template<typename Scalar>
    std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_full(const TensorsFinite &tensors) {
        tools::log->trace("Finding subspace -- full");
        auto t_full = tid::tic_scope("full");

        // Generate the Hamiltonian matrix
        double energy_shift = 0.0;
        if(settings::precision::use_reduced_mpo_energy and settings::precision::use_shifted_mpo_energy)
            energy_shift = tools::finite::measure::energy_minus_energy_reduced(tensors);

        MatrixType<Scalar> H_local       = tools::finite::opt::internal::get_multisite_hamiltonian_matrix<Scalar>(*tensors.model, *tensors.edges, energy_shift);
        const auto        &multisite_mps = tensors.state->get_multisite_mps();
        const auto         multisite_vec = Eigen::Map<const Eigen::VectorXcd>(multisite_mps.data(), multisite_mps.size());

        if(H_local.rows() != H_local.cols())
            throw std::runtime_error(fmt::format("H_local is not square: rows [{}] | cols [{}]", H_local.rows(), H_local.cols()));
        eig::solver solver;
        solver.eig<eig::Form::SYMM>(H_local.data(), H_local.rows(), eig::Vecs::ON, eig::Dephase::OFF);
        auto eigvals = eig::view::get_eigvals<double>(solver.result);
        auto eigvecs = eig::view::get_eigvecs<Scalar>(solver.result);
        tools::log->debug("Finished eigensolver -- reason: Full diagonalization");
        Eigen::VectorXd overlaps = (multisite_vec.adjoint() * eigvecs).cwiseAbs().real();
        int             idx;
        double          max_overlap    = overlaps.maxCoeff(&idx);
        double          min_overlap    = overlaps.minCoeff();
        double          sq_sum_overlap = overlaps.cwiseAbs2().sum();
        double          subspace_error = 1.0 - sq_sum_overlap;
        long            nev            = eigvecs.cols();
        auto            time_eig       = tid::get("eig").get_last_interval();
        auto            time_ham       = tid::get("ham").get_time();
        reports::eigs_add_entry(nev, max_overlap, min_overlap, std::log10(subspace_error), time_eig, time_ham, 0, 1);
        return {eigvecs, eigvals};
    }

    template std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_full<cplx>(const TensorsFinite &tensors);

    template std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_full<real>(const TensorsFinite &tensors);

}
