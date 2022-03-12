#include "math/tenx.h"
// -- (textra first)
#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "math/eig.h"
#include "math/eig/matvec/matvec_dense.h"
#include "math/eig/matvec/matvec_mpo.h"
#include "tensors/edges/EdgesFinite.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/state/StateFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/contraction.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include "tools/finite/opt/opt-internal.h"
#include "tools/finite/opt/report.h"
#include "tools/finite/opt_meta.h"
#include "tools/finite/opt_mps.h"

//
#include "tools/finite/opt/bfgs_callback.h"
#include "tools/finite/opt/bfgs_simps_functor.h"
#include <ceres/gradient_problem.h>
#include <Eigen/QR>
#include <primme/primme.h>

using namespace tools::finite::opt;
using namespace tools::finite::opt::internal;

std::vector<int> subspace::generate_nev_list(int rows) {
    std::vector<int> nev_list = {4, 8};
    if(32 < rows and rows <= 64) nev_list = {16, 32};
    if(64 < rows and rows <= 128) nev_list = {32, 64};
    if(128 < rows and rows <= 512) nev_list = {64, 256};
    if(512 < rows and rows <= 1024) nev_list = {64, 256, 512};
    if(1024 < rows and rows <= 2048) nev_list = {16, 256, 512};
    if(2048 < rows and rows <= 3072) nev_list = {16, 256};
    if(3072 < rows and rows <= 4096) nev_list = {16};
    if(4096 < rows) nev_list = {4};

    while(nev_list.size() > 1 and (nev_list.back() * 2 > rows or static_cast<size_t>(nev_list.back()) > settings::precision::max_subspace_size))
        nev_list.pop_back();
    if(nev_list.empty()) throw std::logic_error("nev_list is empty");
    return nev_list;
}

template<typename Scalar>
std::vector<opt_mps> subspace::find_subspace(const TensorsFinite &tensors, double target_subspace_error, const OptMeta &meta) {
    const auto &state = *tensors.state;
    const auto &model = *tensors.model;
    tools::log->trace("Finding subspace");
    auto t_find     = tid::tic_scope("find");
    auto dbl_length = static_cast<double>(state.get_length());

    Eigen::MatrixXcd eigvecs;
    Eigen::VectorXd  eigvals;

    // If the mps is small enough you can afford full diag.
    if(tensors.state->active_problem_size() <= settings::precision::max_size_full_diag) {
        std::tie(eigvecs, eigvals) = find_subspace_full<Scalar>(tensors);
    } else {
        double eigval_target;
        double energy_target = tools::finite::measure::energy(tensors);
        if(model.is_shifted()) {
            eigval_target = tools::finite::measure::energy_minus_energy_shift(tensors);
            tools::log->trace("Energy shift  = {:.16f} | per site = {:.16f}", model.get_energy_shift(), model.get_energy_shift_per_site());
            tools::log->trace("Energy target = {:.16f} | per site = {:.16f}", energy_target, energy_target / dbl_length);
            tools::log->trace("Eigval target = {:.16f} | per site = {:.16f}", eigval_target, eigval_target / dbl_length);
            tools::log->trace("Eigval target + Energy shift = Energy: {:.16f} + {:.16f} = {:.16f}", eigval_target / dbl_length,
                              model.get_energy_shift_per_site(), energy_target / dbl_length);
        } else {
            eigval_target = energy_target;
        }
        std::tie(eigvecs, eigvals) = find_subspace_part<Scalar>(tensors, eigval_target, target_subspace_error, meta);
    }
    /* clang-format off */
    tools::log->trace("Eigval range         : {:.16f} --> {:.16f}", eigvals.minCoeff(), eigvals.maxCoeff());
    tools::log->trace("Energy range         : {:.16f} --> {:.16f}", eigvals.minCoeff() + model.get_energy_shift(), eigvals.maxCoeff() + model.get_energy_shift());
    tools::log->trace("Energy range per site: {:.16f} --> {:.16f}", eigvals.minCoeff() / dbl_length + model.get_energy_shift_per_site(), eigvals.maxCoeff() / dbl_length + model.get_energy_shift_per_site());
    /* clang-format on */
    reports::print_subs_report();

    if constexpr(std::is_same<Scalar, double>::value) {
        tenx::subtract_phase(eigvecs);
        double trunc = eigvecs.imag().cwiseAbs().sum();
        if(trunc > 1e-12) tools::log->warn("truncating imag of eigvecs, sum: {}", trunc);
        eigvecs = eigvecs.real();
    }
    const auto     &multisite_mps = state.get_multisite_mps();
    const auto      multisite_vec = Eigen::Map<const Eigen::VectorXcd>(multisite_mps.data(), multisite_mps.size());
    auto            energy_shift  = model.get_energy_shift();
    Eigen::VectorXd overlaps      = (multisite_vec.adjoint() * eigvecs).cwiseAbs().real();

    double eigvec_time = 0;
    for(const auto &item : reports::subs_log) { eigvec_time += item.ham_time + item.lu_time + item.eig_time; }

    std::vector<opt_mps> subspace;
    subspace.reserve(static_cast<size_t>(eigvals.size()));
    for(long idx = 0; idx < eigvals.size(); idx++) {
        // Important to normalize the eigenvectors that we get from the solver: they are not always well normalized when we get them!
        auto eigvec_i = tenx::TensorCast(eigvecs.col(idx).normalized(), multisite_mps.dimensions());
        subspace.emplace_back(fmt::format("eigenvector {}", idx), eigvec_i, tensors.active_sites, eigvals(idx), energy_shift, std::nullopt, overlaps(idx),
                              tensors.get_length());
        subspace.back().is_basis_vector = true;
        subspace.back().set_time(eigvec_time);
        subspace.back().set_mv(reports::subs_log.size());
        subspace.back().set_iter(reports::subs_log.size());
        subspace.back().set_eigs_idx(idx);
        subspace.back().set_eigs_eigval(eigvals(idx));
        subspace.back().set_eigs_ritz(enum2sv(meta.optRitz));
        subspace.back().set_optmode(meta.optMode);
        subspace.back().set_optsolver(meta.optSolver);
        subspace.back().validate_basis_vector();
    }
    return subspace;
}

template std::vector<opt_mps> subspace::find_subspace<cplx>(const TensorsFinite &tensors, double target_subspace_error, const OptMeta &meta);

template std::vector<opt_mps> subspace::find_subspace<real>(const TensorsFinite &tensors, double target_subspace_error, const OptMeta &meta);

template<typename Scalar>
std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_part(const TensorsFinite &tensors, double energy_target, double target_subspace_error,
                                                                          const OptMeta &meta) {
    tools::log->trace("Finding subspace -- partial");
    auto  t_iter = tid::tic_scope("part");
    auto &t_lu   = tid::get("lu_decomp");

    // Initial mps and a vector map
    const auto                        &multisite_mps = tensors.state->get_multisite_mps();
    Eigen::Map<const Eigen::VectorXcd> multisite_vector(multisite_mps.data(), multisite_mps.size());
    const auto                         problem_size = multisite_mps.size();

    // Mutable initial mps vector used for initial guess in arpack
    Eigen::Tensor<Scalar, 3> init;
    if constexpr(std::is_same_v<Scalar, real>) init = tensors.state->get_multisite_mps().real();
    if constexpr(std::is_same_v<Scalar, cplx>) init = tensors.state->get_multisite_mps();

    // Get the local effective Hamiltonian as a matrix
    const auto &effective_hamiltonian = tensors.get_effective_hamiltonian<Scalar>();
    double      time_ham              = tid::get("ham").get_last_interval();

    // Create the dense matrix object for the eigenvalue solver
    MatVecDense<Scalar> hamiltonian(effective_hamiltonian.data(), effective_hamiltonian.dimension(0), false, eig::Form::SYMM, eig::Side::R);

    // Create a reusable config for multiple nev trials
    eig::settings config;
    config.tol             = settings::precision::eigs_tolerance;
    config.sigma           = energy_target;
    config.shift_invert    = eig::Shinv::ON;
    config.compute_eigvecs = eig::Vecs::ON;
    config.ritz            = eig::Ritz::LM;
    config.initial_guess.push_back({init.data(), 0});
    std::string reason = "exhausted";

    // Initialize eigvals/eigvecs containers that store the results
    Eigen::VectorXd  eigvals;
    Eigen::MatrixXcd eigvecs;
    for(auto nev : generate_nev_list(static_cast<int>(problem_size))) {
        eig::solver solver;
        solver.config        = config;
        solver.config.maxNev = nev;
        // Set the new initial guess if we are doing one more round
        if(eigvecs.cols() != 0) {
            solver.config.initial_guess.clear();
            for(long n = 0; n < eigvecs.cols(); n++) { solver.config.initial_guess.push_back({eigvecs.col(n).data(), n}); }
        }

        solver.eigs(hamiltonian);
        t_lu += *hamiltonian.t_factorOP;

        eigvals = eig::view::get_eigvals<eig::real>(solver.result);
        eigvecs = eig::view::get_eigvecs<eig::cplx>(solver.result, eig::Side::R);

        // Check the quality of the subspace
        Eigen::VectorXd overlaps       = (multisite_vector.adjoint() * eigvecs).cwiseAbs().real();
        double          max_overlap    = overlaps.maxCoeff();
        double          min_overlap    = overlaps.minCoeff();
        double          sq_sum_overlap = overlaps.cwiseAbs2().sum();
        double          subspace_error = 1.0 - sq_sum_overlap;
        reports::subs_add_entry(nev, max_overlap, min_overlap, subspace_error, solver.result.meta.time_total, time_ham, t_lu.get_last_interval(),
                                solver.result.meta.iter, solver.result.meta.num_mv, solver.result.meta.num_pc);
        time_ham = 0;
        if(max_overlap > 1.0 + 1e-6) throw std::runtime_error(fmt::format("max_overlap larger than one: {:.16f}", max_overlap));
        if(sq_sum_overlap > 1.0 + 1e-6) throw std::runtime_error(fmt::format("eps larger than one: {:.16f}", sq_sum_overlap));
        if(min_overlap < 0.0) throw std::runtime_error(fmt::format("min_overlap smaller than zero: {:.16f}", min_overlap));
        if(subspace_error < target_subspace_error) {
            reason = fmt::format("subspace error is low enough: {:.3e} < threshold {:.3e}", subspace_error, target_subspace_error);
            break;
        }
        if(meta.optMode == OptMode::OVERLAP and sq_sum_overlap >= 1.0 / std::sqrt(2.0)) {
            reason = fmt::format("Overlap is sufficient:  {:.16f} >= threshold {:.16f}", max_overlap, 1.0 / std::sqrt(2.0));
            break;
        }
    }
    tools::log->debug("Finished iterative eigensolver -- reason: {}", reason);
    return {eigvecs, eigvals};
}

template std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_part<cplx>(const TensorsFinite &tensors, double energy_target,
                                                                                         double target_subspace_error, const OptMeta &meta);
template std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_part<real>(const TensorsFinite &tensors, double energy_target,
                                                                                         double target_subspace_error, const OptMeta &meta);

template<typename MatrixProductType>
void subs_bfgs_preconditioner(void *x, int *ldx, void *y, int *ldy, int *blockSize, [[maybe_unused]] primme_params *primme, int *ierr) {
    if(x == nullptr) return;
    if(y == nullptr) return;
    if(primme == nullptr) return;
    using T               = typename MatrixProductType::Scalar;
    auto        H_ptr     = static_cast<MatrixProductType *>(primme->matrix);
    auto        shape_mps = H_ptr->get_shape_mps();
    const auto &mpo       = H_ptr->get_mpo();
    const auto &envL      = H_ptr->get_envL();
    const auto &envR      = H_ptr->get_envR();

    for(int i = 0; i < *blockSize; i++) {
        auto  xt                    = Eigen::TensorMap<const Eigen::Tensor<T, 3>>(static_cast<T *>(x) + *ldx * i, shape_mps);
        auto  yt                    = Eigen::TensorMap<Eigen::Tensor<T, 3>>(static_cast<T *>(y) + *ldy * i, shape_mps);
        auto  summary               = ceres::GradientProblemSolver::Summary();
        auto *functor               = new tools::finite::opt::internal::bfgs_simps_functor<T>(xt, envL, envR, mpo);
        auto  problem               = ceres::GradientProblem(functor);
        auto  options               = tools::finite::opt::internal::bfgs_default_options;
        options.max_num_iterations  = 1000;
        options.function_tolerance  = 1e-6;
        options.gradient_tolerance  = 1e-6;
        options.parameter_tolerance = 1e-6;
        if constexpr(std::is_same_v<T, double>) {
            ceres::Solve(options, problem, yt.data(), &summary);
        } else {
            auto y_cplx_as_2x_real = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double *>(yt.data()), 2 * yt.size());
            ceres::Solve(options, problem, y_cplx_as_2x_real.data(), &summary);
        }

        tools::log->info("LBFGS Preconditioner: i {} | f {:8.5e} | size {} | time {:8.5e} | it {} | exit: {}. Message: {}", i, summary.final_cost,
                         summary.num_parameters, summary.total_time_in_seconds, summary.iterations.size(),
                         ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());
    }
    *ierr = 0;
}

template<typename Scalar>
std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_prec(const TensorsFinite &tensors, double energy_target, double target_subspace_error,
                                                                          const OptMeta &meta) {
    tools::log->trace("Finding subspace -- partial");
    auto  t_iter = tid::tic_scope("part");
    auto &t_lu   = tid::get("lu_decomp");

    // Initial mps and a vector map
    const auto &multisite_mps    = tensors.state->get_multisite_mps();
    const auto  multisite_vector = Eigen::Map<const Eigen::VectorXcd>(multisite_mps.data(), multisite_mps.size());
    const auto  problem_size     = multisite_mps.size();

    // Mutable initial mps vector used for initial guess in arpack
    Eigen::Tensor<Scalar, 3> init;
    if constexpr(std::is_same_v<Scalar, real>) init = multisite_mps.real();
    if constexpr(std::is_same_v<Scalar, cplx>) init = multisite_mps;

    // Get the local effective Hamiltonian as a matrix
    //    const auto &effective_hamiltonian = tensors.get_effective_hamiltonian<Scalar>();
    //    double      time_ham              = tid::get("ham").get_last_interval();
    double time_ham = 0;

    // Create the dense matrix object for the eigenvalue solver
    //    MatVecDense<Scalar> hamiltonian(effective_hamiltonian.data(), effective_hamiltonian.dimension(0), false, eig::Form::SYMM, eig::Side::R);
    const auto       &mpo = tensors.get_multisite_mpo();
    const auto       &env = tensors.get_multisite_env_ene_blk();
    MatVecMPO<Scalar> hamiltonian(env.L, env.R, mpo);

    // Create a reusable config for multiple nev trials
    eig::settings config;
    config.tol                   = 1e-10; // settings::precision::eigs_tolerance;
    config.maxIter               = 5000;
    config.maxTime               = 60 * 60;
    config.compute_eigvecs       = eig::Vecs::ON;
    config.compress              = settings::precision::use_compressed_mpo_squared_otf;
    config.lib                   = eig::Lib::PRIMME;
    config.ritz                  = eig::Ritz::primme_closest_abs;
    config.primme_target_shifts  = {energy_target};
    config.primme_method         = eig::PrimmeMethod::PRIMME_GD_Olsen_plusK;
    config.primme_preconditioner = subs_bfgs_preconditioner<MatVecMPO<Scalar>>;
    config.primme_projection     = "primme_proj_refined";
    config.primme_locking        = true;

    //    config.maxNcv               = static_cast<eig::size_type>(16); // arpack needs ncv ~512 to handle all cases. Primme seems content with 16.
    config.loglevel = 2;
    config.initial_guess.push_back({init.data(), 0});

    std::string reason = "exhausted";

    // Initialize eigvals/eigvecs containers that store the results
    Eigen::VectorXd  eigvals;
    Eigen::MatrixXcd eigvecs;
    for(auto nev : generate_nev_list(static_cast<int>(problem_size))) {
        eig::solver solver;
        solver.config        = config;
        solver.config.maxNev = nev;
        solver.config.maxNcv = std::clamp(8 * nev, nev + 1, static_cast<int>(problem_size));

        // Set the new initial guess if we are doing one more round
        if(eigvecs.cols() != 0) {
            solver.config.initial_guess.clear();
            for(long n = 0; n < eigvecs.cols(); n++) { solver.config.initial_guess.push_back({eigvecs.col(n).data(), n}); }
        }

        solver.eigs(hamiltonian);
        t_lu += *hamiltonian.t_factorOP;

        eigvals = eig::view::get_eigvals<eig::real>(solver.result, false);
        eigvecs = eig::view::get_eigvecs<eig::cplx>(solver.result, eig::Side::R, false);

        // Check the quality of the subspace
        Eigen::VectorXd overlaps       = (multisite_vector.adjoint() * eigvecs).cwiseAbs().real();
        double          max_overlap    = overlaps.maxCoeff();
        double          min_overlap    = overlaps.minCoeff();
        double          sq_sum_overlap = overlaps.cwiseAbs2().sum();
        double          subspace_error = 1.0 - sq_sum_overlap;
        reports::subs_add_entry(nev, max_overlap, min_overlap, subspace_error, solver.result.meta.time_total, time_ham, t_lu.get_last_interval(),
                                solver.result.meta.iter, solver.result.meta.num_mv, solver.result.meta.num_pc);
        time_ham = 0;
        if(max_overlap > 1.0 + 1e-6) throw std::runtime_error(fmt::format("max_overlap larger than one: {:.16f}", max_overlap));
        if(sq_sum_overlap > 1.0 + 1e-6) throw std::runtime_error(fmt::format("eps larger than one: {:.16f}", sq_sum_overlap));
        if(min_overlap < 0.0) throw std::runtime_error(fmt::format("min_overlap smaller than zero: {:.16f}", min_overlap));
        if(subspace_error < target_subspace_error) {
            reason = fmt::format("subspace error is low enough: {:.3e} < threshold {:.3e}", subspace_error, target_subspace_error);
            break;
        }
        if(meta.optMode == OptMode::OVERLAP and sq_sum_overlap >= 1.0 / std::sqrt(2.0)) {
            reason = fmt::format("Overlap is sufficient:  {:.16f} >= threshold {:.16f}", max_overlap, 1.0 / std::sqrt(2.0));
            break;
        }
        hamiltonian.reset();
    }
    tools::log->debug("Finished iterative eigensolver -- reason: {}", reason);
    return {eigvecs, eigvals};
}

template std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_prec<cplx>(const TensorsFinite &tensors, double energy_target,
                                                                                         double target_subspace_error, const OptMeta &meta);
template std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_prec<real>(const TensorsFinite &tensors, double energy_target,
                                                                                         double target_subspace_error, const OptMeta &meta);

template<typename Scalar>
std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_full(const TensorsFinite &tensors) {
    tools::log->trace("Finding subspace -- full diag");
    auto t_full = tid::tic_scope("full");
    // Generate the Hamiltonian matrix
    const auto &effective_hamiltonian = tensors.get_effective_hamiltonian<Scalar>();

    // Create a solver and diagonalize the local effective Hamiltonian
    eig::solver solver;
    solver.eig<eig::Form::SYMM>(effective_hamiltonian.data(), effective_hamiltonian.dimension(0), eig::Vecs::ON, eig::Dephase::OFF);
    tools::log->debug("Finished eigensolver -- reason: Full diagonalization");

    const auto &multisite_mps = tensors.state->get_multisite_mps();
    const auto  multisite_vec = Eigen::Map<const Eigen::VectorXcd>(multisite_mps.data(), multisite_mps.size());

    auto            eigvals  = eig::view::get_eigvals<double>(solver.result);
    auto            eigvecs  = eig::view::get_eigvecs<Scalar>(solver.result);
    Eigen::VectorXd overlaps = (multisite_vec.adjoint() * eigvecs).cwiseAbs().real();
    int             idx;
    double          max_overlap    = overlaps.maxCoeff(&idx);
    double          min_overlap    = overlaps.minCoeff();
    double          sq_sum_overlap = overlaps.cwiseAbs2().sum();
    double          subspace_error = 1.0 - sq_sum_overlap;
    long            nev            = eigvecs.cols();
    auto            time_eig       = tid::get("eig").get_last_interval();
    auto            time_ham       = tid::get("ham").get_last_interval();
    reports::subs_add_entry(nev, max_overlap, min_overlap, subspace_error, time_eig, time_ham, 0, 1, 0, 0);
    return {eigvecs, eigvals};
}

template std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_full<cplx>(const TensorsFinite &tensors);
template std::pair<Eigen::MatrixXcd, Eigen::VectorXd> subspace::find_subspace_full<real>(const TensorsFinite &tensors);

template<typename T>
MatrixType<T> subspace::get_hamiltonian_squared_in_subspace(const ModelFinite &model, const EdgesFinite &edges, const std::vector<opt_mps> &eigvecs) {
    // First, make sure every candidate is actually a basis vector, otherwise this computation would turn difficult if we have to skip rows and columns
    auto t_ham = tid::tic_scope("ham²_sub");
    for(const auto &eigvec : eigvecs)
        if(not eigvec.is_basis_vector)
            throw std::runtime_error("One eigvec is not a basis vector. When constructing a hamiltonian subspace matrix, make sure the candidates are all "
                                     "eigenvectors/basis vectors");

    const auto &env2 = edges.get_multisite_env_var_blk();
    const auto &mpo2 = model.get_multisite_mpo_squared();

    tools::log->trace("Contracting subspace hamiltonian squared new");
    long dim0   = mpo2.dimension(2);
    long dim1   = env2.L.dimension(0);
    long dim2   = env2.R.dimension(0);
    long eignum = static_cast<long>(eigvecs.size()); // Number of eigenvectors

    Eigen::Tensor<std::complex<double>, 0> H2_ij;
    Eigen::Tensor<std::complex<double>, 3> H2_mps(dim0, dim1, dim2); // The local hamiltonian multiplied by mps at column j.
    Eigen::MatrixXcd                       H2_sub(eignum, eignum);   // The local hamiltonian projected to the subspace (spanned by eigvecs)
    for(auto col = 0; col < eignum; col++) {
        const auto &mps_j = std::next(eigvecs.begin(), col)->get_tensor();
        tools::common::contraction::matrix_vector_product(H2_mps, mps_j, mpo2, env2.L, env2.R);
        for(auto row = col; row < eignum; row++) {
            const auto &mps_i                    = std::next(eigvecs.begin(), row)->get_tensor();
            H2_ij.device(tenx::omp::getDevice()) = mps_i.conjugate().contract(H2_mps, tenx::idx({0, 1, 2}, {0, 1, 2}));
            H2_sub(row, col)                     = H2_ij(0);
            H2_sub(col, row)                     = std::conj(H2_ij(0));
        }
    }
    if constexpr(std::is_same_v<T, double>)
        return H2_sub.real();
    else
        return H2_sub;
}

// Explicit instantiations
template MatrixType<real> subspace::get_hamiltonian_squared_in_subspace(const ModelFinite &model, const EdgesFinite &edges,
                                                                        const std::vector<opt_mps> &eigvecs);
template MatrixType<cplx> subspace::get_hamiltonian_squared_in_subspace(const ModelFinite &model, const EdgesFinite &edges,
                                                                        const std::vector<opt_mps> &eigvecs);