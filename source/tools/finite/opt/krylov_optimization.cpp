
#include "../opt_meta.h"
#include "../opt_mps.h"
#include "opt-internal.h"
#include "report.h"
#include <algorithms/AlgorithmStatus.h>
#include <config/settings.h>
#include <general/iter.h>
#include <math/eig.h>
#include <math/eig/matvec/matvec_mps.h>
#include <math/tenx.h>
#include <tensors/model/ModelFinite.h>
#include <tensors/state/StateFinite.h>
#include <tensors/TensorsFinite.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>

// Temporary
//#include <h5pp/h5pp.h>

void tools::finite::opt::internal::krylov_extract_solutions(const TensorsFinite &tensors, const opt_mps &initial_mps, const eig::solver &solver,
                                                            std::vector<opt_mps> &results, const OptMeta &meta, bool converged_only) {
    auto dims_mps = initial_mps.get_tensor().dimensions();
    if(solver.result.meta.eigvals_found and solver.result.meta.eigvecsR_found) {
        auto eigvecs = eig::view::get_eigvecs<cplx>(solver.result, eig::Side::R, converged_only);
        auto eigvals = eig::view::get_eigvals<real>(solver.result, converged_only);
        if(eigvecs.cols() == eigvals.size()) {
            double overlap_sq_sum = 0;
            for(long idx = 0; idx < eigvals.size(); idx++) {
                // It's important to normalize the eigenvectors - they are not always well normalized when we get them from the eig::solver
                auto eigvec_i = tenx::TensorCast(eigvecs.col(idx).normalized(), dims_mps);
                auto overlap  = std::abs(initial_mps.get_vector().dot(eigvecs.col(idx)));
                auto energy   = tools::finite::measure::energy(eigvec_i, tensors);
                auto eigval   = energy - initial_mps.get_energy_reduced();
                auto variance = tools::finite::measure::energy_variance(eigvec_i, tensors);
                results.emplace_back(fmt::format("{:<8} eigenvector {}", solver.config.tag, idx), eigvec_i, tensors.active_sites, eigval,
                                     initial_mps.get_energy_reduced(), variance, overlap, tensors.get_length());
                auto &mps = results.back();
                mps.set_time(solver.result.meta.time_total);
                mps.set_counter(static_cast<size_t>(solver.result.meta.matvecs));
                mps.set_iter(static_cast<size_t>(solver.result.meta.iter));
                mps.set_max_grad(tools::finite::measure::max_gradient(eigvec_i, tensors));
                mps.is_basis_vector = true;
                mps.validate_basis_vector();
                mps.set_krylov_idx(idx);
                mps.set_krylov_nev(solver.result.meta.nev_converged);
                mps.set_krylov_ncv(solver.result.meta.ncv);
                mps.set_krylov_tol(solver.result.meta.tol);
                mps.set_krylov_eigval(eigvals(idx));
                mps.set_krylov_ritz(solver.result.meta.ritz);
                mps.set_krylov_shift(solver.result.meta.sigma);
                mps.set_optmode(meta.optMode);
                mps.set_optspace(meta.optSpace);
                if(not solver.result.meta.residual_norms.empty()) mps.set_krylov_resid(solver.result.meta.residual_norms.at(static_cast<size_t>(idx)));

                // When doing full diagonalization, getting all the solutions is costly, and we basically never need them all. Let's instead break after a small
                // number of solutions
                if(eigvals.size() == eigvecs.rows() /* Detects full diag*/) {
                    overlap_sq_sum += overlap * overlap;
                    if(overlap_sq_sum > 0.5) break;
                }
            }
        }
    }
}

namespace tools::finite::opt::internal {

    template<typename Scalar>
    struct opt_init_t {
        Eigen::Tensor<Scalar, 3> mps = {};
        long                     idx = 0;
    };

    // Make a handy variance comparator
    auto comp_variance = [](const opt_mps &lhs, const opt_mps &rhs) {
        if(lhs.get_krylov_idx() != rhs.get_krylov_idx()) return lhs.get_krylov_idx() < rhs.get_krylov_idx();
        return lhs.get_variance() < rhs.get_variance();
    };
    auto comp_gradient = [](const opt_mps &lhs, const opt_mps &rhs) {
        if(lhs.get_krylov_idx() != rhs.get_krylov_idx()) return lhs.get_krylov_idx() < rhs.get_krylov_idx();
        return lhs.get_max_grad() < rhs.get_max_grad();
    };
    auto comp_gradient_ref = [](const std::reference_wrapper<const opt_mps> &lhs_ref, const std::reference_wrapper<const opt_mps> &rhs_ref) {
        const auto &lhs = lhs_ref.get();
        const auto &rhs = rhs_ref.get();
        if(lhs.get_krylov_idx() != rhs.get_krylov_idx()) return lhs.get_krylov_idx() < rhs.get_krylov_idx();
        return lhs.get_max_grad() < rhs.get_max_grad();
    };

    auto comp_overlap = [](const opt_mps &lhs, const opt_mps &rhs) {
        return lhs.get_overlap() > rhs.get_overlap();
    };

    //    auto min_gradient_idx = [](const std::vector<opt_mps> &elems, long idx) {
    //        double found_grad_norm = 1e20;
    //        long   found_elem_idx  = -1;
    //        long   elem_idx        = 0;
    //        for(const auto &e : elems) {
    //            if(e.get_krylov_idx() == idx and e.get_max_grad() < found_grad_norm) {
    //                found_grad_norm = e.get_max_grad();
    //                found_elem_idx  = elem_idx;
    //            }
    //            elem_idx++;
    //        }
    //        return found_elem_idx;
    //    };

    template<typename Scalar>
    Eigen::Tensor<Scalar, 3> get_initial_guess(const opt_mps &initial_mps, const std::vector<opt_mps> &results) {
        if(results.empty()) {
            if constexpr(std::is_same_v<Scalar, real>)
                return initial_mps.get_tensor().real();
            else
                return initial_mps.get_tensor();
        } else {
            // Return whichever of initial_mps or results that has the lowest variance or gradient
            auto it = std::min_element(results.begin(), results.end(), comp_gradient);
            if(it == results.end()) return get_initial_guess<Scalar>(initial_mps, {});

            if(it->get_max_grad() < initial_mps.get_max_grad()) {
                tools::log->debug("Previous result is a good initial guess: {} | var {:8.2e}  ∇fₘₐₓ {:8.2e}", it->get_name(), it->get_variance(),
                                  it->get_max_grad());
                return get_initial_guess<Scalar>(*it, {});
            } else
                return get_initial_guess<Scalar>(initial_mps, {});
        }
    }

    template<typename Scalar>
    std::vector<opt_init_t<Scalar>> get_initial_guesses(const opt_mps &initial_mps, const std::vector<opt_mps> &results, long nev) {
        std::vector<opt_init_t<Scalar>> init;
        if(results.empty()) {
            if constexpr(std::is_same_v<Scalar, real>)
                init.push_back({initial_mps.get_tensor().real(), 0});
            else
                init.push_back({initial_mps.get_tensor(), 0});
        } else {
            for(long n = 0; n < nev; n++) {
                // Take the latest result with idx == n

                // Start by collecting the results with the correct index
                std::vector<std::reference_wrapper<const opt_mps>> results_idx_n;
                for(const auto &r : results) {
                    if(r.get_krylov_idx() == n) results_idx_n.emplace_back(r);
                }
                if(not results_idx_n.empty()) {
                    if constexpr(std::is_same_v<Scalar, real>) {
                        init.push_back({results_idx_n.back().get().get_tensor().real(), n});
                    } else {
                        init.push_back({results_idx_n.back().get().get_tensor(), n});
                    }
                }
                //
                //
                //                // Return whichever of initial_mps or results that has the lowest variance or gradient
                //
                //
                //
                //                long min_idx = min_gradient_idx(results, n);
                //                if(min_idx >= 0) {
                //                    auto &res = results.at(static_cast<size_t>(min_idx));
                //                    tools::log->debug("Found good initial guess for nev {}: idx {} {:<34} | lg var {:.8f}  ∇fₘₐₓ {:8.2e}", n, min_idx,
                //                    res.get_name(),
                //                                      std::log10(res.get_variance()), res.get_max_grad());
                //
                //                    if constexpr(std::is_same_v<Scalar, real>)
                //                        init.push_back({res.get_tensor().real(), idx});
                //                    else
                //                        init.push_back({res.get_tensor(), idx});
                //                }
            }
        }
        if(init.size() > static_cast<size_t>(nev)) throw std::logic_error("Found too many initial guesses");
        return init;
    }

    bool try_harder(const std::vector<opt_mps> &results) {
        if(results.empty()) {
            tools::log->debug("Try harder: true | first run");
            return true;
        }
        if(results.back().get_name().find("lapack") != std::string::npos) {
            tools::log->debug("Try harder: false | full diagonalization is good enough");
            return false; // Full diag. You are finished.
        }
        // It's important that we only consider the "eigenvector 0" results, i.e. idx == 0.
        std::vector<std::reference_wrapper<const opt_mps>> results_eig0;
        for(const auto &r : results) {
            if(r.get_krylov_idx() == 0) results_eig0.emplace_back(r);
        }

        const auto &back0 = results_eig0.back().get();
        if(back0.get_krylov_tol() <= 1.1 * settings::precision::eig_tolerance) { // 1.1 so to get a non-exact match as well
            tools::log->debug("Try harder: false | minimum tolerance reached: {:8.2e}", back0.get_krylov_tol());
            return false;
        }

        auto results_eig0_min_it = std::min_element(results_eig0.begin(), results_eig0.end(), comp_gradient_ref);
        if(results_eig0_min_it == results_eig0.end()) {
            tools::log->debug("Try harder: true | could not find minimum");
            return true;
        }
        const auto &result_eig0_min_grad_norm = results_eig0_min_it->get();
        bool        tryharder                 = result_eig0_min_grad_norm.get_max_grad() >= 1e0;
        tools::log->debug("Try harder: {} | grad norm {:8.2e} threshold = 1", tryharder, result_eig0_min_grad_norm.get_max_grad());
        return tryharder;
        //        if(result_eig0_min_grad_norm.get_max_grad() >= 1e0) {
        //            return true; // Very large gradient, just try harder
        //        }
        //        if(result_eig0_min_grad_norm.get_max_grad() < 1e0) {
        //            tools::log->debug("Try harder: false | grad norm {:8.2e} < 1", result_eig0_min_grad_norm.get_max_grad());
        //            return false; // We rarely make progress once grad max norm is smaller than one
        //        }
        //        tools::log->debug("Try harder: {} | last resort: check if lowest grad norm {:8.2e} > 1", result_eig0_min_grad_norm.get_max_grad() > 1e-0,
        //                          result_eig0_min_grad_norm.get_max_grad());
        //        return result_eig0_min_grad_norm.get_max_grad() > 1e-0; // We rarely improve variance once grad max norm is smaller than one
    }

    template<typename Scalar>
    double get_largest_eigenvalue_hamiltonian_squared(const TensorsFinite &tensors) {
        const auto &env2                = tensors.get_multisite_env_var_blk();
        auto        hamiltonian_squared = MatVecMps<Scalar>(env2.L, env2.R, tensors.get_multisite_mpo_squared());
        tools::log->trace("Finding largest-magnitude eigenvalue");
        eig::solver solver; // Define a solver just to find the maximum eigenvalue
        solver.config.tol             = settings::precision::eig_tolerance;
        solver.config.maxIter         = 200;
        solver.config.maxNev          = 1;
        solver.config.maxNcv          = 16;
        solver.config.compute_eigvecs = eig::Vecs::OFF;
        solver.config.sigma           = std::nullopt;
        solver.config.compress        = settings::precision::use_compressed_mpo_squared_otf;
        solver.config.ritz            = eig::Ritz::LM;
        solver.setLogLevel(2);
        solver.eigs(hamiltonian_squared);
        return eig::view::get_eigvals<double>(solver.result).cwiseAbs().maxCoeff();
    }

    template<typename Scalar>
    void krylov_executor(eig::solver &solver, MatVecMps<Scalar> &hamiltonian_squared, const TensorsFinite &tensors, const opt_mps &initial_mps,
                         std::vector<opt_mps> &results, const OptMeta &meta) {
        if(std::is_same_v<Scalar, cplx> and meta.optType == OptType::REAL)
            throw std::logic_error("krylov_launcher error: Mixed Scalar:cplx with OptType::REAL");
        if(std::is_same_v<Scalar, real> and meta.optType == OptType::CPLX)
            throw std::logic_error("krylov_launcher error: Mixed Scalar:real with OptType::CPLX");

        const auto       &mpo = tensors.get_multisite_mpo();
        const auto       &env = tensors.get_multisite_env_ene_blk();
        MatVecMps<Scalar> hamiltonian(env.L, env.R, mpo);
        solver.config.primme_extra = &hamiltonian;

//        h5pp::File  h5file("../output/primme_mps.h5", h5pp::FilePermission::READWRITE);
//        long        number = 0;
//        std::string groupname;
//        while(true) {
//            groupname = fmt::format("mps-{}", number++);
//            if(not h5file.linkExists(groupname)) {
//                h5file.writeDataset(hamiltonian.get_mpo(), groupname + "/mpo", H5D_CHUNKED);
//                h5file.writeDataset(hamiltonian.get_envL(), groupname + "/envL", H5D_CHUNKED);
//                h5file.writeDataset(hamiltonian.get_envR(), groupname + "/envR", H5D_CHUNKED);
//                h5file.writeDataset(hamiltonian_squared.get_mpo(), groupname + "/mpo2", H5D_CHUNKED);
//                h5file.writeDataset(hamiltonian_squared.get_envL(), groupname + "/envL2", H5D_CHUNKED);
//                h5file.writeDataset(hamiltonian_squared.get_envR(), groupname + "/envR2", H5D_CHUNKED);
//                auto init = get_initial_guesses<Scalar>(initial_mps, results, solver.config.maxNev.value()); // Init holds the data in memory for this scope
//                for(auto &i : init) h5file.writeDataset(i.mps, fmt::format("{}/mps_init_{}", groupname, i.idx), H5D_CHUNKED);
//                h5file.writeAttribute(tensors.active_problem_dims(), "dimensions", groupname);
//                break;
//            }
//        }

        hamiltonian_squared.reset();
        auto size = tensors.active_problem_size();

        if(size <= 1024) {
            tools::log->trace("Full diagonalization of (H-E)²");
            auto matrix       = tools::finite::opt::internal::get_multisite_hamiltonian_squared_matrix<double>(*tensors.model, *tensors.edges);
            solver.config.tag = "lapack";
            solver.eig(matrix.data(), matrix.rows());
        } else {
            auto init = get_initial_guesses<Scalar>(initial_mps, results, solver.config.maxNev.value()); // Init holds the data in memory for this scope
            for(auto &i : init) solver.config.initial_guess.push_back({i.mps.data(), i.idx});
            tools::log->trace("Defining reduced Hamiltonian-squared matrix-vector product");

            if(tensors.model->is_compressed_mpo_squared()) {
                tools::log->warn("Finding excited state using H² with ritz SM because all MPO²'s are compressed!");
                solver.config.ritz = eig::Ritz::SM;
                solver.eigs(hamiltonian_squared);
            } else {
                if(solver.config.lib == eig::Lib::ARPACK) {
                    if(not solver.config.ritz) solver.config.ritz = eig::Ritz::LM;
                    if(not solver.config.sigma) {
                        double largest_eigval = get_largest_eigenvalue_hamiltonian_squared<Scalar>(tensors) + 1.0; // Add one just to make sure we shift enough
                        solver.config.sigma   = largest_eigval;
                    }
                    tools::log->debug("Finding excited state using shifted operator [(H-E)²-σ], with σ = {:.16f} | arpack {} | init on ...",
                                      std::real(solver.config.sigma.value()), eig::RitzToString(solver.config.ritz.value()));
                    solver.eigs(hamiltonian_squared);
                } else if(solver.config.lib == eig::Lib::PRIMME) {
                    if(not solver.config.ritz) solver.config.ritz = eig::Ritz::SA;
                    if(not solver.config.sigma) solver.config.sigma = 1.0;
                    tools::log->debug("Finding excited state using shifted operator [(H-E)²-σ], with σ = {:.16f} | primme {} | init on ...",
                                      std::real(solver.config.sigma.value()), eig::RitzToString(solver.config.ritz.value()));
                    solver.eigs(hamiltonian_squared);
                }
            }
        }
        krylov_extract_solutions(tensors, initial_mps, solver, results, meta, false);
    }

    template<typename Scalar>
    void krylov_manager(const TensorsFinite &tensors, const opt_mps &initial_mps, std::vector<opt_mps> &results, const OptMeta &meta) {
        eig::settings config_primme;
        config_primme.tol             = settings::precision::eig_tolerance;
        config_primme.maxIter         = 80000;
        config_primme.maxTime         = 60 * 60;
        config_primme.maxNev          = 1;
        config_primme.maxNcv          = 16; // 128
        config_primme.compress        = settings::precision::use_compressed_mpo_squared_otf;
        config_primme.lib             = eig::Lib::PRIMME;
        config_primme.ritz            = eig::Ritz::SA;
        config_primme.compute_eigvecs = eig::Vecs::ON;
        config_primme.tag             = "primme";
        config_primme.primme_method   = eig::PrimmeMethod::PRIMME_GD_Olsen_plusK;
        config_primme.loglevel        = 2;
        std::vector<eig::settings> configs(2);

        config_primme.primme_grad_tol  = 1e0;
        config_primme.primme_grad_iter = 100;
        config_primme.primme_grad_time = 5;

        configs[0]                             = config_primme;
        configs[0].tol                         = 1e3 * settings::precision::eig_tolerance;
        configs[0].maxNcv                      = 16;
        configs[0].primme_max_inner_iterations = -1;

        configs[1]                             = config_primme;
        configs[1].tol                         = 1e0 * settings::precision::eig_tolerance;
        configs[1].maxNcv                      = 32;
        configs[1].primme_max_inner_iterations = -1;

        const auto                      &env2 = tensors.get_multisite_env_var_blk();
        std::optional<MatVecMps<Scalar>> hamiltonian_squared;
        for(const auto &config : configs) {
            eig::solver solver;
            solver.config = config;
            if(not hamiltonian_squared) hamiltonian_squared = MatVecMps<Scalar>(env2.L, env2.R, tensors.get_multisite_mpo_squared());
            krylov_executor<Scalar>(solver, hamiltonian_squared.value(), tensors, initial_mps, results, meta);

            if(not try_harder(results)) break;
            //            if(solver.config.lib == eig::Lib::ARPACK) hamiltonian_squared.reset();
        }
    }

    template<typename Scalar>
    void krylov_manager2(const TensorsFinite &tensors, const opt_mps &initial_mps, std::vector<opt_mps> &results, const OptMeta &meta) {
        eig::settings config_primme;
        config_primme.tol             = settings::precision::eig_tolerance;
        config_primme.maxIter         = 80000;
        config_primme.maxTime         = 60 * 60;
        config_primme.maxNev          = 1;
        config_primme.maxNcv          = 16; // 128
        config_primme.compress        = settings::precision::use_compressed_mpo_squared_otf;
        config_primme.lib             = eig::Lib::PRIMME;
        config_primme.ritz            = eig::Ritz::SA;
        config_primme.compute_eigvecs = eig::Vecs::ON;
        config_primme.tag             = "primme";
        config_primme.primme_method   = eig::PrimmeMethod::PRIMME_GD_Olsen_plusK;
        config_primme.loglevel        = 2;
        std::vector<eig::settings> configs(2);

        config_primme.primme_grad_tol  = 1e0;
        config_primme.primme_grad_iter = 0;
        config_primme.primme_grad_time = 0;

        configs[0]                             = config_primme;
        configs[0].tol                         = 1e3 * settings::precision::eig_tolerance;
        configs[0].maxNcv                      = 16;
        configs[0].primme_max_inner_iterations = 100;

        configs[1]                             = config_primme;
        configs[1].tol                         = 1e0 * settings::precision::eig_tolerance;
        configs[1].maxNcv                      = 32;
        configs[1].primme_max_inner_iterations = 100;

        const auto                      &env2 = tensors.get_multisite_env_var_blk();
        std::optional<MatVecMps<Scalar>> hamiltonian_squared;
        for(const auto &config : configs) {
            eig::solver solver;
            solver.config = config;
            if(not hamiltonian_squared) hamiltonian_squared = MatVecMps<Scalar>(env2.L, env2.R, tensors.get_multisite_mpo_squared());
            krylov_executor<Scalar>(solver, hamiltonian_squared.value(), tensors, initial_mps, results, meta);

            if(not try_harder(results)) break;
            if(solver.config.lib == eig::Lib::ARPACK) hamiltonian_squared.reset();
        }
    }

    template<typename Scalar>
    void krylov_manager3(const TensorsFinite &tensors, const opt_mps &initial_mps, std::vector<opt_mps> &results, const OptMeta &meta) {
        eig::settings config_primme;
        config_primme.tol             = settings::precision::eig_tolerance;
        config_primme.maxIter         = 80000;
        config_primme.maxTime         = 60 * 60;
        config_primme.maxNev          = 1;
        config_primme.maxNcv          = 16; // 128
        config_primme.compress        = settings::precision::use_compressed_mpo_squared_otf;
        config_primme.lib             = eig::Lib::PRIMME;
        config_primme.ritz            = eig::Ritz::SA;
        config_primme.compute_eigvecs = eig::Vecs::ON;
        config_primme.tag             = "primme";
        config_primme.primme_method   = eig::PrimmeMethod::PRIMME_GD_Olsen_plusK;
        config_primme.loglevel        = 2;
        std::vector<eig::settings> configs(2);

        config_primme.primme_grad_tol          = 1e0;
        config_primme.primme_grad_iter         = 10;
        config_primme.primme_grad_time         = 1;
        configs[0]                             = config_primme;
        configs[0].tol                         = 1e3 * settings::precision::eig_tolerance;
        configs[0].maxNcv                      = 16;
        configs[0].primme_max_inner_iterations = -1;

        configs[1]                             = config_primme;
        configs[1].tol                         = 1e0 * settings::precision::eig_tolerance;
        configs[1].maxNcv                      = 32;
        configs[1].primme_max_inner_iterations = -1;

        const auto                      &env2 = tensors.get_multisite_env_var_blk();
        std::optional<MatVecMps<Scalar>> hamiltonian_squared;
        for(const auto &config : configs) {
            eig::solver solver;
            solver.config = config;
            if(not hamiltonian_squared) hamiltonian_squared = MatVecMps<Scalar>(env2.L, env2.R, tensors.get_multisite_mpo_squared());
            krylov_executor<Scalar>(solver, hamiltonian_squared.value(), tensors, initial_mps, results, meta);
            if(not try_harder(results)) break;
            if(solver.config.lib == eig::Lib::ARPACK) hamiltonian_squared.reset();
        }
    }
}

tools::finite::opt::opt_mps tools::finite::opt::internal::krylov_optimization(const TensorsFinite &tensors, const opt_mps &initial_mps,
                                                                              [[maybe_unused]] const AlgorithmStatus &status, OptMeta &meta) {
    using namespace internal;
    using namespace settings::precision;
    initial_mps.validate_basis_vector();
    if(not tensors.model->is_reduced()) throw std::runtime_error("krylov_optimization requires energy-reduced MPO²");
    if(tensors.model->is_compressed_mpo_squared()) throw std::runtime_error("krylov_optimization requires non-compressed MPO²");

    auto t_var    = tid::tic_scope("krylov");
    auto dims_mps = initial_mps.get_tensor().dimensions();
    auto size     = initial_mps.get_tensor().size();

    tools::log->debug("Finding excited state: minimum eigenstate of (H-E)² | dims {} = {}", dims_mps, size);
    std::vector<opt_mps> results;
    auto                 t_var1 = tid::tic_scope("krylov1");
    switch(meta.optType) {
        case OptType::REAL: krylov_manager<real>(tensors, initial_mps, results, meta); break;
        case OptType::CPLX: krylov_manager<cplx>(tensors, initial_mps, results, meta); break;
    }
    t_var1.toc();
    //    tools::log->debug("Krylov optimization time 1: {:.1f}", 1000 * t_var1->get_last_interval());
    //    auto                 t_var2 = tid::tic_scope("krylov2");
    //    std::vector<opt_mps> results2;
    //    switch(meta.optType) {
    //        case OptType::REAL: krylov_manager2<real>(tensors, initial_mps, results2, meta); break;
    //        case OptType::CPLX: krylov_manager2<cplx>(tensors, initial_mps, results2, meta); break;
    //    }
    //    t_var2.toc();
    //    tools::log->debug("Krylov optimization time 2: {:.1f}", 1000 * t_var2->get_last_interval());
    //    auto                 t_var3 = tid::tic_scope("krylov3");
    //    std::vector<opt_mps> results3;
    //    switch(meta.optType) {
    //        case OptType::REAL: krylov_manager3<real>(tensors, initial_mps, results3, meta); break;
    //        case OptType::CPLX: krylov_manager3<cplx>(tensors, initial_mps, results3, meta); break;
    //    }
    //    t_var3.toc();

    tools::log->debug("Krylov 1 optimization time: {:.1f}", 1000 * t_var1->get_last_interval());
    //    tools::log->debug("Krylov 2 optimization time: {:.1f}", 1000 * t_var2->get_last_interval());
    //    tools::log->debug("Krylov 3 optimization time: {:.1f}", 1000 * t_var3->get_last_interval());

    if(results.empty()) {
        meta.optExit = OptExit::FAIL_ERROR;
        return initial_mps; // The solver failed
    }
    bool fulldiag = results.back().get_name().find("lapack") != std::string::npos;
    for(const auto &[num, mps] : iter::enumerate(results)) {
        if(fulldiag and num >= 8) break;
        reports::krylov_add_entry(mps);
    }
    //    for(const auto &[num, mps] : iter::enumerate(results2)) {
    //        if(fulldiag and num >= 8) break;
    //        reports::krylov_add_entry(mps);
    //    }
    //    for(const auto &[num, mps] : iter::enumerate(results3)) {
    //        if(fulldiag and num >= 8) break;
    //        reports::krylov_add_entry(mps);
    //    }

    if(results.size() >= 2) {
#pragma message "Sorting according to gradient"
        if(results.back().get_name().find("lapack") != std::string::npos)
            std::sort(results.begin(), results.end(), comp_variance);
        else if(meta.optMode == OptMode::VARIANCE)
            std::sort(results.begin(), results.end(), comp_gradient);
        else if(meta.optMode == OptMode::OVERLAP)
            std::sort(results.begin(), results.end(), comp_overlap);
    }

    return results.front();
}
