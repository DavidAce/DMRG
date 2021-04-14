
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <general/nmspc_iter.h>
#include <general/nmspc_tensor_extra.h>
#include <math/eig.h>
#include <math/eig/arpack_solver/matrix_product_hamiltonian.h>
#include <math/linalg/tensor.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/contraction.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt-internal/opt-internal.h>
#include <tools/finite/opt-internal/report.h>
#include <tools/finite/opt_mps.h>

void tools::finite::opt::extract_solutions(const opt_mps &initial_mps,const class_tensors_finite &tensors, eig::solver &solver,
                                           std::vector<tools::finite::opt::opt_mps> &eigvecs_mps,  const std::string & tag ) {
    if(solver.result.meta.eigvals_found){
        auto dims_mps = initial_mps.get_tensor().dimensions();
        auto eigvecs = eig::view::get_eigvecs<Scalar>(solver.result,eig::Side::R, false);
        auto eigvals = eig::view::get_eigvals<double>(solver.result,false);
        eigvecs_mps.reserve(eigvecs_mps.size() + static_cast<size_t>(eigvals.size()));
        auto logLevel = tools::log->level();
        tools::Logger::setLogLevel(tools::log,2);
        for(long idx = 0; idx < eigvals.size(); idx++) {
            auto eigvec_i = Textra::TensorCast(eigvecs.col(idx), dims_mps);
            auto overlap  = std::abs(initial_mps.get_vector().dot(eigvecs.col(idx)));
            auto energy   = tools::finite::measure::energy(eigvec_i, tensors);
            auto eigval   = energy - initial_mps.get_energy_reduced();
            auto variance = tools::finite::measure::energy_variance(eigvec_i, tensors);
            eigvecs_mps.emplace_back(fmt::format("eigenvector {} {:6<}", idx,tag), eigvec_i, tensors.active_sites, eigval, initial_mps.get_energy_reduced(), variance,
                                     overlap, tensors.get_length());

            eigvecs_mps.back().set_time(tools::common::profile::get_default_prof()["t_eig"]->get_last_interval());
            eigvecs_mps.back().set_counter(static_cast<size_t>(solver.result.meta.counter));
            eigvecs_mps.back().set_iter(static_cast<size_t>(solver.result.meta.iter));
            eigvecs_mps.back().is_basis_vector = true;
            eigvecs_mps.back().validate_candidate();
        }
        tools::Logger::setLogLevel(tools::log,logLevel);
    }
}


tools::finite::opt::opt_mps tools::finite::opt::internal::arpack_variance_optimization(const class_tensors_finite &tensors, const opt_mps &initial_mps,
                                                                                       const class_algorithm_status &status, OptType optType, OptMode optMode,
                                                                                       OptSpace optSpace) {
    using namespace internal;
    using namespace settings::precision;
    auto t_eig  = tools::common::profile::get_default_prof()["t_eig"]->tic_token();
    auto dims_mps = initial_mps.get_tensor().dimensions();
    auto size     = initial_mps.get_tensor().size();
//    auto nev      = std::clamp<eig::size_type>(size / 8, std::min(size, 2l), 4l);
    auto nev      = 2l;
//    auto ncv_max  = static_cast<long>(std::pow(size,0.4));
    auto ncv      = std::min(size, std::max(2l*nev + 1l, 32l));

    tools::log->info("Excited state variance optimization with arpack | dims {} = {}", dims_mps, size);

    tools::log->trace("Defining eigenvalue solver");
    // NOTE:
    // We can't use the MPO-shift trick here!
    // The lower-left corner of the MPO is NOT where the energy is subtracted on a double-layer MPO,
    // especially when using MPO compression.

    eig::solver solver_full, solver_part,solver_extr, solver_test, solver_test2;
    solver_extr.config.eigThreshold = 1e-12;
    solver_extr.config.eigMaxIter   = 1000;
    solver_extr.config.eigMaxNev    = 1;
    solver_extr.config.eigMaxNcv    = 16;
    solver_extr.setLogLevel(0);

    solver_part.config.eigThreshold = 1e-5;
    solver_part.config.eigMaxIter   = 200000;
    solver_part.config.eigMaxNev    = 1;
    solver_part.config.eigMaxNcv    = 4;
    solver_part.setLogLevel(0);

    solver_test.config.eigThreshold = 1e-5;
    solver_test.config.eigMaxIter   = 1000;
    solver_test.config.eigMaxNev    = 2;
    solver_test.config.eigMaxNcv    = 32;
    solver_test.setLogLevel(0);

    solver_test2.config.eigThreshold = 1e-5;
    solver_test2.config.eigMaxIter   = 1000;
    solver_test2.config.eigMaxNev    = 2;
    solver_test2.config.eigMaxNcv    = 32;
    solver_test2.setLogLevel(0);



    if(optType == OptType::REAL) {
        tools::log->trace("- Generating real-valued multisite components");
        const auto &             env2      = tensors.get_multisite_var_blk();
        Eigen::Tensor<double, 3> env2L     = env2.L.real();
        Eigen::Tensor<double, 3> env2R     = env2.R.real();
        Eigen::Tensor<double, 4> mpo2      = tensors.get_multisite_mpo_squared().real();
        auto                     dims_mpo2 = mpo2.dimensions();
        Eigen::Tensor<double, 3> residual = initial_mps.get_tensor().real();


//        {
//            MatrixProductHamiltonian<double> matrix_test2(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
//
//
//            double eigval_target2 = initial_mps.get_eigval()*initial_mps.get_eigval();
//            tools::log->trace("Finding excited state with squared inverse | residual on ... ");
//            solver_test2.eigs(matrix_test2, 1, 16, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, eigval_target2, eig::Shinv::ON, eig::Vecs::ON, eig::Dephase::OFF, residual.data());
//            tools::log->trace("Finding excited state with squared inverse | residual on ... OK: iters {} | counter {} | time {:.3f} | matprod {:.3f}s",
//                              solver_test2.result.meta.iter, solver_test2.result.meta.counter, solver_test2.result.meta.time_total,solver_test2.result.meta.time_matprod);
//
//        }



        if(size <= settings::precision::max_size_full_diag){
//        if(size <= 256){
            tools::log->trace("Finding full spectrum");
            auto matrix = tools::finite::opt::internal::get_multisite_hamiltonian_squared_matrix<double>(*tensors.model,*tensors.edges);
            solver_full.eig(matrix.data(), matrix.rows());
        }else{
            tools::log->trace("Defining reduced Hamiltonian-squared matrix-vector product");
            MatrixProductHamiltonian<double> matrix_extr(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
            MatrixProductHamiltonian<double> matrix_part(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);

            if(not tensors.model->is_compressed_mpo_squared()){

                tools::log->trace("Finding largest-magnitude state");
                solver_extr.eigs(matrix_extr, 1, 16, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::OFF, eig::Dephase::OFF);
                auto eigvals_extr = eig::view::get_eigvals<double>(solver_extr.result);
                auto shift = eigvals_extr.cwiseAbs().maxCoeff() + 1.0; // Add one just to make sure we shift enough


                tools::log->trace("Finding excited state with shift {:.16f} | residual on ...",shift );
                solver_part.eigs(matrix_part, nev, ncv, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, shift, eig::Shinv::OFF, eig::Vecs::ON, eig::Dephase::OFF, residual.data());
                tools::log->trace("Finding excited state with shift {:.16f} | residual on ... OK: iters {} | counter {} | time {:.3f} | matprod {:.3f}s",
                                  shift, solver_part.result.meta.iter, solver_part.result.meta.counter,solver_part.result.meta.time_total,solver_part.result.meta.time_matprod);

            }else{
                tools::log->trace("Finding excited state");
                solver_part.eigs(matrix_part, nev, ncv, eig::Ritz::SM, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::ON, eig::Dephase::OFF, residual.data());
            }
        }

//        {
//            tools::log->trace("Testing matrix-free solver");
//            MatrixReplacement<double> matRepl;
//            matRepl.attachTensors(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo);
//
//            Eigen::VectorXd b(size), x;
//            b.setZero();
//
//            // Solve Ax = b using various iterative solver with matrix-free version:
//            {
//                Eigen::ConjugateGradient<MatrixReplacement<double>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg;
//                cg.setMaxIterations(500);
//                cg.compute(matRepl);
//                x = cg.solve(b);
//                std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
//            }
//        }

    }
    if(optType == OptType::CPLX) {
        tools::log->trace("- Generating cplx-valued multisite components");
        const auto &env2     = tensors.get_multisite_var_blk();
        const auto &env2L    = env2.L;
        const auto &env2R    = env2.R;
        const auto &mpo2     = tensors.get_multisite_mpo_squared();
        const auto  dims_mpo = mpo2.dimensions();
        tools::log->trace("Defining reduced Hamiltonian-squared matrix-vector product");
        MatrixProductHamiltonian<Scalar> matrix(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo);

        // Since we use reduced energy mpo's, we set a sigma shift = 1.0 to move the smallest eigenvalue away from 0, which
        // would otherwise cause trouble for Arpack. This equates to subtracting sigma * identity from the bottom corner of the mpo.
        // The resulting eigenvalue will be shifted by the same amount, but the eigenvector will be the same, and that's what we keep.
        Eigen::Tensor<Scalar, 3> residual = initial_mps.get_tensor();
        tools::log->trace("Finding excited state");
        solver_part.eigs(matrix, nev, ncv,  eig::Ritz::SR, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::ON, eig::Dephase::OFF,residual.data());
    }
    t_eig.toc();

    std::vector<opt_mps> eigvecs_mps;
    extract_solutions(initial_mps, tensors, solver_part, eigvecs_mps, "part");
    extract_solutions(initial_mps, tensors, solver_full, eigvecs_mps, "full");
    extract_solutions(initial_mps, tensors, solver_test, eigvecs_mps, "bic1");
    extract_solutions(initial_mps, tensors, solver_test2, eigvecs_mps,"bic2");
    if (eigvecs_mps.empty()) return internal::ceres_direct_optimization(tensors, initial_mps, status, optType, OptMode::VARIANCE, OptSpace::DIRECT); // The solver failed


    auto comp_variance =  [](const opt_mps &lhs, const opt_mps &rhs){return lhs.get_variance() < rhs.get_variance();};
    auto comp_overlap  =  [](const opt_mps &lhs, const opt_mps &rhs){return lhs.get_overlap() > rhs.get_overlap();};

    if(eigvecs_mps.size() >= 2){
        if(optMode == OptMode::VARIANCE)
            std::sort(eigvecs_mps.begin(), eigvecs_mps.end(), comp_variance);
        if(optMode == OptMode::OVERLAP)
            std::sort(eigvecs_mps.begin(), eigvecs_mps.end(), comp_overlap);
    }


    for(const auto & [i,mps] : iter::enumerate(eigvecs_mps)){
        tools::log->info("{:<24} overlap {:.16f} | E {:>20.16f} | E/L {:>20.16f} | eigval {:>20.16f} {:>20.16f} | Er {:>20.16f} | variance: "
                         "{:>20.16f} | norm {:.16f} | iters {} | counter {}",
                         mps.get_name(), mps.get_overlap(), mps.get_energy(), mps.get_energy_per_site(),
                         mps.get_eigval(), mps.get_eigval(), mps.get_energy_reduced(), std::log10(mps.get_variance()),
                         mps.get_norm(),mps.get_iter(), mps.get_counter());
        if(i > 32) break;
    }

    // Try optimizing the results with l-bfgs
    std::vector<opt_mps> optimized_mps;

    // Try optimizing the first eigenvectors
    for(const auto & mps : eigvecs_mps){
        if(&mps != &eigvecs_mps.front()){
            // Break secondary states that are too far away in energy from the first one
            if(std::abs(mps.get_energy() - eigvecs_mps.front().get_energy()) > 1e2 * eigvecs_mps.front().get_variance()) break;
            // Break secondary states that are too far away in overlap from the first one
            if(mps.get_overlap() < 0.7) break;

        }
        optimized_mps.emplace_back(internal::ceres_direct_optimization(tensors, mps, status, optType, optMode, optSpace));
    }

    // If there are multiple solutions, try building a linear combination and optimize that super_mps with l-bfgs
//    if(eigvecs_mps.size() >= 2){
//        std::sort(eigvecs_mps.begin(), eigvecs_mps.end(), comp_variance);
//        auto super_vector = (eigvecs_mps[0].get_vector() + eigvecs_mps[1].get_vector()).normalized();
//        auto super_tensor = Textra::TensorCast(super_vector, dims_mps);
//        auto overlap  = std::abs(initial_mps.get_vector().dot(super_vector));
//        auto energy   = tools::finite::measure::energy(super_tensor, tensors);
//        auto eigval   = energy - initial_mps.get_energy_reduced();
//        auto variance = tools::finite::measure::energy_variance(super_tensor, tensors);
//        auto super_mps = opt_mps(fmt::format("superposition 0+1"), super_tensor, tensors.active_sites, eigval, initial_mps.get_energy_reduced(), variance,
//                           overlap, tensors.get_length());
//        optimized_mps.emplace_back(internal::ceres_direct_optimization(tensors, super_mps, status, optType, optMode, optSpace));
//
//    }

    if(optimized_mps.size() >= 2) // Sort the results in order of increasing variance (again!)
        std::sort(optimized_mps.begin(), optimized_mps.end(), comp_variance);

    if(optimized_mps.empty()) return internal::ceres_direct_optimization(tensors, initial_mps, status, optType, OptMode::VARIANCE, OptSpace::DIRECT); // The optimization strategy failed
    else return optimized_mps.front();


}
