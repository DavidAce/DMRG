//
// Created by david on 2019-07-15.
//

#include "ceres_subspace.h"
#include <Eigen/Core>
#include <iostream>
#include <simulation/class_simulation_status.h>
#include <math/class_eigsolver.h>
#include <math/arpack_extra/matrix_product_stl.h>
#include <math/arpack_extra/matrix_product_sparse.h>
#include <spdlog/spdlog.h>
#include <tools/finite/opt.h>
#include <state/class_finite_state.h>
#include <simulation/nmspc_settings.h>
#include <general/nmspc_random_numbers.h>
#include <spdlog/fmt/bundled/ranges.h>
#include <ceres/ceres.h>

using namespace tools::finite::opt;
using namespace tools::finite::opt::internals;

template<typename T> using MatrixType = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;


std::vector<int> tools::finite::opt::internals::generate_size_list(size_t shape){
    int max_nev ;
    if      (shape <= 512)  {max_nev = shape/4;}
    else if (shape > 512  and shape <= 1024) {max_nev = 128;}
    else if (shape > 1024 and shape <= 2048) {max_nev = 64;}
    else if (shape > 2048 and shape <= 4096) {max_nev = 32;}
    else if (shape > 4096 and shape <= 8192) {max_nev = 16;}
    else                                     {max_nev = 8;}

    int min_nev = std::min(std::min(8,(int)shape),max_nev);

    std::vector<int> nev_list = {min_nev};
    int tmp_nev = min_nev;
    while (tmp_nev < max_nev){
        tmp_nev = std::min(4*tmp_nev, max_nev);
        nev_list.push_back(tmp_nev);
    }
    return nev_list;
}


std::tuple<Eigen::MatrixXcd,Eigen::VectorXd,double> filter_states(const Eigen::MatrixXcd &eigvecs,const Eigen::VectorXd& eigvals, Eigen::VectorXd &overlaps, double quality_threshold, size_t max_accept){

    size_t min_accept = std::min(8ul,(size_t)eigvals.size());
    max_accept        = std::min(max_accept,(size_t)eigvals.size());
    if(min_accept == max_accept) return std::make_tuple(eigvecs,eigvals,1.0 - overlaps.cwiseAbs2().sum());
    tools::log->debug("Filtering states keeping between {} to {}, log10 quality threshold {}", min_accept,max_accept, std::log10(quality_threshold));
    double subspace_quality  = 1.0;
    Eigen::VectorXd overlaps_filtered = overlaps;
    std::vector<int>    overlaps_accepted_idx;
    std::vector<double> overlaps_accepted;

    while(true){
        int idx;
        double overlap = overlaps_filtered.maxCoeff(&idx);
        overlaps_accepted_idx.push_back(idx);
        overlaps_accepted    .push_back(overlap);
        Eigen::Map<Eigen::VectorXd> overlaps_map(overlaps_accepted.data(),overlaps_accepted.size());
        subspace_quality  = 1.0 - overlaps_map.cwiseAbs2().sum();
        if(overlaps_accepted.size() >= min_accept){
            if(subspace_quality < quality_threshold) break;
            if(overlaps_accepted.size() >= max_accept) break;
        }
        overlaps_filtered(idx) = 0;
        if(overlaps_filtered.sum() == 0) break;
    }

    Eigen::MatrixXcd eigvecs_filtered(eigvecs.rows(),overlaps_accepted.size());
    Eigen::VectorXd  eigvals_filtered(overlaps_accepted.size());

    int col_num = 0;
    for (auto &idx : overlaps_accepted_idx){
        eigvecs_filtered.col(col_num) = eigvecs.col(idx);
        eigvals_filtered    (col_num) = eigvals(idx);
        col_num++;
    }
    tools::log->debug("Filtered from {} down to {} states", eigvals.size(),eigvals_filtered.size());
    tools::log->debug("Filtered quality: log10(1-eps) = {}", std::log10(subspace_quality));
    return std::make_tuple(eigvecs_filtered,eigvals_filtered, subspace_quality);
}

int idx_best_overlap_in_window(const Eigen::VectorXd &overlaps, const Eigen::VectorXd & eigvals, double lbound, double ubound){
    assert(overlaps.size() == eigvals.size() and "idx_best_overlap_in_window: Mismatch in overlaps and eigvals sizes");
    Eigen::VectorXd overlaps_in_window = overlaps;
    for (int i = 0; i < overlaps.size(); i++){
        if (eigvals(i) > ubound) overlaps_in_window(i) = 0.0;
        if (eigvals(i) < lbound) overlaps_in_window(i) = 0.0;
    }
    if (overlaps_in_window.isZero(0.0)){overlaps_in_window = overlaps;}

    int idx;
    [[maybe_unused]] double max_overlap = overlaps_in_window.maxCoeff(&idx);
    return idx;

}

template<typename Scalar>
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
find_subspace_full(const MatrixType<Scalar> & H_local, Eigen::Tensor<std::complex<double>,3> &theta, std::vector<reports::eig_tuple> &eig_log){
    tools::log->trace("Finding subspace -- full");
    using namespace eigutils::eigSetting;
    t_eig->tic();

    Eigen::VectorXd   eigvals;
    Eigen::MatrixXcd  eigvecs;
    class_eigsolver solver;

    if constexpr (!std::is_same<Scalar, double>::value)
    {
        solver.eig<Type::CPLX, Form::SYMMETRIC>(H_local, true, false);
        eigvals = Eigen::Map<const Eigen::VectorXd>(solver.solution.get_eigvals<Form::SYMMETRIC>().data(),
                                                    solver.solution.meta.cols);
        eigvecs = Eigen::Map<const Eigen::MatrixXcd>(solver.solution.get_eigvecs<Type::CPLX, Form::SYMMETRIC>().data(),
                                                     solver.solution.meta.rows, solver.solution.meta.cols);
    }
    else
    {
        solver.eig<Type::REAL, Form::SYMMETRIC>(H_local, true, false);
        eigvals = Eigen::Map<const Eigen::VectorXd>(solver.solution.get_eigvals<Form::SYMMETRIC>().data(),
                                                    solver.solution.meta.cols);
        eigvecs = Eigen::Map<const Eigen::MatrixXd>(solver.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),
                                                    solver.solution.meta.rows, solver.solution.meta.cols);
    }
    t_eig->toc();
    tools::log->debug("Finished eigensolver -- condition: Full diagonalization");
    Eigen::Map<const Eigen::VectorXcd> theta_vec   (theta.data(),theta.size());
    Eigen::VectorXd overlaps = (theta_vec.adjoint() * eigvecs).cwiseAbs().real();
    int idx;
    double max_overlap       = overlaps.maxCoeff(&idx);
    double min_overlap       = overlaps.minCoeff();
    double sq_sum_overlap    = overlaps.cwiseAbs2().sum();
    double subspace_quality  = 1.0 - sq_sum_overlap;
    int nev = eigvecs.cols();
    eig_log.emplace_back(nev, max_overlap, min_overlap, sq_sum_overlap, std::log10(subspace_quality), t_eig->get_last_time_interval(), 0);

    return std::make_tuple(eigvecs,eigvals);
}



template<typename Scalar>
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
find_subspace_part(const MatrixType<Scalar> & H_local, Eigen::Tensor<std::complex<double>,3> &theta, double energy_target, std::vector<reports::eig_tuple> &eig_log,OptMode optMode){
    using namespace eigutils::eigSetting;
    tools::log->trace("Finding subspace -- partial");


    t_eig->tic();
    // You need to copy the data into StlMatrixProduct, because the PartialPivLU will overwrite the data in H_local otherwise.
    StlMatrixProduct<Scalar> hamiltonian(H_local.data(),H_local.rows(),Form::SYMMETRIC,Side::R, true);
    hamiltonian.set_shift(energy_target);
    hamiltonian.FactorOP();
    double t_lu = hamiltonian.t_factorOp.get_last_time_interval();
    t_eig->toc();

    double max_overlap_threshold      = optMode == OptMode::OVERLAP ? 0.9 : 1; //1.0/std::sqrt(2); //Slightly less than 1/sqrt(2), in case that the choice is between cat states.
    double subspace_quality_threshold = settings::precision::MinSubspaceQuality;

    class_eigsolver solver;
    std::string reason = "exhausted";
    Eigen::VectorXd  eigvals;
    Eigen::MatrixXcd eigvecs;
    Eigen::Map<const Eigen::VectorXcd> theta_vec   (theta.data(),theta.size());
    for (auto nev : generate_size_list(theta.size())){
        t_eig->tic();
        solver.eigs_stl(hamiltonian,nev,-1, energy_target,Form::SYMMETRIC,Ritz::LM,Side::R, true,false);
        t_eig->toc();

        eigvals = Eigen::Map<const Eigen::VectorXd > (solver.solution.get_eigvals<Form::SYMMETRIC>().data()      ,solver.solution.meta.cols);
        if constexpr (std::is_same<std::complex<double>, Scalar >::value){
            eigvecs = Eigen::Map<const Eigen::MatrixXcd> (solver.solution.get_eigvecs<Type::CPLX, Form::SYMMETRIC>().data(),solver.solution.meta.rows,solver.solution.meta.cols);
        }else{
            eigvecs = Eigen::Map<const Eigen::MatrixXd> (solver.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solver.solution.meta.rows,solver.solution.meta.cols);
        }

        Eigen::VectorXd overlaps = (theta_vec.adjoint() * eigvecs).cwiseAbs().real();
        double max_overlap       = overlaps.maxCoeff();
        double min_overlap       = overlaps.minCoeff();
        double sq_sum_overlap    = overlaps.cwiseAbs2().sum();
        double subspace_quality  = 1.0 - sq_sum_overlap;
        eig_log.emplace_back(nev, max_overlap, min_overlap, sq_sum_overlap, std::log10(subspace_quality), t_eig->get_last_time_interval(), t_lu);
        t_lu = 0;
        if(max_overlap            > 1.0 + 1e-6)                  throw std::runtime_error("max_overlap larger than one : "  + std::to_string(max_overlap));
        if(sq_sum_overlap         > 1.0 + 1e-6)                  throw std::runtime_error("eps larger than one : "          + std::to_string(sq_sum_overlap));
        if(min_overlap            < 0.0)                         throw std::runtime_error("min_overlap smaller than zero: " + std::to_string(min_overlap));
        if(max_overlap            >= max_overlap_threshold )    {reason = "overlap is good enough"; break;}
        if(subspace_quality       < subspace_quality_threshold) {reason = "subspace quality is good enough"; break;}
    }
    tools::log->debug("Finished partial eigensolver -- condition: {}",reason);
    return std::make_tuple(eigvecs,eigvals);
}





template<typename Scalar>
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
find_subspace(const class_finite_state & state,OptMode optMode){
    tools::log->trace("Finding subspace");

    using namespace eigutils::eigSetting;
    t_ham->tic();

    MatrixType<Scalar> H_local;
    if constexpr(std::is_same<Scalar,double>::value){
        H_local = state.get_multi_hamiltonian_matrix().real();
    }
    if constexpr(std::is_same<Scalar,std::complex<double>>::value){
        H_local = state.get_multi_hamiltonian_matrix();
    }

    if(not H_local.isApprox(H_local.adjoint(), 1e-14)){
        throw std::runtime_error(fmt::format("H_local is not hermitian: {:.16f}", (H_local - H_local.adjoint()).cwiseAbs().sum()));
    }
    double sparcity = (H_local.array().cwiseAbs2() != 0.0).count()/(double)H_local.size();
    tools::log->debug("H_local nonzeros: {:.8f} %", sparcity*100);

    t_ham->toc();
    auto theta = state.get_multitheta();


    Eigen::MatrixXcd eigvecs;
    Eigen::VectorXd  eigvals;
    std::vector<reports::eig_tuple> eig_log;
    double energy_target = tools::finite::measure::multisite::energy(state,theta);

    // If theta is small enough you can afford full diag.
    if   ((size_t)theta.size() <= settings::precision::MaxSizeFullDiag) {
        std::tie(eigvecs, eigvals) = find_subspace_full(H_local, theta, eig_log);
    }else{
        std::tie(eigvecs, eigvals) = find_subspace_part(H_local,theta,energy_target,eig_log,optMode);
    }
    reports::print_report(eig_log);

    if constexpr(std::is_same<Scalar,double>::value){
        Textra::subtract_phase(eigvecs);
        tools::log->trace("truncating imag of eigvecs, sum: {}", eigvecs.imag().cwiseAbs().sum() );
        eigvecs = eigvecs.real();
    }

    eigvecs.colwise().normalize();
    return std::make_tuple(eigvecs, eigvals);
}



Eigen::Tensor<class_finite_state::Scalar,3>
tools::finite::opt::internals::ceres_subspace_optimization(const class_finite_state &state,
                                                         const class_simulation_status &sim_status, OptType optType,
                                                         OptMode optMode){
    tools::log->trace("Optimizing in SUBSPACE mode");
    using Scalar = class_finite_state::Scalar;
    using namespace eigutils::eigSetting;

    auto theta             = state.get_multitheta();
    auto theta_old         = Eigen::Map<const Eigen::VectorXcd>  (theta.data(),theta.size());
    Eigen::MatrixXcd eigvecs;
    Eigen::VectorXd  eigvals;
    switch(optType){
        case OptType::CPLX:     std::tie (eigvecs,eigvals)  = find_subspace<Scalar>(state,optMode); break;
        case OptType::REAL:     std::tie (eigvecs,eigvals)  = find_subspace<double>(state,optMode); break;
    }
    tools::log->trace("Subspace found with {} eigenvectors", eigvecs.cols());
    Eigen::VectorXd overlaps = (theta_old.adjoint() * eigvecs).cwiseAbs().real();



    switch (optMode){
        case OptMode::OVERLAP:
        {
            int idx = idx_best_overlap_in_window(overlaps,eigvals,sim_status.energy_lbound,sim_status.energy_ubound);
            return Textra::Matrix_to_Tensor(eigvecs.col(idx), state.active_dimensions());
        }

        case OptMode::VARIANCE:
        {
            double sq_sum_overlap    = overlaps.cwiseAbs2().sum();
            double subspace_quality  = 1.0 - sq_sum_overlap;
            if(subspace_quality > settings::precision::MinSubspaceQuality) {
                tools::log->debug("Log subspace quality is poor: {} > {}. Switching to direct mode", std::log10(subspace_quality), std::log10(settings::precision::MinSubspaceQuality));
                return ceres_direct_optimization(state, sim_status, optType);
            }else{
                std::tie(eigvecs,eigvals,subspace_quality) = filter_states(eigvecs,eigvals,overlaps,settings::precision::MinSubspaceQuality, 128);
            }
        }

    }


    Eigen::VectorXcd theta_new;
    double overlap_new  = 0;
    double energy_new,variance_new,norm;
    //Should really use theta_start as the projection towards the previous theta, not best overlapping!
    // Note that alpha_i = <theta_old | theta_new_i> is not supposed to be squared! The overlap
    // Between theta_start and theta_old should be
    Eigen::VectorXcd theta_start      = (eigvecs.adjoint()  * theta_old).normalized()  ;

    std::vector<reports::subspc_opt_tuple> opt_log;


    t_opt->tic();
    state.unset_measurements();
    double energy_0   = tools::finite::measure::energy_per_site(state);
    double variance_0 = tools::finite::measure::energy_variance_per_site(state);
    t_opt->toc();
    Eigen::VectorXcd theta_0 = (eigvecs * theta_start.asDiagonal()).rowwise().sum().normalized();

    int iter_0 = 0;
    double overlap_0 = std::abs(theta_old.dot(theta_0));

    opt_log.emplace_back("Initial",theta.size(), energy_0, std::log10(variance_0), overlap_0, iter_0,0, t_opt->get_last_time_interval());

    ceres::GradientProblemSolver::Options options;
    options.line_search_type = ceres::LineSearchType::WOLFE;
    options.line_search_interpolation_type = ceres::LineSearchInterpolationType::CUBIC;
    options.line_search_direction_type = ceres::LineSearchDirectionType::LBFGS;
    options.nonlinear_conjugate_gradient_type = ceres::NonlinearConjugateGradientType::POLAK_RIBIERE;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = 150;
    options.max_lbfgs_rank     = 250;
    options.use_approximate_eigenvalue_bfgs_scaling = true;
    options.max_line_search_step_expansion = 100.0;
    options.min_line_search_step_size = 1e-8;
    options.max_line_search_step_contraction = 1e-3;
    options.min_line_search_step_contraction = 0.6;
    options.max_num_line_search_step_size_iterations  = 20;
    options.max_num_line_search_direction_restarts    = 2;
    options.line_search_sufficient_function_decrease  = 1e-2;
    options.line_search_sufficient_curvature_decrease = 0.5;
    options.max_solver_time_in_seconds = 60*2;
    options.function_tolerance = 1e-4;
    options.gradient_tolerance = 1e-8;
    options.parameter_tolerance = 1e-12;
    ceres::GradientProblemSolver::Summary summary;

    t_opt->tic();
    using namespace LBFGSpp;
    using namespace tools::finite::opt::internals;
    int counter,iter;
    LBFGSpp::LBFGSSolver<double> solver(params);
    switch (optType){
        case OptType::CPLX:{
            Eigen::VectorXd  theta_start_cast = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double*> (theta_start.data()), 2*theta_start.size());
            auto * functor = new ceres_subspace_functor<std::complex<double>>(state, sim_status,eigvecs,eigvals);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running L-BFGS");
            ceres::Solve(options, problem, theta_start_cast.data(), &summary);
            iter         = (int)summary.iterations.size();
            counter      = functor->get_count();
            norm         = functor->get_norm();
            energy_new   = functor->get_energy();
            variance_new = functor->get_variance();
            theta_start  = Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<Scalar*> (theta_start_cast.data()), theta_start_cast.size()/2).normalized();
            theta_new    = (eigvecs * theta_start.asDiagonal()).rowwise().sum().normalized();
            break;
        }
        case OptType::REAL:{
            Eigen::VectorXd  theta_start_cast = theta_start.real();
            auto * functor = new ceres_subspace_functor<double>(state, sim_status,eigvecs.real(),eigvals);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running LBFGS");
            ceres::Solve(options, problem, theta_start_cast.data(), &summary);
            iter         = (int)summary.iterations.size();
            counter      = functor->get_count();
            norm         = functor->get_norm();
            energy_new   = functor->get_energy();
            variance_new = functor->get_variance();
            theta_start  = theta_start_cast.normalized().cast<Scalar>();
            theta_new    = (eigvecs.real() * theta_start.real().asDiagonal()).rowwise().sum().normalized();
            break;
        }
    }
    t_opt->toc();

    overlap_new = (theta_old.adjoint() * theta_new).cwiseAbs().sum();
    opt_log.emplace_back("Ceres L-BFGS",theta.size(), energy_new, std::log10(variance_new), overlap_new, iter,counter, t_opt->get_last_time_interval());
    tools::log->trace("Finished Ceres. Exit status: {}. Message: {}", ceres::TerminationTypeToString(summary.termination_type) , summary.message.c_str());
    //    std::cout << summary.FullReport() << "\n";
    reports::print_report(opt_log);


    if (variance_new < variance_0 * 2.0){
        state.unset_measurements();
        tools::log->debug("Returning new theta");
        return  Textra::Matrix_to_Tensor(theta_new, state.active_dimensions());

    }else{
        tools::log->debug("Returning old theta");
        return  theta;
    }


}






