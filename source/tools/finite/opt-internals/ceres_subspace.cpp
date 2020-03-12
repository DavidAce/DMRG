//
// Created by david on 2019-07-15.
//

#include "ceres_subspace_functor.h"
#include <iostream>
#include <simulation/class_simulation_status.h>
#include <math/class_eigsolver.h>
#include <math/arpack_extra/matrix_product_stl.h>
#include <state/class_state_finite.h>
#include <simulation/nmspc_settings.h>
#include <math/nmspc_random.h>
#include <tools/finite/measure.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
using namespace tools::finite::opt;
using namespace tools::finite::opt::internal;

template<typename T> using MatrixType = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;


std::vector<int> tools::finite::opt::internal::generate_size_list(size_t shape){
    int max_nev ;
    if      (shape <= 512)  {max_nev = shape/2;}
    else if (shape > 512  and shape <= 1024) {max_nev = shape/4;} // should do full diag
    else if (shape > 1024 and shape <= 2048) {max_nev = shape/4;} // should do full diag
    else if (shape > 2048 and shape <= 4096) {max_nev = 64;}
    else if (shape > 4096 and shape <= 8192) {max_nev = 32;}
    else                                     {max_nev = 16;}

    int min_nev = std::min(std::min(8,(int)shape),max_nev);

    std::vector<int> nev_list = {min_nev};
    int tmp_nev = min_nev;
    while (tmp_nev < max_nev){
        tmp_nev = std::min(4*tmp_nev, max_nev);
        nev_list.push_back(tmp_nev);
    }
    return nev_list;
}


std::tuple<Eigen::MatrixXcd,Eigen::VectorXd,Eigen::VectorXd,double> filter_states(const Eigen::MatrixXcd &eigvecs, const Eigen::VectorXd& eigvals, Eigen::VectorXd &overlaps, double maximum_subspace_error, size_t max_accept){

    size_t min_accept = std::min(8ul,(size_t)eigvals.size());
    max_accept        = std::min(max_accept,(size_t)eigvals.size());
    if(min_accept == max_accept) return std::make_tuple(eigvecs, eigvals,overlaps, 1.0 - overlaps.cwiseAbs2().sum());
    Eigen::VectorXd     overlaps_filtered = overlaps;
    std::vector<int>    overlaps_accepted_idx;
    std::vector<double> overlaps_accepted;
    double epsilon           = std::numeric_limits<double>::epsilon();
    double subspace_error    = 1.0 - overlaps.cwiseAbs2().sum();
    maximum_subspace_error   = epsilon + std::min(subspace_error, maximum_subspace_error); //Make sure you don't actually increase the allowed subspace error


    while(true){
        int idx;
        double overlap = overlaps_filtered.maxCoeff(&idx);
        overlaps_accepted_idx.push_back(idx);
        overlaps_accepted    .push_back(overlap);
        Eigen::Map<Eigen::VectorXd> overlaps_map(overlaps_accepted.data(),overlaps_accepted.size());
        subspace_error  = 1.0 - overlaps_map.cwiseAbs2().sum();
        if(overlaps_accepted.size() >= min_accept){
            if(subspace_error < maximum_subspace_error) break;
            if(overlaps_accepted.size() >= max_accept) break;
        }
        overlaps_filtered(idx) = 0;
        if(overlaps_filtered.sum() == 0) break;
    }

    Eigen::MatrixXcd eigvecs_filtered(eigvecs.rows(),overlaps_accepted.size());
    Eigen::VectorXd  eigvals_filtered(overlaps_accepted.size());
    overlaps_filtered.resize(overlaps_accepted.size());
    int col_num = 0;
    for (auto &idx : overlaps_accepted_idx){
        eigvecs_filtered.col(col_num) = eigvecs.col(idx);
        eigvals_filtered    (col_num) = eigvals(idx);
        overlaps_filtered   (col_num) = overlaps(idx);
        col_num++;
    }
    tools::log->trace("Filtered from {} down to {} states", eigvals.size(), eigvals_filtered.size());
    tools::log->trace("Subspace error after filter log10(1-eps) = {}", std::log10(epsilon + subspace_error));
    return std::make_tuple(eigvecs_filtered, eigvals_filtered,overlaps_filtered, subspace_error);
}


std::pair<double,int>
get_best_variance_in_window(const class_state_finite &state, const Eigen::MatrixXcd &eigvecs, const Eigen::VectorXd & energies_per_site, double lbound, double ubound){
    Eigen::VectorXd variances(eigvecs.cols());
    for(long idx = 0; idx < eigvecs.cols(); idx++){
        if (energies_per_site(idx) <=  ubound and energies_per_site(idx) >= lbound ) {
            auto multitheta = Textra::MatrixTensorMap(eigvecs.col(idx), state.active_dimensions());
            variances(idx)  = tools::finite::measure::energy_variance_per_site(state, multitheta);
        }else{
            variances(idx) = std::numeric_limits<double>::infinity();
        }
    }

    if (variances.minCoeff() == std::numeric_limits<double>::infinity()) {
        tools::log->debug("No eigenstates with good variance in given energy window {} to {}.", lbound,ubound);
        tools::log->debug("Subspace energy range is {} to {}.", energies_per_site.minCoeff(), energies_per_site.maxCoeff());
        return std::make_pair(std::numeric_limits<double>::quiet_NaN(), -1);
    }
    int    min_variance_idx;
    double min_variance_val = variances.minCoeff(&min_variance_idx);
    return std::make_pair(min_variance_val, min_variance_idx);
}

std::pair<double,int> get_best_overlap_in_window(const Eigen::VectorXd &overlaps, const Eigen::VectorXd & energies_per_site, double lbound, double ubound){
    assert(overlaps.size() == energies_per_site.size() and "get_best_overlap_in_window: Mismatch in overlaps and energies_per_site sizes");
    Eigen::VectorXd overlaps_in_window = overlaps;
    for (long i = 0; i < overlaps.size(); i++){
        if (energies_per_site(i) > ubound) overlaps_in_window(i) = 0.0;
        if (energies_per_site(i) < lbound) overlaps_in_window(i) = 0.0;
    }
    if (overlaps_in_window.maxCoeff() == 0.0){
        tools::log->debug("No overlapping eigenstates in given energy window {} to {}.", lbound,ubound);
        tools::log->debug("Subspace energy range is {} to {}.", energies_per_site.minCoeff(), energies_per_site.maxCoeff());
        return std::make_pair(std::numeric_limits<double>::quiet_NaN() , -1);
    }

    int    max_overlap_idx;
    double max_overlap          = overlaps_in_window.maxCoeff(&max_overlap_idx);
    return std::make_pair(max_overlap,max_overlap_idx);

}



std::vector<std::pair<double,int>> get_best_candidates_in_window(const Eigen::VectorXd &overlaps, const Eigen::VectorXd & energies_per_site, double lbound, double ubound){
    assert(overlaps.size() == energies_per_site.size() and "get_best_overlap_in_window: Mismatch in overlaps and energies_per_site sizes");
    std::vector<std::pair<double,int>> overlaps_in_window;
    for (long i = 0; i < overlaps.size(); i++){
        if (energies_per_site(i) < ubound and energies_per_site(i) > lbound)
            overlaps_in_window.emplace_back(std::make_pair(overlaps(i),i));
    }
    if (overlaps_in_window.empty()){
        tools::log->debug("No candidate eigenstates in given energy window {} to {}.", lbound,ubound);
        tools::log->debug("Subspace energy range is {} to {}.", energies_per_site.minCoeff(), energies_per_site.maxCoeff());
        return overlaps_in_window;
    }

    std::sort(overlaps_in_window.begin(),overlaps_in_window.end()); // Sort in ascending order of overlap
    std::vector<std::pair<double,int>> candidates;

    // We have a list of states in the energy window. We collect "candidates"
    // from this list until the squared sum of their overlaps goes above a certain threshold.
    // Note that the squared sum overlap = 1 if the states in the list form a complete basis for the current state.
    tools::log->debug("Collecting initial guesses");
    auto lambda_sq_sum = [&](double acc, std::pair<double,int> & p){return acc + p.first * p.first; };
    while(true){
        if(overlaps_in_window.empty()) break;
        double sq_sum_overlap = std::accumulate(candidates.begin(),candidates.end(), 0.0, lambda_sq_sum);
        tools::log->debug("Sq_sum_overlap:  {:.16f}",sq_sum_overlap);
        if(sq_sum_overlap  > 0.6) break; // Half means cat state.
        else {
            candidates.emplace_back(overlaps_in_window.back());
            overlaps_in_window.pop_back();
        }
    }
    tools::log->debug("Found {} candidates with good overlap in energy window",candidates.size());
    return candidates;

}


template<typename Scalar>
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
find_subspace_full(const MatrixType<Scalar> & H_local, Eigen::Tensor<std::complex<double>,3> &theta){
    using namespace eigutils::eigSetting;
    using namespace tools::common::profile;
    tools::log->trace("Finding subspace -- full");
    Eigen::VectorXd   eigvals;
    Eigen::MatrixXcd  eigvecs;
    class_eigsolver   solver;

    t_eig->tic();
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
    tools::log->debug("Finished eigensolver -- reason: Full diagonalization");

    Eigen::Map<const Eigen::VectorXcd> theta_vec   (theta.data(),theta.size());
    Eigen::VectorXd overlaps = (theta_vec.adjoint() * eigvecs).cwiseAbs().real();
    int idx;
    double max_overlap       = overlaps.maxCoeff(&idx);
    double min_overlap       = overlaps.minCoeff();
    double sq_sum_overlap    = overlaps.cwiseAbs2().sum();
    double subspace_error    = 1.0 - sq_sum_overlap;
    int nev = eigvecs.cols();
    reports::eigs_log.emplace_back(nev, max_overlap, min_overlap, std::log10(subspace_error),t_ham->get_last_time_interval(), t_eig->get_last_time_interval(), 0);
    return std::make_tuple(eigvecs,eigvals);
}



template<typename Scalar>
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
find_subspace_part(const MatrixType<Scalar> & H_local, Eigen::Tensor<std::complex<double>,3> &theta, double energy_target, double subspace_error_threshold, OptMode optMode, OptSpace optSpace){
    using namespace eigutils::eigSetting;
    using namespace tools::common::profile;
    tools::log->trace("Finding subspace -- partial");
    // You need to copy the data into StlMatrixProduct, because the PartialPivLU will overwrite the data in H_local otherwise.
    StlMatrixProduct<Scalar> hamiltonian(H_local.data(),H_local.rows(),true,Form::SYMMETRIC,Side::R);
    hamiltonian.set_shift(energy_target);
    hamiltonian.FactorOP();
    double time_lu  = hamiltonian.t_factorOp.get_last_time_interval();
    double time_ham = tools::common::profile::t_ham->get_last_time_interval();
    tools::common::profile::t_eig->toc();

    class_eigsolver solver;
    solver.solverConf.eigThreshold = settings::precision::eig_threshold;
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
        double subspace_error    = 1.0 - sq_sum_overlap;
        reports::eigs_log.emplace_back(nev, max_overlap, min_overlap, std::log10(subspace_error), time_ham, t_eig->get_last_time_interval(), time_lu);
        time_lu  = 0;
        time_ham = 0;
        if(max_overlap            > 1.0 + 1e-6)                                  throw std::runtime_error("max_overlap larger than one : "  + std::to_string(max_overlap));
        if(sq_sum_overlap         > 1.0 + 1e-6)                                  throw std::runtime_error("eps larger than one : "          + std::to_string(sq_sum_overlap));
        if(min_overlap            < 0.0)                                         throw std::runtime_error("min_overlap smaller than zero: " + std::to_string(min_overlap));
        if(subspace_error < subspace_error_threshold)                           {reason = "subspace error is low enough"; break;}
        if(optSpace == OptSpace::SUBSPACE_AND_DIRECT and subspace_error < 1e-3) {reason = "subspace error sufficient for SUBSPACE_AND_DIRECT"; break;}
        if(optSpace == OptSpace::SUBSPACE_ONLY and optMode == OptMode::OVERLAP  and max_overlap >= 1.0/std::sqrt(2.0))
            {reason = "Overlap sufficient for OVERLAP and SUBSPACE_ONLY"; break;}
    }
    tools::log->debug("Finished partial eigensolver -- reason: {}",reason);
    tools::common::profile::t_eig->toc();
    return std::make_tuple(eigvecs,eigvals);
}





template<typename Scalar>
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
find_subspace(const class_state_finite & state, double subspace_error_threshold, OptMode optMode,OptSpace optSpace){
    tools::log->trace("Finding subspace");

    using namespace eigutils::eigSetting;
    tools::common::profile::t_ham->tic();
    MatrixType<Scalar> H_local = tools::finite::opt::internal::get_multi_hamiltonian_matrix<Scalar>(state);
    tools::common::profile::t_ham->toc();
    auto multitheta = state.get_multitheta();


    Eigen::MatrixXcd eigvecs;
    Eigen::VectorXd  eigvals;


    // If multitheta is small enough you can afford full diag.
    if   ((size_t)multitheta.size() <= settings::precision::max_size_full_diag) {
        std::tie(eigvecs, eigvals) = find_subspace_full(H_local, multitheta);
    }else{
        double energy_target;
        if (state.isReduced()) energy_target = tools::finite::measure::energy_minus_energy_reduced(state, multitheta);
        else                   energy_target = tools::finite::measure::energy(state, multitheta);
//        tools::log->debug("Energy target, per site: {}",energy_target/state.get_length());
        tools::log->trace("Energy target + energy reduced = energy per site: {} + {} = {}",
                energy_target/state.get_length(),
                state.get_energy_reduced()/state.get_length(),
                (energy_target + state.get_energy_reduced())/state.get_length());

        std::tie(eigvecs, eigvals) = find_subspace_part(H_local, multitheta, energy_target, subspace_error_threshold, optMode,optSpace);
    }
    tools::log->trace("Eigenvalue range: {} --> {}",
            (eigvals.minCoeff() + state.get_energy_reduced())/state.get_length(),
            (eigvals.maxCoeff() + state.get_energy_reduced())/state.get_length());
//    eigvals = eigvals.array() + state.get_energy_reduced();
//    tools::log->debug("Eigenvalue range: {} --> {}", eigvals.minCoeff()/state.get_length(),eigvals.maxCoeff()/state.get_length());
    reports::print_eigs_report();

    if constexpr(std::is_same<Scalar,double>::value){
        Textra::subtract_phase(eigvecs);
        tools::log->trace("truncating imag of eigvecs, sum: {}", eigvecs.imag().cwiseAbs().sum() );
        eigvecs = eigvecs.real();
    }

    return std::make_tuple(eigvecs, eigvals);
}




Eigen::Tensor<class_state_finite::Scalar,3>
tools::finite::opt::internal::ceres_subspace_optimization(const class_state_finite &state, const class_simulation_status &sim_status, OptType optType, OptMode optMode, OptSpace optSpace){
    tools::log->trace("Optimizing in SUBSPACE mode");
    tools::common::profile::t_opt->tic();
    using Scalar = class_state_finite::Scalar;
    using namespace eigutils::eigSetting;
    auto options = ceres_default_options;
    double theta_old_variance    = tools::finite::measure::energy_variance_per_site(state);
    double subspace_error_threshold     = settings::precision::min_subspace_error;


    auto & theta_old               = state.get_multitheta();
    auto theta_old_vec             = Eigen::Map<const Eigen::VectorXcd>  (theta_old.data(), theta_old.size());

    Eigen::MatrixXcd eigvecs;
    Eigen::VectorXd  eigvals;
    switch(optType){
        case OptType::CPLX:     std::tie (eigvecs,eigvals)  = find_subspace<Scalar>(state,subspace_error_threshold,optMode,optSpace); break;
        case OptType::REAL:     std::tie (eigvecs,eigvals)  = find_subspace<double>(state,subspace_error_threshold,optMode,optSpace); break;
    }
    Eigen::VectorXd eigvals_per_site_unreduced = (eigvals.array() + state.get_energy_reduced())/state.get_length(); // Remove energy reduction for energy window comparisons
    tools::log->trace("Subspace found with {} eigenvectors", eigvecs.cols());
    Eigen::VectorXd overlaps = (theta_old_vec.adjoint() * eigvecs).cwiseAbs().real();


    int     idx;
    double max_overlap          = overlaps.maxCoeff(&idx);
    double max_overlap_energy   = eigvals_per_site_unreduced(idx);
    bool   max_overlap_inwindow = sim_status.energy_lbound < max_overlap_energy and max_overlap_energy < sim_status.energy_ubound;
    tools::log->trace("Max overlap: {} -- Energy per site: {} -- Idx: {} -- inside of window: {}", max_overlap, max_overlap_energy, idx,max_overlap_inwindow );






    // For options LC - E we need to filter down the set of states in case we do subspace optimization, otherwise we can easily run out of memory. 64 candidates should do it.
//    double subspace_error_unfiltered = 1.0 - overlaps.cwiseAbs2().sum();
    double subspace_error_filtered;

    std::tie(eigvecs,eigvals,overlaps,subspace_error_filtered) = filter_states(eigvecs, eigvals, overlaps, subspace_error_threshold, 128);
    eigvals_per_site_unreduced = (eigvals.array() + state.get_energy_reduced())/state.get_length(); // Remove energy reduction for energy window comparisons
//    bool force_accept = false;


    tools::log->trace("Current energy          : {:.16f}", tools::finite::measure::energy_per_site(state));
    tools::log->trace("Current energy (2site)  : {:.16f}", tools::finite::measure::twosite::energy_per_site(state,state.get_theta()));
    tools::log->trace("Current energy (multi)  : {:.16f}", tools::finite::measure::multisite::energy_per_site(state,state.get_multitheta()));
    tools::log->trace("Current variance: {:.16f}", std::log10(theta_old_variance) );
//    auto [best_overlap,best_overlap_idx] = get_best_overlap_in_window(overlaps, eigvals_per_site_unreduced, sim_status.energy_lbound, sim_status.energy_ubound);
    if (optMode == OptMode::OVERLAP){
        auto [best_overlap,best_overlap_idx]   = get_best_overlap_in_window(overlaps, eigvals_per_site_unreduced, sim_status.energy_lbound, sim_status.energy_ubound);
        if (best_overlap_idx  < 0 ){
            //Option A
//            tools::log->trace("No overlapping states in energy range. Returning old theta");
            tools::log->info("No overlapping states in energy range. Returning best overlap out of window");
//            auto   best_overlap_theta              = Textra::MatrixTensorMap(eigvecs.col(0), state.active_dimensions());
//            state.tag_active_sites_have_been_updated(false);
//            return best_overlap_theta;
            return theta_old;
        }else if (best_overlap < 0.0){
            //Overlap is too, bad, just go to the next site and hope for something better to come along
            //Turn this option off with best_overlap < 0.0
            tools::log->info("Overlap too low, returning old theta");
            return theta_old;
        }else{
            auto   best_overlap_theta              = Textra::MatrixTensorMap(eigvecs.col(best_overlap_idx), state.active_dimensions());
            double best_overlap_energy             = eigvals_per_site_unreduced(best_overlap_idx);
            double best_overlap_variance           = tools::finite::measure::energy_variance_per_site(state, best_overlap_theta);
            tools::log->trace("Candidate {:<2} has highest overlap: Overlap: {:.16f} Energy: {:>20.16f} Variance: {:>20.16f}",
                    best_overlap_idx ,overlaps(best_overlap_idx) ,best_overlap_energy  ,std::log10(best_overlap_variance) );
            return best_overlap_theta;
        }

    }

    auto list_of_candidates = get_best_candidates_in_window(overlaps, eigvals_per_site_unreduced, sim_status.energy_lbound, sim_status.energy_ubound);
    if (list_of_candidates.empty()){
        //Option A
        tools::log->warn("Went for option A -- No overlapping states in energy range. Returning old theta");
        state.tag_active_sites_have_been_updated(false);
        return theta_old;
    }



    // We should try different initial guesses when doing subspace optimization
    // All the candidates with significant overlap should be considered, including the current theta
    std::vector<Eigen::Tensor<Scalar,3>> theta_candidates = {theta_old};
    for(auto &candidate : list_of_candidates){
        auto   candidate_theta              = Textra::MatrixTensorMap(eigvecs.col(candidate.second), state.active_dimensions());
        double candidate_energy             = eigvals_per_site_unreduced(candidate.second);
        double candidate_variance           = tools::finite::measure::energy_variance_per_site(state, candidate_theta);
        tools::log->trace("Candidate {:<2} has good overlap: Overlap: {:.16f} Energy: {:>20.16f} Variance: {:>20.16f}",candidate.second ,candidate.first ,candidate_energy  ,std::log10(candidate_variance) );
        theta_candidates.emplace_back(candidate_theta);
        if(theta_candidates.size() > 5) break;
    }

    if (tools::log->level() <= spdlog::level::debug ){
        // Initial sanity check
        tools::common::profile::t_opt->tic();
        double energy_old        = tools::finite::measure::energy_per_site(state,theta_old);
        double variance_old      = tools::finite::measure::energy_variance_per_site(state,theta_old);
        tools::common::profile::t_opt->toc();
        reports::bfgs_log.emplace_back("Current state", theta_old.size(), 0, energy_old, std::log10(variance_old), 1.0 , theta_old_vec.norm(), 1,1, tools::common::profile::t_opt->get_last_time_interval());
    }

    tools::log->debug("Optimizing mode {} space {} using subspace size {} with {} initial guesses", optMode, optSpace, eigvecs.cols(), theta_candidates.size() );
    tools::common::profile::t_opt->tic();
    Eigen::MatrixXcd H2_subspace = tools::finite::opt::internal::get_multi_hamiltonian_squared_subspace_matrix_new<Scalar>(state, eigvecs);
    if(optType == OptType::REAL) H2_subspace = H2_subspace.real();


    tools::common::profile::t_opt->toc();
    double t_H2_subspace = tools::common::profile::t_opt->get_last_time_interval();
    std::vector<std::pair<double,Eigen::Tensor<Scalar,3>>> optimized_results;
    size_t candidate_count = 0;
    for(auto &theta_candidate : theta_candidates){
        auto theta_candidate_map = Eigen::Map<const Eigen::VectorXcd>  (theta_candidate.data(), theta_candidate.size());
        Eigen::VectorXcd theta_new;
        double overlap_new  = 0;
        double energy_new,variance_new;
        [[maybe_unused]] double norm;
        // Note that alpha_i = <theta_initial | theta_new_i> is not supposed to be squared!
        Eigen::VectorXcd theta_start      = (eigvecs.adjoint() * theta_candidate_map).normalized()  ;

        if (tools::log->level() <= spdlog::level::trace){
            // Check current candidate
            tools::common::profile::t_opt->tic();
            Eigen::VectorXcd theta_0 = (eigvecs * theta_start.asDiagonal() ).rowwise().sum().normalized();
            auto theta_0_tensor      = Textra::MatrixTensorMap(theta_0,state.active_dimensions());
            double energy_0          = tools::finite::measure::energy_per_site(state,theta_0_tensor);
            double variance_0        = tools::finite::measure::energy_variance_per_site(state,theta_0_tensor);
            double overlap_0         = std::abs(theta_old_vec.dot(theta_0));
            tools::common::profile::t_opt->toc();
            reports::bfgs_log.emplace_back("Candidate " + std::to_string(candidate_count) , theta_old.size(),0, energy_0, std::log10(variance_0), overlap_0, theta_0.norm(), 1,1, tools::common::profile::t_opt->get_last_time_interval());
        }
        if (tools::log->level() <= spdlog::level::trace and settings::debug){
            // Check using explicit matrix
            tools::common::profile::t_opt->tic();
            Eigen::VectorXcd theta_0 = (eigvecs * theta_start.asDiagonal() ).rowwise().sum().normalized();
            double overlap_0         = std::abs(theta_old_vec.dot(theta_0));
            Eigen::VectorXcd Hv  = eigvals.asDiagonal() * theta_start;
            Eigen::VectorXcd H2v = H2_subspace.template selfadjointView<Eigen::Upper>()*theta_start;
            Scalar vHv  = theta_start.dot(Hv);
            Scalar vH2v = theta_start.dot(H2v);
            double vv   = theta_start.squaredNorm();
            Scalar ene  = vHv/vv;
            Scalar var  = vH2v/vv - ene*ene;
            double ene_init_san = std::real(ene+state.get_energy_reduced())/state.get_length();
            double var_init_san = std::real(var)/state.get_length();
            tools::common::profile::t_opt->toc();
            reports::bfgs_log.emplace_back("Candidate " + std::to_string(candidate_count) + " (matrix)",theta_start.size(),0, ene_init_san, std::log10(var_init_san), overlap_0,theta_start.norm(), 1,1, t_H2_subspace+tools::common::profile::t_opt->get_last_time_interval());
//            if(not H2_subspace.isApprox(H2_subspace_old,1e-4)){
//                std::cout << "H2 new = \n" << H2_subspace.topLeftCorner(6,6) << std::endl;
//                std::cout << "H2 old = \n" << H2_subspace_old.topLeftCorner(6,6) << std::endl;
//                tools::log->warn("H2 subspace mismatch: {:.16f}", (H2_subspace - H2_subspace_old).cwiseAbs().sum());
//            }
        }




        ceres::GradientProblemSolver::Summary summary;
        tools::common::profile::t_opt->tic();
        using namespace tools::finite::opt::internal;
        int counter,iter;
        switch (optType){
            case OptType::CPLX:{
                Eigen::VectorXd  theta_start_cast = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double*> (theta_start.data()), 2*theta_start.size());
                auto * functor = new ceres_subspace_functor<std::complex<double>>(state, sim_status,H2_subspace,eigvals);
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
                Eigen::MatrixXd H2_subspace_real = H2_subspace.real();
                auto * functor = new ceres_subspace_functor<double>(state, sim_status,H2_subspace_real,eigvals);
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
        tools::common::profile::t_opt->toc();
        reports::time_log.emplace_back(tools::common::profile::t_vH2v->get_measured_time_and_reset(),
                                       tools::common::profile::t_vHv->get_measured_time_and_reset(),
                                       tools::common::profile::t_vH2->get_measured_time_and_reset(),
                                       tools::common::profile::t_vH->get_measured_time_and_reset(),
                                       tools::common::profile::t_op->get_measured_time_and_reset());


        if (tools::log->level() <= spdlog::level::debug){
            // Print results of Ceres LBFGS
            overlap_new = (theta_old_vec.adjoint() * theta_new).cwiseAbs().sum();
            reports::bfgs_log.emplace_back("LBFGS subspace", theta_old.size(),options.max_lbfgs_rank, energy_new, std::log10(variance_new), overlap_new, theta_new.norm(), iter, counter, tools::common::profile::t_opt->get_last_time_interval());
        }
        if (tools::log->level() <= spdlog::level::trace and settings::debug){
            // Sanity check
            tools::common::profile::t_opt->tic();
            auto theta_san      = Textra::MatrixTensorMap(theta_new, state.active_dimensions());
            double energy_san   = tools::finite::measure::energy_per_site(state,theta_san);
            double variance_san = tools::finite::measure::energy_variance_per_site(state,theta_san);
            tools::common::profile::t_opt->toc();
            if(std::abs((variance_san - variance_new) / variance_san ) > 0.01 ) tools::log->warn("Variance mismatch in sanity check: {:.16f} != {:.16f}", variance_san, variance_new);
            reports::bfgs_log.emplace_back("Sanity check",theta_san.size(),0, energy_san, std::log10(variance_san), overlap_new,theta_new.norm(), 1,1, tools::common::profile::t_opt->get_last_time_interval());
        }

        tools::log->trace("Finished LBFGS after {} seconds ({} iters). Exit status: {}. Message: {}",summary.total_time_in_seconds, summary.iterations.size(), ceres::TerminationTypeToString(summary.termination_type) , summary.message.c_str());
        //    std::cout << summary.FullReport() << "\n";

        if(optSpace == OptSpace::SUBSPACE_ONLY)
        {
            auto optimized_theta    = Textra::MatrixTensorMap(theta_new, state.active_dimensions());
            auto optimized_energy   = tools::finite::measure::energy_per_site(state,optimized_theta);
            auto optimized_variance = tools::finite::measure::energy_variance_per_site(state,optimized_theta);
            auto optimized_vec      = Eigen::Map<const Eigen::VectorXcd>  (optimized_theta.data(),optimized_theta.size());
            auto optimized_overlap  = std::abs(theta_old_vec.dot(optimized_vec));
            optimized_results.emplace_back(std::make_pair(optimized_variance,optimized_theta));
            reports::bfgs_log.emplace_back("LBFGS subspace",optimized_theta.size(),0, optimized_energy, std::log10(optimized_variance), optimized_overlap,optimized_vec.norm(), 1,1, tools::common::profile::t_opt->get_last_time_interval());

        }
        if(optSpace == OptSpace::SUBSPACE_AND_DIRECT)
        {
            tools::log->trace("Fine tuning new theta after SUBSPACE optimization");
            auto optimized_theta    = ceres_direct_optimization(state, Textra::MatrixTensorMap(theta_new, state.active_dimensions()) ,sim_status, optType,optMode,optSpace);
            auto optimized_variance = tools::finite::measure::energy_variance_per_site(state,optimized_theta);
            optimized_results.emplace_back(std::make_pair(optimized_variance,optimized_theta));
        }
        candidate_count++;
    }

    // Finish up and print reports
    reports::print_bfgs_report();
    reports::print_time_report();


    //Sort thetas in ascending order in variance
    std::sort(optimized_results.begin(), optimized_results.end(), [](auto &left, auto &right) {return left.first < right.first;});
    //Return the best theta
    return optimized_results.front().second;




//    auto theta_direct_fine_tuned    =  ceres_direct_optimization(state, Textra::Matrix_to_Tensor(theta_new, state.active_dimensions()) ,sim_status, optType);
//    auto theta_direct_best_overlap  =  ceres_direct_optimization(state, best_overlap_theta ,sim_status, optType);
//    auto theta_direct_current_theta =  ceres_direct_optimization(state, theta_old ,sim_status, optType);
//
//    double variance_direct_fine_tuned       =  tools::finite::measure::energy_variance_per_site(state,theta_direct_fine_tuned   );
//    double variance_direct_best_overlap     =  tools::finite::measure::energy_variance_per_site(state,theta_direct_best_overlap );
//    double variance_direct_current_theta    =  tools::finite::measure::energy_variance_per_site(state,theta_direct_current_theta);
//    auto spin_components = tools::finite::measure::spin_components(state);
//    tools::log->debug("spin component x              = {:.16f}", spin_components[0] );
//    tools::log->debug("spin component y              = {:.16f}", spin_components[1] );
//    tools::log->debug("spin component z              = {:.16f}", spin_components[2] );
//    tools::log->debug("suspace_error                 = {:.16f}", std::log10(std::numeric_limits<double>::epsilon() + subspace_error_filtered) );
//    tools::log->debug("best_overlap                  = {:.16f}", best_overlap);
//    tools::log->debug("variance_original             = {:.16f}", std::log10(theta_old_variance           ));
//    tools::log->debug("variance_subspace_optimized   = {:.16f}", std::log10(variance_new                 ));
//    tools::log->debug("variance_direct_fine_tuned    = {:.16f}", std::log10(variance_direct_fine_tuned   ));
//    tools::log->debug("variance_direct_best_overlap  = {:.16f}", std::log10(variance_direct_best_overlap ));
//    tools::log->debug("variance_direct_current_theta = {:.16f}", std::log10(variance_direct_current_theta));




}







//
//
//    if (variance_new < 1.0 * tools::finite::measure::energy_variance_per_site(state)){
//        tools::log->debug("Returning new (better) theta");
//        state.tag_active_sites_have_been_updated(true);
//    }else{
//        tools::log->debug("Made it worse. Sending initial guess to DIRECT");
//        state.tag_active_sites_have_been_updated(false);
//
//    }
// Perhaps send theta initial to direct if worse?

//    return  Textra::Matrix_to_Tensor(theta_new, state.active_dimensions());
//
//    if (variance_new < 1.0 * tools::finite::measure::energy_variance_per_site(state)){
//        // Only an improvement of 1% is considered to be an actual improvement
//        tools::log->debug("Returning new (better) theta");
//        state.tag_active_sites_have_been_updated(true);
//        return  Textra::Matrix_to_Tensor(theta_new, state.active_dimensions());
//
//    }
//    else if (variance_new < 10.0 * tools::finite::measure::energy_variance_per_site(state)) {
//        // Allow for variance to increase a bit to come out of local minima
//        tools::log->debug("Returning new (but not good enough) theta");
//        state.tag_active_sites_have_been_updated(false);
//        return  Textra::Matrix_to_Tensor(theta_new, state.active_dimensions());
//    }
//    else{
//        tools::log->debug("Subspace optimization didn't improve variance.");
//        tools::log->debug("Returning old theta");
//        if (variance_new <= settings::precision::variance_convergence_threshold)
//              state.tag_active_sites_have_been_updated(true);
//        else  state.tag_active_sites_have_been_updated(false);
//        return  theta_old;
//
//    }











// Explanation:
// theta_initial: The starting point , or initial guess, for the gradient descent (L-BFGS) optimization routine.
//                By default theta_initial = theta_old, i.e. the current state.
// candidate    : One of the eigenvectors obtained from either full or partial diagonalization, i.e. lapack or arpack.
// relevant candidate : Eigenvectors inside of the energy window with high enough overlap with the old theta_old.
// subspace_error = 1 - Σ_i |<theta_new_i|theta_old>|^2
//      If == 0, it means that the set of candidate theta_old's span the old theta_old, i.e. the set can describe the current state.
//      The subspace error is "low enough" when subspace_error < subspace_error_threshold
//
// best_overlap : The highest overlap to the old theta_old achieved by any candidate inside of the energy window, i.e. |<candidate_i | theta_old>|_max
// best_overlap_idx: The index of the best overlapping candidate in the energy window. If -1, it means that no state is in the window.
// best_variance:

// overlap_high = 0.9
// overlap_cat  = 1/sqrt(2) = 0.707.. (cat state or worse)


// New Decision tree
// Step 1)  Start by filtering eigenvectors down to a  smaller set of "relevant" candidates for
//          doing subspace optimization. Allowing a maximum of 64 candidates keeps ram below 2GB
//          when theta_old.size() == 4096. This means that we filter out
//              * candidates outside of the energy window,
//              * candidates with little or no overlap to the current state.
//          Compute subspace_error_filtered = 1 - Σ_i |<candidate_i|theta_old>|^2
//          If subspace_error_filtered > subspace_error_threshold, set optSpace = DIRECT
//          Else, set optSpace = SUBSPACE.
//
// Step 2)  Find the best overlapping state among the relevant candidates.
//
// Step 3)  We can now make different decisions based on the overlap.
//          A)  If best_overlap_idx == -1
//              No state is in energy window -> discard! Return old theta_old.
//          B)  If overlap_high <= best_overlap.
//              This can happen if the environments have been modified just slightly since the last time considered
//              these sites, but the signal is still clear -- we are still targeting the same state.
//              However we can't be sure that the contributions from nearby states is just noise. Instead of just
//              keeping the state we should optimize its variance. This is important in the later stages when variance
//              is low and we don't want to ruin those last decimals.
//              We just need to decide which initial guess to use.
//                  B1) If best_overlap_variance <= theta_variance: set theta_initial = best_overlap_theta.
//                  B2) Else, set theta_initial = theta_old.
//          LC)  If overlap_cat <= best_overlap and best_overlap < overlap_high
//              This can happen for one reasons:
//                  1) There are a few candidate states with significant overlap (superposition)
//              It's clear that we need to optimize, but we have to think carefully about the initial guess.
//              Right now it makes sense to always choose best overlap theta, since that forces the algorithm to
//              choose a particular state and not get stuck in superposition. Choosing the old theta may not always
//              may just entrench the algorithm into a local minima.
//          D)  If 0 <= best_overlap and best_overlap < overlap_cat
//              This can happen for three reasons, most often early in the simulation.
//                  1) There are several candidate states with significant overlap (superposition)
//                  2) The highest overlapping states were outside of the energy window, leaving just these candidates.
//                  3) The energy targeting of states has failed for some reason, perhaps the spectrum is particularly dense_lu.
//              In any case, it is clear we are lost Hilbert space.
//              Also, the subspace_error is no longer a good measure of how useful the subspace is to us, since it's only
//              measuring how well the old state can be described, but the old state is likely very different from what
//              we're looking for.
//              So to address all three cases, do DIRECT optimization with best_overlap_theta as initial guess.
//
// In particular, notice that we never use the candidate that happens to have the best variance.











// Decision tree:
// Step 1)  Start by filtering eigenvectors down to a  smaller set of "relevant" candidates for
//          doing subspace optimization. Allowing a maximum of 64 candidates keeps ram below 2GB
//          when theta_old.size() == 4096. This means that we filter out
//              * candidates outside of the energy window,
//              * candidates with little or no overlap to the current state.
//              * TODO: Filter states which have high variance?
//          Compute subspace_error_filtered = 1 - Σ_i |<candidate_i|theta_old>|^2
// Step 2)  Find the best overlapping state among the relevant candidates.
// Step 3)  We can now make decisions A-F based on the overlaps.
//          A)  If best_overlap_idx == -1
//              No state is in energy window -> discard! Return old theta_old.
//          B)  If best_overlap >= overlap_high.
//              NOTE: When variance is low we need to be more careful about defining overlap_high.
//              A good estimate may be overlap_high = 1 - variance.
//              This happens when the environments haven't changed and we basically just found the old theta_old among the new eigenvectors,
//              and other eigenvectors just contribute to negligible noise. We could essentially just go ahead and keep it, but sometimes
//              when variance is low we don't want to ruin those last decimals.
//                  B1) If  best_overlap_variance <= theta_old_variance: keep
//
//          LC) If overlap_good <= best_overlap < overlap_high, do variance optimization.
//               This can happen if the environments have been modified slightly since the last time we visited this site,
//               but the signal is still clear -- we are still targeting the same state. However we can't be sure that
//               the contributions from nearby states is just noise anymore.
//               First, set theta_initial = theta_best_overlap_candidate
//               //////First, check the variance ONLY of the best overlapping relevant candidate.
//               //////If the candidate has lower variance than the current one, set theta_initial = theta_best_overlap_candidate.
//               TODO: Which makes the most sense in the two options above?
//               The best course of action now is to:
//                   C1) If the subspace error is low enough, do subspace optimization, with initial guess theta_initial.
//                   C2) Else, send theta_initial as a starting guess for DIRECT optimization.
//                   TODO) Think about what theta_initial is supposed to be. Either it can be the eigenvector with best overlap, or just the old theta_old.
//          D) If overlap_ok < best_overlap < overlap_good
//               This happens if the environments have changed some more since the last time we visited this site,
//               for instance when some other site got optimized a lot.
//               The signal is less clear, but we are probably still targeting the same state.
//               This time we need to be more careful though.
//               First, check the variance of ALL relevant candidates.
//               If a candidate theta_j has lower variance than the current one, set theta_initial = theta_j.
//               Now:
//                   D1) If the subspace quality is good enough, do subspace optimization with initial guess theta_initial
//                   D2) Else, send theta_initial as a starting guess for DIRECT optimization
//          E) If overlap_low < best_overlap < overlap_ok
//               This happens if the environments have changed a lot since the last time we visited this site,
//               which is usually the case early in the simulation.
//               The signal is not clear anymore, in fact there are many candidates with significant overlap to the old theta_old.
//               The subspace_error is probably not a good measure anymore, since we're not trying to find a new theta which
//               is only a fine tuning away from of the old one: we would just get stuck in a local minima far away.
//               First, check the variance of ALL relevant candidates.
//                   E1) If any candidate state has better variance than the current one, send it as a starting guess for DIRECT optimization.
//                   E2) Else, send the best overlapping state as a starting guess for DIRECT optimization.
//          F) If best_overlap < overlap_low
//               Mayday! We are lost in Hilbert space!
//               Send the old theta as a starting guess for DIRECT optimization, and brace for impact.















//
//
//
//
//
//    auto [best_overlap,best_overlap_idx] = get_best_overlap_in_window(overlaps, eigvals_per_site_unreduced, sim_status.energy_lbound, sim_status.energy_ubound);
//
//    if (best_overlap_idx < 0){
//        //Option A
//        tools::log->info("Went for option A");
//        tools::log->debug("No overlapping states in energy range. Returning old theta");
//        state.tag_active_sites_have_been_updated(false);
//        return theta;
//    }
//    else if(best_overlap > settings::precision::overlap_high){
//        //Option B
//        tools::log->info("Went for option B");
//        tools::log->debug("... Overlap of candidate {} is great: {} . Keeping it.", best_overlap_idx, best_overlap);
//        state.tag_active_sites_have_been_updated(true);
//        state.clear_measurements();
//        return Textra::Matrix_to_Tensor(eigvecs.col(best_overlap_idx), state.active_dimensions());
//    }
//    else if(best_overlap > settings::precision::overlap_good and best_overlap < settings::precision::overlap_high ){
//        if(subspace_error < subspace_error_threshold){
//            //Option C1
//            tools::log->info("Went for option C1");
//            tools::log->debug("... We can try subspace optimization anyway then");
//            std::tie(eigvecs,eigvals,subspace_error) = filter_states(eigvecs, eigvals, overlaps, subspace_error_threshold, 64);
//            eigvals_per_site_unreduced = (eigvals.array() + state.get_energy_reduced())/state.get_length(); // Remove energy reduction for energy window comparisons
//        }else{
//            //Option C2
//            tools::log->info("Went for option C2");
//            tools::log->debug("... Switching to DIRECT mode");
//            return ceres_direct_optimization(state, Textra::Matrix_to_Tensor(eigvecs.col(best_overlap_idx), state.active_dimensions()) ,sim_status, optType);
//        }
//
//    }
//    else if(best_overlap > settings::precision::lowOverlap and best_overlap < settings::precision::overlap_good ){
//        if(subspace_error < subspace_error_threshold){
//            //Option D1
//            tools::log->info("Went for option D1");
//            tools::log->debug("... We can try subspace optimization anyway then");
//            std::tie(eigvecs,eigvals,subspace_error) = filter_states(eigvecs, eigvals, overlaps, subspace_error_threshold, 64);
//            eigvals_per_site_unreduced = (eigvals.array() + state.get_energy_reduced())/state.get_length(); // Remove energy reduction for energy window comparisons
//        }else{
//            //Option D2
//            tools::log->info("Went for option D2");
//            tools::log->debug("... Overlap of candidate {} is low: {} . Keeping it.", best_overlap_idx, best_overlap);
//            state.tag_active_sites_have_been_updated(true);
//            state.clear_measurements();
//            return Textra::Matrix_to_Tensor(eigvecs.col(best_overlap_idx), state.active_dimensions());
//        }
//    }
//    else if(best_overlap < settings::precision::lowOverlap and best_overlap > settings::precision::badOverlap ){
//        if(subspace_error < subspace_error_threshold){
//            //Option E1
//            tools::log->info("Went for option E1");
//            tools::log->debug("... We can try subspace optimization anyway then");
//            std::tie(eigvecs,eigvals,subspace_error) = filter_states(eigvecs, eigvals, overlaps, subspace_error_threshold, 64);
//            eigvals_per_site_unreduced = (eigvals.array() + state.get_energy_reduced())/state.get_length(); // Remove energy reduction for energy window comparisons
//        }else{
//            //Option E2
//            tools::log->info("Went for option E2");
//            tools::log->debug("... Overlap of candidate {} is low: {} . Keeping it.", best_overlap_idx, best_overlap);
//            state.tag_active_sites_have_been_updated(true);
//            state.clear_measurements();
//            return Textra::Matrix_to_Tensor(eigvecs.col(best_overlap_idx), state.active_dimensions());
//        }
//    }else{
//        throw std::runtime_error("Nothing matched");
//    }





//    if(best_overlap > settings::precision::overlap_good and best_overlap < settings::precision::overlap_high ){
//        //Option C1
//        tools::log->info("Went for option C1");
//        tools::log->debug("... Candidate {} has fair overlap {} and variance (log10): {}", best_variance_idx, best_overlap, std::log10(best_overlap_variance));
//        state.tag_active_sites_have_been_updated(true);
//        state.clear_measurements();
//        return best_overlap_theta;
//    }
//    if(best_overlap > high_overlap and subspace_error < 100*subspace_error_threshold ){
//        //Option D1
//        tools::log->info("Went for option D1");
//        tools::log->debug("... We can try subspace anyway then", best_variance_idx, best_overlap, std::log10(best_overlap_variance));
//        std::tie(eigvecs,eigvals,subspace_error) = filter_states(eigvecs, eigvals, overlaps, subspace_error_threshold, 64);
//        eigvals_per_site_unreduced = (eigvals.array() + state.get_energy_reduced())/state.get_length(); // Remove energy reduction for energy window comparisons
//        break;
//    }


//
//
//    switch (optMode){
//        case OptMode::OVERLAP:
//        {
//            auto [best_overlap,best_overlap_idx] = get_best_overlap_in_window(overlaps, eigvals_per_site_unreduced, sim_status.energy_lbound, sim_status.energy_ubound);
//            if (best_overlap_idx < 0){
//                tools::log->debug("No overlapping states in energy range. Returning old theta");
//                state.tag_active_sites_have_been_updated(false);
//                return theta;
//            }
//            if(best_overlap < 0.1){
//                tools::log->debug("Overlap of state {} is too low: {}. Checking for candidates with lower variance", best_overlap_idx, best_overlap);
//                double old_variance = tools::finite::measure::energy_variance_per_site(state);
//                auto [best_variance, best_variance_idx] = get_best_variance_in_window(state, eigvecs, eigvals_per_site_unreduced, sim_status.energy_lbound, sim_status.energy_ubound);
//                if (best_variance_idx < 0 or overlaps(best_variance_idx) < 0.01){
//                    tools::log->debug("No better variance states (with sufficient overlap > 0.01) found in energy range. Returning old theta");
//                    state.tag_active_sites_have_been_updated(false);
//                    return theta;
//                }
//                if(best_variance < old_variance){
//                    tools::log->debug("... Eigenstate {} had better (log10) variance: {} < {}. Energy: {}, overlap: {}.",
//                                      best_variance_idx, std::log10(best_variance), std::log10(old_variance), eigvals_per_site_unreduced(best_variance_idx), overlaps(best_variance_idx));
//                    state.tag_active_sites_have_been_updated(true);
//                    state.clear_measurements();
//                    return Textra::Matrix_to_Tensor(eigvecs.col(best_variance_idx), state.active_dimensions());
//                }else{
//                    tools::log->debug("... No found state had good enough overlap or varaince, returning old theta");
//                    state.tag_active_sites_have_been_updated(false);
//                    return theta;
////                    tools::log->debug("... No eigenstate was better, keeping badly overlapping state");
////                    state.tag_active_sites_have_been_updated(true);
////                    state.clear_measurements();
////                    return Textra::Matrix_to_Tensor(eigvecs.col(best_overlap_idx), state.active_dimensions());
//
//                }
//            }else{
//                tools::log->debug("Candidate theta {} has good overlap {}", best_overlap_idx, best_overlap);
//                auto   new_theta     = Textra::Matrix_to_Tensor(eigvecs.col(best_overlap_idx), state.active_dimensions());
//                double old_variance  = tools::finite::measure::energy_variance_per_site(state);
//                double new_variance  = tools::finite::measure::energy_variance_per_site(state,new_theta);
//                // Check that the new state is smaller than at least twice the old one
//                if (new_variance <= 2*old_variance){
//                    tools::log->debug("Kept candidate {} -- it has good enough overlap {} and variance {}", best_overlap_idx, best_overlap, std::log10(new_variance));
//                    state.tag_active_sites_have_been_updated(true);
//                    state.clear_measurements();
//                    return new_theta;
//                }else{
//                    tools::log->debug("The candidate theta has worse variance than before [ idx = {} | overlap = {} | variance = {} ]...", best_overlap_idx, best_overlap, std::log10(new_variance));
//                    tools::log->debug("Looking for a candidate with lower variance...");
//                    double subspace_error;
//                    std::tie(eigvecs,eigvals,subspace_error) = filter_states(eigvecs, eigvals, overlaps, subspace_error_threshold, 64);
//                    eigvals_per_site_unreduced = (eigvals.array() + state.get_energy_reduced())/state.get_length(); // Remove energy reduction for energy window comparisons
//                    auto [best_variance, best_variance_idx] = get_best_variance_in_window(state, eigvecs, eigvals_per_site_unreduced, sim_status.energy_lbound, sim_status.energy_ubound);
//                    if(best_variance < old_variance){
//                        tools::log->debug("... Candidate {} has better variance: {} < {}. Energy: {}, overlap: {}.",
//                                          best_variance_idx, std::log10(best_variance), std::log10(old_variance), eigvals_per_site_unreduced(best_variance_idx), overlaps(best_variance_idx));
//                        state.tag_active_sites_have_been_updated(true);
//                        state.clear_measurements();
//                        return Textra::Matrix_to_Tensor(eigvecs.col(best_variance_idx), state.active_dimensions());
//                    }
//                    else{
//                        tools::log->debug("... No candidate has good enough overlap or variance, returning old theta");
//                        state.tag_active_sites_have_been_updated(false);
//                        return theta;
//                    }
//                }
//            }
//            break;
//        }
//
//
//
//        case OptMode::VARIANCE:
//        {
//            bool preferKeepingLowOverlapCandidate = true;
//            auto  [best_overlap,best_overlap_idx] = get_best_overlap_in_window(overlaps, eigvals_per_site_unreduced, sim_status.energy_lbound, sim_status.energy_ubound);
//
//            double sq_sum_overlap    = overlaps.cwiseAbs2().sum();
//            double subspace_error    = 1.0 - sq_sum_overlap;
//
//
//
//            // LC) If low_overlap < best_overlap < medium_overlap
//            // Ooops, the subspace error is too high. We still have some options. In order of priority:
//            // a) If no state is inside the energy window, discard all and return old theta
//            // b) If any state inside the energy window has lower variance, keep it
//
//            // Now we can do different things depending on the variable preferKeepingLowOverlapCandidate
//            // The reasoning is that if the maximum overlap is too low, perhaps we're at a local minima and
//            // it's not worth it to keep optimizing there. Better then to keep that candidate and escape the minima.
//            // So then, if preferKeepingLowOverlapCandidate == true
//            // c1) If the maximum overlap in energy window is intermediate, say between 0.1 and 0.9, and its variance isn't too bad, then keep it.
//            // d1) If the maximum overlap in energy window is high, say higher than 0.9, and the the subspace error isn't too bad, try subspace optimization anyway
//            // but if preferKeepingLowOverlapCandidate == false
//            // c2) If the maximum overlap in energy window is high enough, say higher than 0.9, and its variance isn't too bad, then keep it
//            // d2) If the maximum overlap in energy window is low, say lower than 0.99, and the the subspace error isn't too bad, try subspace optimization anyway
//            // e) If any state is inside the energy window, but none of the above applies, switch to DIRECT
//
//
//
//
//            if(subspace_error > subspace_error_threshold) {
//
//                tools::log->debug("Subspace error is too high (log10): {} > {}. Deciding what to do...", std::log10(subspace_error), std::log10(subspace_error_threshold));
//                double prev_variance = tools::finite::measure::energy_variance_per_site(state);
//                auto [best_variance, best_variance_idx] = get_best_variance_in_window(state, eigvecs, eigvals_per_site_unreduced, sim_status.energy_lbound, sim_status.energy_ubound);
//                if (best_variance_idx < 0){
//                    // Option A
//                    tools::log->info("Went for option A");
//                    tools::log->debug("... No candidate in energy window, returning old theta");
//                    state.tag_active_sites_have_been_updated(false);
//                    return theta;
//                }
//                tools::log->debug("... Candidate {} has lowest variance (log10): {}", best_variance_idx, std::log10(best_variance));
//                if(best_variance < prev_variance){
//                    // Option B
//                    tools::log->info("Went for option B");
//                    tools::log->debug("... Candidate {} has better variance (log10): {} < {}", best_variance_idx, std::log10(best_variance), std::log10(prev_variance));
//                    state.tag_active_sites_have_been_updated(true);
//                    state.clear_measurements();
//                    return Textra::Matrix_to_Tensor(eigvecs.col(best_variance_idx), state.active_dimensions());
//                }
//
//                auto [best_overlap,best_overlap_idx] = get_best_overlap_in_window(overlaps, eigvals_per_site_unreduced, sim_status.energy_lbound, sim_status.energy_ubound);
//                auto best_overlap_theta = Textra::Matrix_to_Tensor(eigvecs.col(best_overlap_idx), state.active_dimensions());
//                double best_overlap_variance = tools::finite::measure::energy_variance_per_site(state, best_overlap_theta);
//                tools::log->debug("... Candidate {} has highest overlap: {} and variance(log10): {}", best_overlap_idx, best_overlap ,std::log10(best_overlap_variance));
//
//                if(preferKeepingLowOverlapCandidate){
//                    if(best_overlap > low_overlap and best_overlap < high_overlap ){
//                        //Option C1
//                        tools::log->info("Went for option C1");
//                        tools::log->debug("... Candidate {} has fair overlap {} and variance (log10): {}", best_variance_idx, best_overlap, std::log10(best_overlap_variance));
//                        state.tag_active_sites_have_been_updated(true);
//                        state.clear_measurements();
//                        return best_overlap_theta;
//                    }
//                    if(best_overlap > high_overlap and subspace_error < 100*subspace_error_threshold ){
//                        //Option D1
//                        tools::log->info("Went for option D1");
//                        tools::log->debug("... We can try subspace anyway then", best_variance_idx, best_overlap, std::log10(best_overlap_variance));
//                        std::tie(eigvecs,eigvals,subspace_error) = filter_states(eigvecs, eigvals, overlaps, subspace_error_threshold, 64);
//                        eigvals_per_site_unreduced = (eigvals.array() + state.get_energy_reduced())/state.get_length(); // Remove energy reduction for energy window comparisons
//                        break;
//                    }
//                }else{
//                    if(best_overlap > high_overlap and best_overlap_variance < 100.0 * prev_variance ){
//                        //Option C2
//                        tools::log->info("Went for option C2");
//                        tools::log->debug("... Candidate {} has fair overlap {} and variance (log10): {}", best_variance_idx, best_overlap, std::log10(best_overlap_variance));
//                        state.tag_active_sites_have_been_updated(true);
//                        state.clear_measurements();
//                        return best_overlap_theta;
//                    }
//                    if(best_overlap < high_overlap and subspace_error < 100*subspace_error_threshold ){
//                        //Option D2
//                        tools::log->info("Went for option D2");
//                        tools::log->debug("... We can try subspace anyway then", best_variance_idx, best_overlap, std::log10(best_overlap_variance));
//                        std::tie(eigvecs,eigvals,subspace_error) = filter_states(eigvecs, eigvals, overlaps, subspace_error_threshold, 64);
//                        eigvals_per_site_unreduced = (eigvals.array() + state.get_energy_reduced())/state.get_length(); // Remove energy reduction for energy window comparisons
//                        break;
//                    }
//                }
//
//
//                // Option E
//                tools::log->info("Went for option E");
//                tools::log->debug("... Switching to DIRECT mode");
//                return ceres_direct_optimization(state, sim_status, optType);
//                //////////////////////////////
////                tools::log->debug("Subspace error is too high (log10): {} > {}. Deciding what to do...", std::log10(subspace_error), std::log10(subspace_error_threshold));
////                double prev_variance = tools::finite::measure::energy_variance_per_site(state);
////                auto [best_variance, best_variance_idx] = get_best_state_in_window(state,eigvecs,eigvals_per_site_unreduced,sim_status.energy_lbound,sim_status.energy_ubound);
////                if (best_variance_idx < 0){
////                    tools::log->debug("Returning old theta");
////                    state.tag_active_sites_have_been_updated(false);
////                    return theta;
////                }
////                else if(best_variance < prev_variance){
////                    tools::log->debug("... Eigenstate {} has better (log10) variance: {} < {}", best_variance_idx, std::log10(best_variance), std::log10(prev_variance));
////                    state.tag_active_sites_have_been_updated(true);
////                    state.clear_measurements();
////                    return Textra::Matrix_to_Tensor(eigvecs.col(best_variance_idx), state.active_dimensions());
////                }
//            }else{
//                std::tie(eigvecs,eigvals,subspace_error) = filter_states(eigvecs, eigvals, overlaps, subspace_error_threshold, 64);
//                eigvals_per_site_unreduced = (eigvals.array() + state.get_energy_reduced())/state.get_length(); // Remove energy reduction for energy window comparisons
//            }
//        }
//        break;
//    }















//
//        double prev_variance = tools::finite::measure::energy_variance_per_site(state);
//        auto [best_variance, idx_variance] = get_best_variance_in_window(state,eigvecs,eigvals_per_site_unreduced,sim_status.energy_lbound,sim_status.energy_ubound);
//        if (idx_variance < 0){
//            tools::log->debug("Returning old theta");
//            state.tag_active_sites_have_been_updated(false);
//            return theta;
//        }
//        else if(best_variance < prev_variance){
//            tools::log->debug("... Eigenstate {} has better (log10) variance: {} < {}",idx_variance, std::log10(best_variance), std::log10(prev_variance));
//            state.tag_active_sites_have_been_updated(true);
//            return Textra::Matrix_to_Tensor(eigvecs.col(idx_variance), state.active_dimensions());
//        }else{
//            tools::log->debug("... discarding subspace and switching to direct mode");
//            return ceres_direct_optimization(state, sim_status, optType);
//        }
//    }
