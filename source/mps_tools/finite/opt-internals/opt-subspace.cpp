//
// Created by david on 2019-03-18.
//
#include <Eigen/Core>
#include <iostream>
#include <algorithms/class_simulation_state.h>
#include <general/class_eigsolver.h>
#include <general/arpack_extra/matrix_product_stl.h>
#include <general/arpack_extra/matrix_product_sparse.h>
#include <spdlog/spdlog.h>
#include <mps_tools/finite/opt.h>
#include <mps_state/class_finite_chain_state.h>
#include <LBFGS.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <general/nmspc_random_numbers.h>

using namespace mpstools::finite::opt;
using namespace mpstools::finite::opt::internals;

template<typename T> using MatrixType = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;




template<typename Scalar>
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
find_subspace_full(const MatrixType<Scalar> & H_local, Eigen::Tensor<std::complex<double>,3> &theta, std::vector<reports::eig_tuple> &eig_log){
    mpstools::log->trace("Finding subspace -- full");
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
    mpstools::log->debug("Finished eigensolver -- condition: Full diagonalization");
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


std::vector<int> mpstools::finite::opt::internals::generate_size_list(size_t shape){
    std::vector<int> nev_list;
    int max_nev = std::max(std::min(8,(int)shape),(int)shape/8);
    max_nev = std::min(max_nev,256);

    int min_nev = std::min(std::min(8,(int)shape),max_nev);

    int tmp_nev = min_nev;
    while (tmp_nev <= max_nev){
        nev_list.push_back(tmp_nev);
        tmp_nev *= 4;
    }
    return nev_list;
}

template<typename Scalar>
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
find_subspace_part(const MatrixType<Scalar> & H_local, Eigen::Tensor<std::complex<double>,3> &theta, double energy_target, std::vector<reports::eig_tuple> &eig_log){
    using namespace eigutils::eigSetting;
    mpstools::log->trace("Finding subspace -- partial");


    t_eig->tic();
    // You need to copy the data into StlMatrixProduct, because the PartialPivLU will overwrite the data in H_local otherwise.
    StlMatrixProduct<Scalar> hamiltonian(H_local.data(),H_local.rows(),Form::SYMMETRIC,Side::R, true);
    hamiltonian.set_shift(energy_target);
    hamiltonian.FactorOP();
    double t_lu = hamiltonian.t_factorOp.get_last_time_interval();
    t_eig->toc();

    double prec                       = settings::precision::VarConvergenceThreshold;
    double max_overlap_threshold      = 1 - prec; //1.0/std::sqrt(2); //Slightly less than 1/sqrt(2), in case that the choice is between cat states.
    double subspace_quality_threshold = prec;

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
        if(max_overlap    > 1.0 + 1e-10) throw std::runtime_error("max_overlap larger than one : "  + std::to_string(max_overlap));
        if(sq_sum_overlap > 1.0 + 1e-10) throw std::runtime_error("eps larger than one : "          + std::to_string(sq_sum_overlap));
        if(min_overlap    < 0.0)         throw std::runtime_error("min_overlap smaller than zero: " + std::to_string(min_overlap));
        if(max_overlap >= max_overlap_threshold )         {reason = "overlap is good enough"; break;}
        if(subspace_quality < subspace_quality_threshold) {reason = "subspace quality is good enough"; break;}
    }
    mpstools::log->debug("Finished partial eigensolver -- condition: {}",reason);
    return std::make_tuple(eigvecs,eigvals);
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
find_subspace(const class_finite_chain_state & state, const class_simulation_state & sim_state, OptMode optMode, OptSpace optSpace){
    mpstools::log->trace("Finding subspace");

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
    t_ham->toc();
    auto theta = state.get_multitheta();


    Eigen::MatrixXcd eigvecs;
    Eigen::VectorXd  eigvals;
    std::vector<reports::eig_tuple> eig_log;
    double energy_target = state.get_length() * sim_state.energy_now;

    // If theta is small enough, do FULL even if it asked for PARTIAL.
    if   (theta.size() <= settings::precision::MaxSizeFullDiag) optSpace = OptSpace::FULL;
    switch (optSpace){
        case OptSpace::FULL    : std::tie(eigvecs,eigvals) = find_subspace_full(H_local,theta,eig_log); break;
        case OptSpace::PARTIAL : std::tie(eigvecs,eigvals) = find_subspace_part(H_local,theta,energy_target,eig_log); break;
        default:
            std::stringstream optspace_sstr;
            optspace_sstr << optSpace;
            throw std::runtime_error("Wrong OptSpace for subspace optimization. Expected FULL or PARTIAL. Got: " + optspace_sstr.str() );
    }
    reports::print_report(eig_log);

    if constexpr(std::is_same<Scalar,double>::value){
        Textra::subtract_phase(eigvecs);
        std::cout << "truncating imag of eigvecs, sum: " << eigvecs.imag().cwiseAbs().sum() << std::endl;
        eigvecs = eigvecs.real();
    }

    eigvecs.colwise().normalize();

    switch (optMode){
        case OptMode::VARIANCE: return std::make_tuple(eigvecs, eigvals);
        case OptMode::OVERLAP :
            Eigen::Map<Eigen::VectorXcd> theta_vec(theta.data(),theta.size());
            Eigen::VectorXd overlaps  = (theta_vec.adjoint() * eigvecs).cwiseAbs().real();
            int idx = idx_best_overlap_in_window(overlaps,eigvals,sim_state.energy_lbound,sim_state.energy_ubound);
            return std::make_tuple(eigvecs.col(idx),eigvals.row(idx));
    }
}



std::tuple<Eigen::Tensor<std::complex<double>,3>, double>
mpstools::finite::opt::internals::subspace_optimization(const class_finite_chain_state & state, const class_simulation_state & sim_state, OptMode optMode, OptSpace optSpace, OptType optType){
    mpstools::log->trace("Optimizing in SUBSPACE mode");
    using Scalar = std::complex<double>;
    using namespace eigutils::eigSetting;

    double chain_length    = state.get_length();
    auto theta             = state.get_multitheta();
    auto theta_old         = Eigen::Map<const Eigen::VectorXcd>  (theta.data(),theta.size());
    Eigen::MatrixXcd eigvecs;
    Eigen::VectorXd  eigvals;
    switch(optType){
        case OptType::CPLX:     std::tie (eigvecs,eigvals)  = find_subspace<Scalar>(state,sim_state, optMode,optSpace); break;
        case OptType::REAL:     std::tie (eigvecs,eigvals)  = find_subspace<double>(state,sim_state, optMode,optSpace); break;
    }

    if (optMode == OptMode::OVERLAP){
        return std::make_tuple(
                Textra::Matrix_to_Tensor(eigvecs, state.active_dimensions()),
                eigvals(0)/state.get_length()
        );
    }else if (optMode == OptMode::VARIANCE){
        Eigen::VectorXd overlaps = (theta_old.adjoint() * eigvecs).cwiseAbs().real();
        double sq_sum_overlap    = overlaps.cwiseAbs2().sum();
        double subspace_quality  = 1.0 - sq_sum_overlap;
        if(subspace_quality >= settings::precision::VarConvergenceThreshold) {
            return direct_optimization(state, sim_state, optType);
        }
    }


    Eigen::VectorXcd theta_new;
    double overlap_new  = 0;
    double energy_new,variance_new,norm;
    //Should really use theta_start as the projection towards the previous theta, not best overlapping!
    // Note that alpha_i = <theta_old | theta_new_i> is not supposed to be squared! The overlap
    // Between theta_start and theta_old should be
    Eigen::VectorXcd theta_start      = (eigvecs.adjoint()  * theta_old)  ;

    std::vector<reports::subspc_opt_tuple> opt_log;


    t_opt->tic();
    state.unset_measurements();
    double energy_0   = mpstools::finite::measure::energy_per_site(state);
    double variance_0 = mpstools::finite::measure::energy_variance_per_site(state);
    t_opt->toc();
    Eigen::VectorXcd theta_0 = (eigvecs * theta_start.asDiagonal()).rowwise().sum().normalized();

    int iter_0 = 0;
    double overlap_0 = std::abs(theta_old.dot(theta_0));

    opt_log.emplace_back("Initial",theta.size(), energy_0, std::log10(variance_0), overlap_0, iter_0,0, t_opt->get_last_time_interval());



//    Eigen::VectorXcd theta_old_check = theta_old;
//    Scalar angle = 3.0;
//    Scalar exp_random_phase = std::exp(Scalar(0.0,-1.0) * angle);
//    theta_old_check *= exp_random_phase;
//    for (int i= 0; i < eigvecs.cols() ; i++){
//        Scalar angle = 6*rn::uniform_double_1();
//        Scalar exp_random_phase = std::exp(Scalar(0.0,-1.0) * angle);
//        eigvecs.col(i) *= exp_random_phase;
//    }


//    auto check_num = (eigvecs.adjoint() * eigvecs).diagonal().sum();
//    Eigen::VectorXcd theta_start_check     =  (eigvecs.adjoint()  * theta_old_check) ;
//    theta_start = (eigvecs.adjoint()  * theta_old_check);

//    Eigen::MatrixXcd H  = eigvecs.adjoint() * state.get_H_local_matrix() * eigvecs;
//    Eigen::MatrixXcd H2 = eigvecs.adjoint() * state.get_H_local_sq_matrix() * eigvecs;
//
//    std::cout << "check_num              : " << check_num << std::endl;
//    std::cout << "Norm theta_start_check : " << theta_start_check.norm() << std::endl;
//    std::cout << "Size theta_start_check : " << theta_start_check.rows() << std::endl;
//    std::cout << "Size eigvals           : " << eigvals.rows()  << std::endl;
//    std::cout << "Size eigvecs           : " << eigvecs.rows() << " x " << eigvecs.cols() << std::endl;
//    std::cout << "Size H2                : " << H2.rows() << " x " << H2.cols() << std::endl;
//
//    Scalar vHv   = theta_start_check.adjoint() * eigvals.asDiagonal()               * theta_start_check;
//    Scalar vH2v  = theta_start_check.adjoint() * H2                                 * theta_start_check;
//
//    Scalar energy_check   = vHv;
//    Scalar variance_check = vH2v - vHv * vHv ;
//
//    Eigen::VectorXcd vH   = theta_start_check.adjoint() * eigvals.asDiagonal();
////    Scalar vHv_   =  (vH.transpose() *  theta_start_check).eval()(0);
//    Scalar vHv_ =  (vH.transpose() *  theta_start_check)(0);
//    Eigen::VectorXcd vH2  = theta_start_check.adjoint() * H2;
//    Scalar vH2v_ = (vH2.transpose()*theta_start_check)(0);
//    Scalar energy_check2   = vHv_;
//    Scalar variance_check2 = vH2v_ - vHv_ * vHv_ ;
//
//    Eigen::VectorXcd theta_check = (eigvecs * theta_start_check.asDiagonal()).rowwise().sum().normalized();
//
//    Scalar overlap_check            = theta_old_check.adjoint() * theta_check ;
//
//    double energy_init   = mpstools::common::measure::energy         (state,Textra::Matrix_to_Tensor(theta_check, state.dimensions()));
//    double variance_init = mpstools::common::measure::energy_variance(state,Textra::Matrix_to_Tensor(theta_check, state.dimensions()),energy_init);
//

//    opt_log.emplace_back("Sanity check",theta.size(), vHv.real()/chain_length, std::log10(variance_check.real()/ chain_length), std::abs(overlap_check), 0,0, 0);


    t_opt->tic();
    mpstools::log->trace("Running LBFGS");
    using namespace LBFGSpp;
    using namespace mpstools::finite::opt::internals;
    double fx;
    int niter,counter;
    if (sim_state.variance_mpo_has_converged or variance_0 < settings::precision::VarConvergenceThreshold){params.max_iterations = 10; params.max_linesearch = 20;}
    else{ params = get_params();}
    LBFGSpp::LBFGSSolver<double> solver(params);
    switch (optType){
        case OptType::CPLX:{
            mpstools::finite::opt::internals::subspace_functor <Scalar>  functor (state,sim_state,eigvecs,eigvals);
            Eigen::VectorXd  theta_start_cast = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double*> (theta_start.data()), 2*theta_start.size());
            niter = solver.minimize(functor, theta_start_cast, fx);
            theta_start = Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<Scalar*> (theta_start_cast.data()), theta_start_cast.size()/2).normalized();
            counter      = functor.get_count();
            norm         = functor.get_norm();
            energy_new   = functor.get_energy() / chain_length;
            variance_new = functor.get_variance()/chain_length;
            theta_new    = (eigvecs * theta_start.asDiagonal()).rowwise().sum().normalized();
            break;
        }
        case OptType::REAL:{
            mpstools::finite::opt::internals::subspace_functor <double> functor (state,sim_state,eigvecs.real(),eigvals);
            Eigen::VectorXd  theta_start_cast = theta_start.real();
            niter = solver.minimize(functor, theta_start_cast, fx);
            theta_start = theta_start_cast.normalized().cast<Scalar>();
            counter      = functor.get_count();
            norm         = functor.get_norm();
            energy_new   = functor.get_energy() / chain_length;
            variance_new = functor.get_variance()/chain_length;
            theta_new    = (eigvecs.real() * theta_start.real().asDiagonal()).rowwise().sum().normalized();

            break;
        }
    }
    t_opt->toc();

    overlap_new = (theta_old.adjoint() * theta_new).cwiseAbs().sum();
    opt_log.emplace_back("LBFGS++",theta.size(), energy_new, std::log10(variance_new), overlap_new, niter,counter, t_opt->get_last_time_interval());
    mpstools::log->trace("Finished LBFGS");

    reports::print_report(opt_log);
    state.unset_measurements();
//    double energy_tmp   = mpstools::finite::measure::energy(state,Textra::Matrix_to_Tensor(theta_new, state.active_dimensions()));
//    double variance_tmp = mpstools::finite::measure::energy_variance(state,Textra::Matrix_to_Tensor(theta_new, state.active_dimensions()),energy_tmp);

//    std::cout << "Energy   check2           : " << energy_check2/chain_length    << std::endl;
//    std::cout << "Variance check2 (complex) : " << std::log10(variance_check2/chain_length)         << std::endl;
//    std::cout << "Energy   init             : " << energy_init/chain_length    << std::endl;
//    std::cout << "Energy   check (complex)  : " << energy_check/chain_length   << std::endl;
//    std::cout << "Overlap  check (complex)  : " << overlap_check  << std::endl;
//    std::cout << "Variance check (complex)  : " << std::log10(variance_check/chain_length)         << std::endl;
//    std::cout << "Variance init             : " << std::log10(variance_init /chain_length)         << "\t \t norm = " << theta_check.norm() << std::endl;
//    std::cout << "Variance check            : " << std::log10(variance_check.real()/chain_length)  << "\t \t norm = " << theta_check.norm() << std::endl;
//    std::cout << std::setprecision(16) << std::fixed;
//    std::cout << "Variance LBFGS            : " << std::log10(variance_new)                        << "\t \t norm = " << norm << std::endl;
//    std::cout << "Variance after            : " << std::log10(variance_tmp/chain_length)           << "\t \t norm = " << theta_new.norm()   << std::endl;

//    std::cout << "imag sum theta_new        : " << theta_new.imag().cwiseAbs().sum() << std::endl;
//    std::cout << "imag sum theta_old        : " << theta_old.imag().cwiseAbs().sum() << std::endl;
//    std::cout << "imag sum eigvecs          : " << eigvecs.imag().cwiseAbs().sum() << std::endl;


    if (variance_new < variance_0){
        mpstools::log->info("Returning new theta");
        return  std::make_tuple(Textra::Matrix_to_Tensor(theta_new, state.active_dimensions()), energy_new);

    }else{
        mpstools::log->info("Returning old theta");
        return  std::make_tuple(theta, energy_0);
    }


}






