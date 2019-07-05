//
// Created by david on 2019-05-25.
//


#include <spdlog/spdlog.h>
#include <LBFGS.h>
#include <mps_tools/finite/opt.h>
#include <general/class_tic_toc.h>
#include <mps_state/class_finite_chain_state.h>
#include <mps_state/class_environment.h>
#include <algorithms/class_simulation_state.h>
#include <general/nmspc_random_numbers.h>
#include <sim_parameters/nmspc_sim_settings.h>

#include <variant>

template<typename Scalar, auto rank>
void warn_if_has_imaginary_part(const Eigen::Tensor<Scalar,rank> &tensor, double threshold = 1e-14) {
    Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>> vector (tensor.data(),tensor.size());
    if constexpr (std::is_same<Scalar, std::complex<double>>::value){
        auto imagSum = vector.imag().cwiseAbs().sum();
        if (imagSum > threshold){
            mpstools::log->warn("Has imaginary part. Sum: " + std::to_string(imagSum));
        }
    }
}


std::tuple<Eigen::Tensor<std::complex<double>,3>, double>
mpstools::finite::opt::internals::direct_optimization(const class_finite_chain_state & state, const class_simulation_state & sim_state, OptType optType){
    mpstools::log->trace("Optimizing in DIRECT mode");
    using Scalar = std::complex<double>;
    auto theta = state.get_multitheta();


    t_opt->tic();
    double chain_length    = state.get_length();
    double energy_new,variance_new,overlap_new;
    Eigen::VectorXcd theta_start  = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size());

    double energy_0   = mpstools::finite::measure::energy_per_site(state);
    double variance_0 = mpstools::finite::measure::energy_variance_per_site(state);
    int iter_0 = 0;

    t_opt->toc();
    std::vector<reports::direct_opt_tuple> opt_log;
    opt_log.emplace_back("Initial",theta.size(), energy_0, std::log10(variance_0), 1.0, iter_0 ,0,t_opt->get_last_time_interval());
    t_opt->tic();
    using namespace LBFGSpp;
    using namespace mpstools::finite::opt::internals;
    double fx;
    int niter,counter;
//    if (sim_state.variance_mpo_has_converged or variance_0 < settings::precision::VarConvergenceThreshold){params.max_iterations = 10; params.max_linesearch = 20;}
//    else{ params = get_params();}
    LBFGSpp::LBFGSSolver<double> solver(params);

    switch (optType){
        case OptType::CPLX:{
            mpstools::finite::opt::internals::direct_functor <Scalar>  functor (state, sim_state);
            Eigen::VectorXd  theta_start_cast = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double*> (theta_start.data()), 2*theta_start.size());
            mpstools::log->trace("Running LBFGS");
            niter = solver.minimize(functor, theta_start_cast, fx);
            theta_start  = Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<Scalar*> (theta_start_cast.data()), theta_start_cast.size()/2).normalized();
            counter      = functor.get_count();
            energy_new   = functor.get_energy() / chain_length;
            variance_new = functor.get_variance()/chain_length;
            break;
        }
        case OptType::REAL:{
            mpstools::finite::opt::internals::direct_functor <double> functor (state, sim_state);
            Eigen::VectorXd  theta_start_cast = theta_start.real();
            mpstools::log->trace("Running LBFGS");
            niter = solver.minimize(functor, theta_start_cast, fx);
            theta_start  = theta_start_cast.normalized().cast<Scalar>();
            counter      = functor.get_count();
            energy_new   = functor.get_energy() / chain_length;
            variance_new = functor.get_variance()/chain_length;
            break;
        }
    }
    t_opt->toc();


    auto theta_old = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size());


    overlap_new  = std::abs(theta_old.dot(theta_start));
    opt_log.emplace_back("LBFGS++",theta.size(), energy_new, std::log10(variance_new), overlap_new, niter,counter, t_opt->get_last_time_interval());
    mpstools::log->trace("Finished LBFGS");

    reports::print_report(opt_log);
    reports::print_report(std::make_tuple(
            mpstools::finite::opt::internals::t_vH2v->get_measured_time(),
            mpstools::finite::opt::internals::t_vHv->get_measured_time(),
            mpstools::finite::opt::internals::t_vH2->get_measured_time(),
            mpstools::finite::opt::internals::t_vH->get_measured_time(),
            mpstools::finite::opt::internals::t_op->get_measured_time()
            ));

    state.unset_measurements();
    if (variance_new < variance_0){
        mpstools::log->info("Returning new theta");
        return  std::make_tuple(Textra::Matrix_to_Tensor(theta_start, state.active_dimensions()), energy_new);

    }else{
        mpstools::log->info("Returning old theta");
        return  std::make_tuple(theta, energy_0);
    }

//    bool outside_of_window = energy_new < sim_state.energy_lbound or energy_new > sim_state.energy_ubound;
//    if (outside_of_window){
//        if(rn::uniform_double_1() < 0.1){
//            return  std::make_tuple(theta, energy_0);
//        }else{
//            std::cout << "Randomizing" << std::endl;
//            auto theta_random = Eigen::VectorXcd::Random(theta.size()).normalized();
//            auto energy_random = mpstools::common::measure::energy(state,theta)/chain_length;
//            outside_of_window = energy_random < sim_state.energy_lbound or energy_random > sim_state.energy_ubound;
//            while(outside_of_window){
//                std::cout << "Randomizing" << std::endl;
//                theta_random = Eigen::VectorXcd::Random(theta.size()).normalized();
//                energy_random = mpstools::common::measure::energy(state,theta)/chain_length;
//                outside_of_window = energy_random < sim_state.energy_lbound or energy_random > sim_state.energy_ubound;
//            }
//            return std::make_tuple(Textra::Matrix_to_Tensor(theta_random, state.dimensions()), energy_random);
//        }
//    }else{
//    }
}


