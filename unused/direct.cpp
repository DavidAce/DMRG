//
// Created by david on 2019-05-25.
//


#include <spdlog/spdlog.h>
#include <LBFGS.h>
#include <state/tools/finite/opt.h>
#include <general/class_tic_toc.h>
#include <state/class_finite_state.h>
#include <state/class_environment.h>
#include <simulation/class_simulation_status.h>
#include <general/nmspc_random_numbers.h>
#include <simulation/nmspc_settings.h>

#include <variant>

template<typename Scalar, auto rank>
void warn_if_has_imaginary_part(const Eigen::Tensor<Scalar,rank> &tensor, double threshold = 1e-14) {
    Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>> vector (tensor.data(),tensor.size());
    if constexpr (std::is_same<Scalar, std::complex<double>>::value){
        auto imagSum = vector.imag().cwiseAbs().sum();
        if (imagSum > threshold){
            tools::log->warn("Has imaginary part. Sum: " + std::to_string(imagSum));
        }
    }
}


Eigen::Tensor<std::complex<double>,3>
tools::finite::opt::internals::old_direct_optimization(const class_finite_state &state,
                                                       const class_simulation_status &sim_status, OptType optType){
    tools::log->trace("Optimizing in DIRECT mode");
    using Scalar = std::complex<double>;
    t_opt->tic();
    auto theta = state.get_multitheta();
    double energy_0   = tools::finite::measure::multisite::energy_per_site(state,theta);
    double variance_0 = tools::finite::measure::multisite::energy_variance_per_site(state,theta);
    t_opt->toc();
    std::vector<reports::direct_opt_tuple> opt_log;
    opt_log.emplace_back("Initial",theta.size(), energy_0, std::log10(variance_0), 1.0, 0 ,0,t_opt->get_last_time_interval());


    double chain_length    = state.get_length();
    double energy_new,variance_new,overlap_new;
    Eigen::VectorXcd theta_start  = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size());

    using namespace LBFGSpp;
    using namespace tools::finite::opt::internals;
    double fx;
    int niter,counter;
    LBFGSpp::LBFGSSolver<double> solver(params);

    t_opt->tic();
    switch (optType){
        case OptType::CPLX:{
            tools::finite::opt::internals::direct_functor <Scalar>  functor (state, sim_status);
            Eigen::VectorXd  theta_start_cast = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double*> (theta_start.data()), 2*theta_start.size());
            tools::log->trace("Running LBFGS");
            niter = solver.minimize(functor, theta_start_cast, fx);
            theta_start  = Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<Scalar*> (theta_start_cast.data()), theta_start_cast.size()/2).normalized();
            counter      = functor.get_count();
            energy_new   = functor.get_energy() / chain_length;
            variance_new = functor.get_variance()/chain_length;
            break;
        }
        case OptType::REAL:{
            tools::finite::opt::internals::direct_functor <double> functor (state, sim_status);
            Eigen::VectorXd  theta_start_cast = theta_start.real();
            tools::log->trace("Running LBFGS");
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
    tools::log->trace("Finished LBFGS");

    reports::print_report(opt_log);
    reports::print_report(std::make_tuple(
            tools::finite::opt::internals::t_vH2v->get_measured_time(),
            tools::finite::opt::internals::t_vHv->get_measured_time(),
            tools::finite::opt::internals::t_vH2->get_measured_time(),
            tools::finite::opt::internals::t_vH->get_measured_time(),
            tools::finite::opt::internals::t_op->get_measured_time()
            ));

    state.unset_measurements();
    if (variance_new < variance_0){
        tools::log->debug("Returning new theta");
        return  Textra::Matrix_to_Tensor(theta_start, state.active_dimensions());

    }else{
        tools::log->debug("Returning old theta");
        return  theta;
    }

}