//
// Created by david on 2019-03-18.
//

#include <spdlog/spdlog.h>
#include <LBFGS.h>
#include <iomanip>
#include <mps_routines/mps_tools/finite/opt.h>
#include <general/class_tic_toc.h>
#include <mps_routines/class_superblock.h>
#include <model/class_hamiltonian_base.h>
#include <mps_routines/class_environment.h>





std::tuple<Eigen::Tensor<std::complex<double>,4>, double>
MPS_Tools::Finite::Opt::internals::direct_optimization(const class_superblock & superblock){
    spdlog::trace("Optimizing in DIRECT mode");
    using Scalar = std::complex<double>;
    class_tic_toc t_lbfgs(true,5,"lbfgs");
    t_opt->tic();
    double chain_length    = superblock.get_length();
    double energy_new,variance_new,overlap_new;
    //Should really use xstart as the projection towards the previous theta, not best overlapping!
    // Note that alpha_i = <theta_old | theta_new_i> is not supposed to be squared! The overlap
    // Between xstart and theta_old should be

    auto theta = superblock.get_theta();
    Eigen::VectorXd xstart = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size()).real();

    std::stringstream problem_report;
    problem_report
            << "Starting optimization \n"
            << std::setprecision(10)
            << "      mode        : "    << OptSpace::DIRECT << '\n'
            << "      position    : "    << superblock.get_position() << '\n'
            << "      chi         : "    << superblock.get_chi() << '\n'
            << "      shape       : "    << theta.size() << " x " << theta.size() << '\n' << '\n' << std::flush;
    spdlog::debug(problem_report.str());


    std::vector<std::tuple<std::string,int,double,Scalar,double,int,int,double,double>> opt_log;
    {
        int iter_0 = 0;
        double energy_0   = MPS_Tools::Common::Measure::energy_mpo(superblock,theta);
        double variance_0 = MPS_Tools::Common::Measure::energy_variance_mpo(superblock,theta,energy_0);
        t_opt->toc();
        opt_log.emplace_back("Start (best overlap)",theta.size(), energy_0/chain_length, std::log10(variance_0/chain_length), 1.0, iter_0 ,0,t_opt->get_last_time_interval(),  t_tot->get_age());
        t_opt->tic();
        using namespace LBFGSpp;
        MPS_Tools::Finite::Opt::internals::direct_functor functor (superblock);
//        functor.set_energy_bounds(sim_state.energy_ubound,sim_state.energy_lbound);
//        double threshold = 1e-5;
        LBFGSpp::LBFGSParam<double> param;
        param.max_iterations = 2000;
        param.max_linesearch = 60; // Default is 20. 5 is really bad, 80 seems better.
        param.m              = 8;
        param.epsilon        = 1e-3;  // Default is 1e-5.
        param.delta          = 1e-6;  // Default is 0. Trying this one instead of ftol.
        param.ftol           = 1e-3;  // Default is 1e-4. this really helped at threshold 1e-8. Perhaps it should be low. Ok..it didn't
        param.past           = 1;     // Or perhaps it was this one that helped.
        param.wolfe          = 5e-1;
        param.min_step       = 1e-40;
        param.max_step       = 1e+40;
        // Create solver and function object
        LBFGSpp::LBFGSSolver<double> solver_3(param);
        // x will be overwritten to be the best point found
        double fx;
        spdlog::trace("Running LBFGS");
        int niter = solver_3.minimize(functor, xstart, fx);
        int counter = functor.get_count();
        t_opt->toc();
        xstart.normalize();
        auto theta_old = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size()).real();

        energy_new   = functor.get_energy() / chain_length;
        variance_new = functor.get_variance()/chain_length;
        overlap_new  = (theta_old.adjoint() * xstart).cwiseAbs().sum();
//        opt_log.emplace_back("LBFGS++",theta.size(), energy_new, std::log10(variance_new), overlap_new, niter,counter, t_opt->get_last_time_interval(), 0);
        opt_log.emplace_back("LBFGS++",theta.size(), energy_new, std::log10(variance_new), overlap_new, niter,counter, t_opt->get_last_time_interval(), t_tot->get_age());
        spdlog::trace("Time in function = {:.3f}", t_opt->get_measured_time()*1000);
        spdlog::trace("Finished LBFGS");

    }
    std::stringstream report;
    report    << std::setprecision(16) << '\n'
              <<"    "<< std::setw(24) << std::left << "Algorithm"
              <<"    "<< std::setw(8)  << std::left << "size"
              <<"    "<< std::setw(24) << std::left << "energy"
              <<"    "<< std::setw(44) << std::left << "variance"
              <<"    "<< std::setw(24) << std::left << "overlap"
              <<"    "<< std::setw(8)  << std::left << "iter"
              <<"    "<< std::setw(8)  << std::left << "counter"
              <<"    "<< std::setw(20) << std::left << "Elapsed time [ms]"
              <<"    "<< std::setw(20) << std::left << "Time per count [ms]"
              <<"    "<< std::setw(20) << std::left << "Wall time [s]"
              << '\n';
    for(auto &log : opt_log){
        report   << std::setprecision(16)
                 << "    " << std::setw(24) << std::left << std::fixed << std::get<0>(log)
                 << "    " << std::setw(8)  << std::left << std::fixed << std::get<1>(log)
                 << "    " << std::setw(24) << std::left << std::fixed << std::get<2>(log)
                 << "    " << std::setw(44) << std::left << std::fixed << std::get<3>(log)
                 << "    " << std::setw(24) << std::left << std::fixed << std::get<4>(log)
                 << "    " << std::setw(8)  << std::left << std::fixed << std::get<5>(log) << std::setprecision(3)
                 << "    " << std::setw(8)  << std::left << std::fixed << std::get<6>(log) << std::setprecision(3)
                 << "    " << std::setw(20) << std::left << std::fixed << std::get<7>(log)*1000
                 << "    " << std::setw(20) << std::left << std::fixed << std::get<7>(log)*1000 / (double)std::get<6>(log)
                 << "    " << std::setw(20) << std::left << std::fixed << std::get<8>(log)
                 << '\n';
    }
    report << '\n';
    spdlog::debug(report.str());

    return  std::make_tuple(Textra::Matrix_to_Tensor(xstart.cast<Scalar>(), superblock.dimensions()), energy_new);
}











MPS_Tools::Finite::Opt::internals::direct_functor::direct_functor(
        const class_superblock & superblock_): base_functor()

{
    superblock.HA_MPO  = superblock_.HA->MPO().real();
    superblock.HB_MPO  = superblock_.HB->MPO().real();
    superblock.Lblock  = superblock_.Lblock->block.real();
    superblock.Rblock  = superblock_.Rblock->block.real();
    superblock.Lblock2 = superblock_.Lblock2->block.real();
    superblock.Rblock2 = superblock_.Rblock2->block.real();
    superblock.dsizes  = superblock_.dimensions();
    superblock.HAHB    = superblock.HA_MPO.contract(superblock.HB_MPO, Textra::idx({1},{0}));
    superblock.HAHB2   = superblock.HAHB.contract(superblock.HAHB, Textra::idx({2,5},{1,4}));
}





double MPS_Tools::Finite::Opt::internals::direct_functor::operator()(const Eigen::VectorXd &v, Eigen::VectorXd &grad) {
    long double vH2v,vHv,vv,var;
    double lambda,lambda2,log10var, fx,energy_penalty;
    Eigen::VectorXd vH, vH2;

    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {vH2v   = get_vH2v(v,superblock);}
            #pragma omp section
            {vHv    = get_vHv(v,superblock);}
            #pragma omp section
            {vv     = v.cwiseAbs2().sum();}
            #pragma omp section
            {vH2    = get_vH2(v,superblock);}
            #pragma omp section
            {vH     = get_vH(v,superblock);}
        }
        #pragma omp barrier
        #pragma omp single
        {
            lambda         = 1.0;
            lambda2        = 1.0;
            var            = std::abs(vH2v*vv - vHv*vHv)/std::pow(vv,2);
            variance       = var;
            var            = var == 0  ? std::numeric_limits<double>::epsilon() : var;
            energy         = vHv/vv;
            double energy_distance = std::abs(energy - energy_target);
            energy_penalty = energy_distance > energy_window/2.0 ?  energy - energy_target : 0.0;
            log10var    = std::log10(var);

        }
        #pragma omp barrier
        #pragma omp for schedule(static,1)
        for (int k = 0; k < v.size(); k++){
            double vi2H2ik         = 2.0*vH2(k);             // 2 * sum_i x_i H2_ik
            double vk4EkvEv        = 4.0*vH(k) * vHv ;       // 4 * x_k * E_k * (sum_i x_i^2 E_i)
            grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vHv*vHv*4.0*v(k)*vv)/(std::pow(vv,4)))/var/std::log(10)
                      + lambda * 4.0 * v(k) * (vv - 1);
            if (have_bounds_on_energy){
//                grad(k) +=  (2*vH(k)/vv - 2*v(k)*vHv/vv/vv);
//                grad(k) +=  sgn(energy-energy_target) *2/vv*(vH(k) - v(k)*energy);
                grad(k) += lambda2* 4 * energy_penalty*(vH(k)/vv - energy*v(k)/vv);

            }
            grad(k) = std::isinf(grad(k)) ? std::numeric_limits<double>::max() * sgn(grad(k)) : grad(k);

//            grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vHv*vHv*4.0*v(k)*vv)/(std::pow(vv,4)))/var/std::log(10)
//                      + lambda * 2.0 * v(k) * sgn(vv - 1);//
//            grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vHv*vHv*4.0*v(k)*vv)/(std::pow(vv,4)))/var/std::log(10)
//                      + lambda * 2.0 * v(k) / vv / std::log(10);
        }
    }
//    grad.normalize(); //Seems to stabilize the step size
    if(std::isnan(log10var) or std::isinf(log10var)){
        spdlog::warn("log10 variance is invalid");
        std::cout << "v: \n" << v << std::endl;
        std::cout << "grad: \n" << grad << std::endl;
        spdlog::warn("vH2v            = {}" , vH2v );
        spdlog::warn("vHv             = {}" , vHv  );
        spdlog::warn("vv              = {}" , vv   );
        spdlog::warn("vH2v/vv         = {}" , vH2v/vv    );
        spdlog::warn("vEv*vEv/vv/vv   = {}" , vHv*vHv/vv/vv    );
        spdlog::warn("var             = {}" , var);
        exit(1);
        log10var    = std::abs(var) == 0  ?  -20.0 : std::log10(std::abs(var));
    }
    fx = log10var  + lambda * std::pow(vv-1.0,2);
    if (have_bounds_on_energy) {
        fx += lambda2 * std::pow(energy_penalty,2);
    }


    counter++;
    return fx;
}





double MPS_Tools::Finite::Opt::internals::get_vH2v(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components & superblock){
    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), superblock.dsizes);
    Eigen::Tensor<double, 0> H2 =
            superblock.Lblock2
                    .contract(theta,                Textra::idx({0}  ,{1}))
                    .contract(superblock.HAHB2,     Textra::idx({2,1,3,4},{4,0,1,3}))
                    .contract(theta.conjugate(),    Textra::idx({0,3,5},{1,0,2}))
                    .contract(superblock.Rblock2,   Textra::idx({0,3,1,2},{0,1,2,3})).real();
    return H2(0);
}

double MPS_Tools::Finite::Opt::internals::get_vHv(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components & superblock){
    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), superblock.dsizes);
    Eigen::Tensor<double, 0>  E =
            superblock.Lblock
                    .contract(theta,                           Textra::idx({0},{1}))
                    .contract(superblock.HAHB,                 Textra::idx({1,2,3},{0,1,4}))
                    .contract(theta.conjugate(),               Textra::idx({0,2,4},{1,0,2}))
                    .contract(superblock.Rblock,               Textra::idx({0,2,1},{0,1,2}));
    return E(0);
}


Eigen::VectorXd MPS_Tools::Finite::Opt::internals::get_vH2 (const Eigen::Matrix<double,Eigen::Dynamic,1> &v,const superblock_components & superblock){
    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), superblock.dsizes);
    Eigen::Tensor<double, 4> vH2 =
            superblock.Lblock2
                    .contract(theta,                            Textra::idx({0}  ,{1}))
                    .contract(superblock.HAHB2,                 Textra::idx({2,1,3,4},{4,0,1,3}))
                    .contract(superblock.Rblock2,               Textra::idx({1,2,4},{0,2,3}))
                    .shuffle(Textra::array4{1,0,2,3});
    return Eigen::Map<Eigen::VectorXd>(vH2.data(),vH2.size());
}

Eigen::VectorXd MPS_Tools::Finite::Opt::internals::get_vH (const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components & superblock){
    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), superblock.dsizes);
    Eigen::Tensor<double, 4> vH =
            superblock.Lblock
                    .contract(theta,                               Textra::idx({0},{1}))
                    .contract(superblock.HAHB,                     Textra::idx({1,2,3},{0,1,4}))
                    .contract(superblock.Rblock ,                  Textra::idx({1,3},{0,2}))
                    .shuffle(Textra::array4{1,0,2,3});
    return Eigen::Map<Eigen::VectorXd>(vH.data(),vH.size());
}

