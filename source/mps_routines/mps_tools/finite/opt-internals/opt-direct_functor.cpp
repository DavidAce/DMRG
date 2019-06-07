//
// Created by david on 2019-05-31.
//

#include <mps_routines/mps_tools/finite/opt.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <model/class_hamiltonian_base.h>
#include <algorithms/class_simulation_state.h>





MPS_Tools::Finite::Opt::internals::direct_functor::direct_functor(
        const class_superblock & superblock_, const class_simulation_state &sim_state):
        base_functor(superblock_,sim_state){}


double MPS_Tools::Finite::Opt::internals::direct_functor::operator()(const Eigen::VectorXd &v, Eigen::VectorXd &grad) {
    t_op->tic();
    long double vH2v,vHv,vv,var;
    double lambda,log10var, fx;
    Eigen::VectorXd vH, vH2;

    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {std::tie(vH2,vH2v)  = get_vH2_vH2v(v,superComponents);}
            #pragma omp section
            {std::tie(vH,vHv)    = get_vH_vHv(v,superComponents);}
            #pragma omp section
            {vv     = v.squaredNorm();}
        }
        #pragma omp barrier

        #pragma omp single
        {
            lambda         = 1.0;
            var            = std::abs(vH2v*vv - vHv*vHv)/std::pow(vv,2);
            variance       = var;
            var            = var == 0  ? std::numeric_limits<double>::epsilon() : var;
            energy         = vHv/vv;
            norm_offset    = vv-1.0;
            log10var       = std::log10(var);
        }
        #pragma omp barrier
        fx = log10var  + lambda * std::pow(norm_offset,2);
        auto vv_1  = std::pow(vv,-1);
        auto vv_2  = std::pow(vv,-2);
        auto var_1 = 1.0/var/std::log(10);

        grad = var_1 * (2.0*(vH2*vv_1 - v * vH2v * vv_2) - 4.0 * energy * (vH * vv_1 - v * vHv * vv_2))
                               + lambda * 2.0 * norm_offset * 2.0 * v;
    }

//    std::cout   << std::setprecision(12) << std::fixed
//                << " Variance: "   << std::setw(18)   << log10var
//                << " Energy : "    << std::setw(18)   << energy
//                << " Energy t : "  << std::setw(18)   << energy_target
//                << " Energy w : "  << std::setw(18)   << energy_density_window
//                << " Energy d : "  << std::setw(18)   << energy_dens
//                << " Energy td : " << std::setw(18)   << energy_target_dens
//                << " Energy o : "  << std::setw(18)   << energy_offset
//                << " norm o : "    << std::setw(18)   << norm_offset
//                << " lambda : "    << std::setw(18)   << lambda
//                << " fx : "        << std::setw(18)   << fx
//                << std::endl;


    if(std::isnan(log10var) or std::isinf(log10var)){
        MPS_Tools::log->warn("log10 variance is invalid");
        MPS_Tools::log->warn("energy offset   = {}" , energy_offset );
        MPS_Tools::log->warn("norm   offset   = {}" , norm_offset );
        MPS_Tools::log->warn("vH2v            = {}" , vH2v );
        MPS_Tools::log->warn("vHv             = {}" , vHv  );
        MPS_Tools::log->warn("vv              = {}" , vv   );
        MPS_Tools::log->warn("vH2v/vv         = {}" , vH2v/vv    );
        MPS_Tools::log->warn("vEv*vEv/vv/vv   = {}" , vHv*vHv/vv/vv    );
        MPS_Tools::log->warn("var             = {}" , var);
        MPS_Tools::log->warn("lambda          = {}" , lambda);
        std::cout << "v: \n" << v << std::endl;
        std::cout << "grad: \n" << grad << std::endl;
        throw std::runtime_error("LBFGS: log10 variance is invalid");
    }


    counter++;
    t_op->toc();
    return fx;
}



