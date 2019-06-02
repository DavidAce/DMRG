//
// Created by david on 2019-05-31.
//



#include <mps_routines/mps_tools/finite/opt.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <model/class_hamiltonian_base.h>
#include <algorithms/class_simulation_state.h>


MPS_Tools::Finite::Opt::internals::guided_functor::guided_functor(
        const class_superblock & superblock, class_simulation_state &sim_state): base_functor(superblock,sim_state)
        {}


double MPS_Tools::Finite::Opt::internals::guided_functor::windowed_func_abs(double x,double window){
    if (std::abs(x) >= window){
        return std::abs(x);
    }else{
        return std::abs(window);
    }
}
double MPS_Tools::Finite::Opt::internals::guided_functor::windowed_grad_abs(double x,double window){
    if (std::abs(x) >= window){
        return sgn(x);
    }else{
        return 0.0;
    }
}



double MPS_Tools::Finite::Opt::internals::guided_functor::windowed_func_pow(double x,double window){
    if (std::abs(x) >= window){
        return std::pow(x,2);
    }else{
        return std::pow(window,2);;
    }
}
double MPS_Tools::Finite::Opt::internals::guided_functor::windowed_grad_pow(double x,double window){
    if (std::abs(x) >= window){
        return 2.0*(x);
    }else{
        return 0.0;
    }
}




double MPS_Tools::Finite::Opt::internals::guided_functor::operator()(const Eigen::VectorXd &v_and_lambdas, Eigen::VectorXd &grad) {
    t_op->tic();
    double vH2v,vHv,vv,var;
    double log10var, fx;
    double energy_func, energy_grad;
    double norm_func,norm_grad;
    Eigen::VectorXd vH, vH2;
    Eigen::Map<const Eigen::VectorXd> v (v_and_lambdas.data(), v_and_lambdas.size()-1);
    auto lambdas = v_and_lambdas.tail(1);
    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {std::tie(vH2,vH2v)  = get_vH2_vH2v(v,superComponents);}
            #pragma omp section
            {std::tie(vH,vHv)    = get_vH_vHv(v,superComponents);}
            #pragma omp section
            {vv     = v.squaredNorm(); }
        }
        #pragma omp barrier

        #pragma omp single
        {
            energy         = vHv/vv;
            energy_dens    = (energy/length - energy_min ) / (energy_max - energy_min);
            energy_offset  = energy_dens - energy_target_dens;
            energy_func    = windowed_func_pow(energy_offset,energy_window);
            energy_grad    = windowed_grad_pow(energy_offset,energy_window);

            var            = vH2v/vv - energy*energy;
            variance       = var;
            var            = var == 0  ? std::numeric_limits<double>::epsilon() : var;

            norm_offset    = vv - 1.0 ;
            norm_func      = windowed_func_pow(norm_offset,1e-1);
            norm_grad      = windowed_grad_pow(norm_offset,1e-1);

            log10var       = std::log10(var);
        }
        #pragma omp barrier

        fx = log10var
             + energy_func * lambdas(0)
             + norm_func;
        auto vv_1  = std::pow(vv,-1);
        auto vv_2  = std::pow(vv,-2);
        auto var_1 = 1.0/var/std::log(10);

        grad.head(v.size())  = var_1 * (2.0*(vH2*vv_1 - v * vH2v * vv_2) - 4.0 * energy * (vH * vv_1 - v * vHv * vv_2))
                               + lambdas(0) * energy_grad * 2.0 * (vH * vv_1 - v * vHv * vv_2)
                               + norm_grad * 2.0 * v;
        grad(grad.size()-1)  = energy_func;
    }


//    std::cout   << std::setprecision(12) << std::fixed
//                << " Variance: "   << std::setw(18)   << log10var
//                << " Energy : "    << std::setw(18)   << energy
//                << " Energy t : "  << std::setw(18)   << energy_target
//                << " Energy w : "  << std::setw(18)   << energy_window
//                << " Energy d : "  << std::setw(18)   << energy_dens
//                << " Energy td : " << std::setw(18)   << energy_target_dens
//                << " Energy o : "  << std::setw(18)   << energy_offset
//                << " norm o : "    << std::setw(18)   << norm_offset
//                << " lambda 0: "   << std::setw(18)   << lambdas(0)
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
        MPS_Tools::log->warn("lambda 0        = {}" , lambdas(0));
        MPS_Tools::log->warn("lambda 1        = {}" , lambdas(1));
        std::cout << "v: \n" << v << std::endl;
        std::cout << "grad: \n" << grad << std::endl;
        throw std::runtime_error("LBFGS: log10 variance is invalid");
    }



    counter++;
    t_op->toc();
    return fx;
}






