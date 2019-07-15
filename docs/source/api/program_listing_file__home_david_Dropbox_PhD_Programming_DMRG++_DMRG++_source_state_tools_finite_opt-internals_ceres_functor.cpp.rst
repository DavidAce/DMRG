
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_tools_finite_opt-internals_ceres_functor.cpp:

Program Listing for File ceres_functor.cpp
==========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_tools_finite_opt-internals_ceres_functor.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/state/tools/finite/opt-internals/ceres_functor.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-07-09.
   //
   
   #include <general/class_tic_toc.h>
   #include <simulation/class_simulation_status.h>
   #include <state/tools/finite/opt.h>
   #include <state/class_finite_state.h>
   #include <ceres/ceres.h>
   
   using namespace tools::finite::opt::internals;
   
   ceres_functor::ceres_functor(
           const class_finite_state & state,
           const class_simulation_status & sim_status)
           : multiComponents(state)
   {
       reset_timers();
       length                   = state.get_length();
   
       //All energies in sim_status are per site!
       energy_target            = sim_status.energy_target;
       energy_max               = sim_status.energy_max;
       energy_min               = sim_status.energy_min;
       energy_lower_bound       = sim_status.energy_lbound;
       energy_upper_bound       = sim_status.energy_ubound;
       energy_target_dens       = sim_status.energy_dens_target;
       energy_window            = sim_status.energy_dens_window;
       iteration                = sim_status.iteration;
   }
   
   double ceres_functor::get_variance()const{return variance;}
   double ceres_functor::get_energy  ()const{return energy  ;}
   size_t ceres_functor::get_count   ()const{return counter;}
   double ceres_functor::get_norm    ()const{return norm;}
   int    ceres_functor::NumParameters()const{return multiComponents.num_parameters;}
   
   
   
   
   bool tools::finite::opt::internals::ceres_functor::Evaluate(const double* v_double_double,
                                                               double* fx,
                                                               double* grad_double_double) const {
       t_op->tic();
       Scalar vH2v,vHv;
       Scalar ene,var;
       double vv,log10var;
       double energy_func, energy_grad;
       double norm_func,norm_grad;
       VectorType Hv, H2v;
       int vecSize = NumParameters();
   //    if constexpr (std::is_same<Scalar,std::complex<double>>::value){vecSize = v_double_double.size()/2;}
       auto v    = Eigen::Map<const VectorType> (reinterpret_cast<const Scalar*>(v_double_double)   , vecSize);
       auto grad = Eigen::Map<      VectorType> (reinterpret_cast<      Scalar*>(grad_double_double), vecSize);
       auto lambdas =    Eigen::VectorXd::Ones(1) ;
       vv    = v.squaredNorm();
       norm  = std::sqrt(vv);
   
   
       #pragma omp parallel
       {
           #pragma omp sections
           {
               #pragma omp section
               { std::tie(H2v, vH2v) = get_H2v_vH2v(v, multiComponents); }
               #pragma omp section
               { std::tie(Hv, vHv) = get_Hv_vHv(v, multiComponents); }
           }
       }
   
   
   
       ene             = vHv/vv;
       var             = vH2v/vv - ene*ene;
   //    double loss_of_precision = std::log10(std::abs(ene*ene));
   //    double expected_error    = std::pow(10, -(13-loss_of_precision));
   //    if (std::imag(ene)      > expected_error) tools::log->warn("Energy has imaginary component              : {:.16f} + i {:.16f}" , std::real(ene)    , std::imag(ene));
   //    if (std::imag(vH2v/vv)  > expected_error) tools::log->warn("Hamiltonian squared has imaginary component : {:.16f} + i {:.16f}" , std::real(vH2v/vv), std::imag(vH2v/vv));
   //    if (std::imag(var)      > expected_error) tools::log->warn("Variance has imaginary component            : {:.16f} + i {:.16f}" , std::real(var)    , std::imag(var));
       if (std::real(var)      < 0.0           ) tools::log->warn("Counter = {}. Variance is negative:  {:.16f} + i {:.16f}" , counter, std::real(var)    , std::imag(var));
   
       energy         = std::real(ene);
       energy_dens    = (energy/length - energy_min ) / (energy_max - energy_min);
       energy_offset  = energy_dens - energy_target_dens;
       energy_func    = windowed_func_pow(energy_offset,energy_window);
       energy_grad    = windowed_grad_pow(energy_offset,energy_window);
   
       variance       = std::abs(var);
       variance       = variance < 1e-15  ? 1e-15 : variance;
   
       norm_offset    = std::abs(vv) - 1.0 ;
       norm_func      = windowed_func_pow(norm_offset,0.0);
       norm_grad      = windowed_grad_pow(norm_offset,0.0);
   
       log10var       = std::log10(variance);
   
   
       fx[0] = log10var
            + energy_func * lambdas(0)
            + norm_func;
       auto vv_1  = std::pow(vv,-1);
       auto var_1 = 1.0/variance/std::log(10);
       if (grad_double_double){
           grad = var_1 * vv_1 * (H2v  - v  * vH2v - 2.0 * ene * (Hv - v * ene))
                  + lambdas(0) * energy_grad * vv_1 * (Hv - v * ene)
                  +  norm_grad * v;
       }
   
   
   //        grad(grad.size()-1)  = energy_func;
   //    if constexpr (std::is_same<Scalar,std::complex<double>>::value){grad*=2; vecSize = grad.size()*2;}
   //    grad_double_double  = Eigen::Map<Eigen::VectorXd> (reinterpret_cast<double*> (grad.data()), vecSize);
   
   
   
   //    std::cout   << std::setprecision(12) << std::fixed
   //                << " Variance: "   << std::setw(18)   << log10var
   //                << " Energy : "    << std::setw(18)   << energy
   //                << " Energy t : "  << std::setw(18)   << energy_target
   //                << " Energy w : "  << std::setw(18)   << energy_density_window
   //                << " Energy d : "  << std::setw(18)   << energy_dens
   //                << " Energy td : " << std::setw(18)   << energy_target_dens
   //                << " Energy o : "  << std::setw(18)   << energy_offset
   //                << " norm o : "    << std::setw(18)   << norm_offset
   //                << " lambda 0: "   << std::setw(18)   << lambdas(0)
   //                << " fx : "        << std::setw(18)   << fx
   //                << std::endl;
   
   
       if(std::isnan(log10var) or std::isinf(log10var)){
           tools::log->warn("log10 variance is invalid");
           tools::log->warn("vv              = {:.16f} + i{:.16f}" , std::real(vv)  , std::imag(vv));
           tools::log->warn("vH2v            = {:.16f} + i{:.16f}" , std::real(vH2v) ,std::imag(vH2v) );
           tools::log->warn("vHv             = {:.16f} + i{:.16f}" , std::real(vHv)  ,std::imag(vHv)  );
           tools::log->warn("var             = {:.16f} + i{:.16f}" , std::real(var)  ,std::imag(var));
           tools::log->warn("ene             = {:.16f} + i{:.16f}" , std::real(ene)  ,std::imag(ene));
           tools::log->warn("log10(var/L)    = {:.16f}" , std::log10(variance/length) );
           tools::log->warn("energy offset   = {:.16f}" , energy_offset );
           tools::log->warn("norm   offset   = {:.16f}" , norm_offset );
           tools::log->warn("lambda 0        = {:.16f}" , lambdas(0));
           throw std::runtime_error("Direct functor failed at counter = " + std::to_string(counter) );
       }
   
   
   
       counter++;
       t_op->toc();
       return true;
   }
   
