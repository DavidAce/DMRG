//
// Created by david on 2019-03-18.
//
#include <mps_routines/mps_tools/finite/opt.h>


double MPS_Tools::Finite::Opt::internals::base_functor::get_variance()const{return variance;}
double MPS_Tools::Finite::Opt::internals::base_functor::get_energy  ()const{return energy  ;}
size_t MPS_Tools::Finite::Opt::internals::base_functor::get_count   ()const{return counter;}
void   MPS_Tools::Finite::Opt::internals::base_functor::set_energy_bounds(double E_lower, double E_upper){
    energy_lower_bound = E_lower;
    energy_upper_bound = E_upper;
    energy_target = 0.5*(energy_upper_bound + energy_lower_bound);
    energy_window = 1.0*(energy_upper_bound - energy_lower_bound);
    have_bounds_on_energy = true;
}
