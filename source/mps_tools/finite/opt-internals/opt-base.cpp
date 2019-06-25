//
// Created by david on 2019-03-18.
//
#include <general/class_tic_toc.h>
#include <algorithms/class_simulation_state.h>
#include <mps_tools/finite/opt.h>
#include <mps_state/class_finite_chain_state.h>
#include <LBFGS.h>

using namespace mpstools::finite::opt::internals;

base_functor::base_functor(
        const class_finite_chain_state & state,
        const class_simulation_state & sim_state)
{
    reset_timers();
    length                        = state.get_length();

    //All energies in sim_state are per site!
    energy_target            = sim_state.energy_target;
    energy_max               = sim_state.energy_max;
    energy_min               = sim_state.energy_min;
    energy_lower_bound       = sim_state.energy_lbound;
    energy_upper_bound       = sim_state.energy_ubound;
    energy_target_dens       = sim_state.energy_dens_target;
    energy_window            = sim_state.energy_dens_window;
    iteration                = sim_state.iteration;
}



double base_functor::get_variance()const{return variance;}
double base_functor::get_energy  ()const{return energy  ;}
size_t base_functor::get_count   ()const{return counter;}
double base_functor::get_norm    ()const{return norm;}





