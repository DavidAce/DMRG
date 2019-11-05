//
// Created by david on 2019-07-09.
//

#include <general/class_tic_toc.h>
#include <simulation/class_simulation_status.h>
#include <tools/finite/opt.h>
#include <state/class_state_finite.h>
#include <ceres/ceres.h>

using namespace tools::finite::opt::internal;

ceres_base_functor::ceres_base_functor(
        const class_state_finite & state,
        const class_simulation_status & sim_status)
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



double ceres_base_functor::get_variance   ()const{return variance;}
double ceres_base_functor::get_energy     ()const{return energy  ;}
size_t ceres_base_functor::get_count      ()const{return counter;}
double ceres_base_functor::get_norm       ()const{return norm;}
int    ceres_base_functor::NumParameters  ()const{return num_parameters;}






