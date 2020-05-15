//
// Created by david on 2019-07-09.
//

#include <ceres/ceres.h>
#include <general/class_tic_toc.h>
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <tensors/state/class_state_finite.h>
#include <tensors/class_tensors_finite.h>
#include <tools/finite/opt.h>

using namespace tools::finite::opt::internal;

ceres_base_functor::ceres_base_functor(
        const class_tensors_finite &tensors,
        const class_algorithm_status & status):
        omp(settings::threading::num_threads)
{
    reset_timers();
    length                   = tensors.get_length();

    //All energies in status are per site!
    energy_target            = status.energy_target;
    energy_max               = status.energy_max;
    energy_min               = status.energy_min;
    energy_lower_bound       = status.energy_lbound;
    energy_upper_bound       = status.energy_ubound;
    energy_target_dens       = status.energy_dens_target;
    energy_window            = status.energy_dens_window;
    iteration                = status.iter;
}



double ceres_base_functor::get_variance   ()const{return variance;}
double ceres_base_functor::get_energy     ()const{return energy  ;}
size_t ceres_base_functor::get_count      ()const{return counter;}
double ceres_base_functor::get_norm       ()const{return norm;}
int    ceres_base_functor::NumParameters  ()const{return num_parameters;}






