//
// Created by david on 2019-07-09.
//

#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <general/class_tic_toc.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/finite/opt.h>

using namespace tools::finite::opt::internal;

ceres_base_functor::ceres_base_functor(const class_tensors_finite &tensors, const class_algorithm_status &status) {
    t_bfgs = std::make_unique<class_tic_toc>(settings::profiling::on, 5, "");
    t_vH2  = std::make_unique<class_tic_toc>(settings::profiling::on, 5, "");
    t_vH2v = std::make_unique<class_tic_toc>(settings::profiling::on, 5, "");
    t_vH   = std::make_unique<class_tic_toc>(settings::profiling::on, 5, "");
    t_vHv  = std::make_unique<class_tic_toc>(settings::profiling::on, 5, "");

    length         = tensors.get_length();
    energy_reduced = tensors.model->get_energy_reduced();

    // All energies in status are per site!
    energy_tgt_per_site  = status.energy_tgt_per_site;
    energy_max_per_site  = status.energy_max_per_site;
    energy_min_per_site  = status.energy_min_per_site;
    energy_llim_per_site = status.energy_llim_per_site;
    energy_ulim_per_site = status.energy_ulim_per_site;
    energy_dens_target   = status.energy_dens_target;
    energy_dens_window   = status.energy_dens_window;
    iteration            = status.iter;
}

double ceres_base_functor::get_variance_per_site() const { return variance_per_site; }
double ceres_base_functor::get_energy_per_site() const { return energy_per_site; }
size_t ceres_base_functor::get_count() const { return counter; }
double ceres_base_functor::get_norm() const { return norm; }
int    ceres_base_functor::NumParameters() const { return num_parameters; }
