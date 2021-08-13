#include "ceres_base.h"
#include <algorithms/AlgorithmStatus.h>
#include <config/settings.h>
#include <tensors/model/ModelFinite.h>
#include <tensors/state/StateFinite.h>
#include <tensors/TensorsFinite.h>
#include <tid/tid.h>

using namespace tools::finite::opt::internal;

ceres_base_functor::ceres_base_functor(const TensorsFinite &tensors, const AlgorithmStatus &status) {
    t_step = std::make_unique<tid::ur>();
    t_H2n  = std::make_unique<tid::ur>();
    t_nH2n = std::make_unique<tid::ur>();
    t_Hn   = std::make_unique<tid::ur>();
    t_nHn  = std::make_unique<tid::ur>();

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

double ceres_base_functor::get_energy() const { return energy; }
double ceres_base_functor::get_energy_per_site() const { return energy_per_site; }
double ceres_base_functor::get_variance() const { return variance; }
double ceres_base_functor::get_variance_per_site() const { return variance_per_site; }
size_t ceres_base_functor::get_count() const { return counter; }
double ceres_base_functor::get_norm() const { return norm; }
double ceres_base_functor::get_norm_offset() const { return norm_offset; }
double ceres_base_functor::get_delta_f() const { return delta_f; }
double ceres_base_functor::get_max_grad_norm() const { return max_grad_norm; }
long   ceres_base_functor::get_ops() const { return ops; }
int    ceres_base_functor::NumParameters() const { return num_parameters; }

void ceres_base_functor::set_delta_f(double delta_f_) const { delta_f = delta_f_; }
void ceres_base_functor::set_max_grad_norm(double max_grad_norm_) const { max_grad_norm = max_grad_norm_; }