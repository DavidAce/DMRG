#pragma once

#include <general/nmspc_tensor_omp.h>
class class_state_infinite;
namespace tools::infinite::measure{
    using Scalar = std::complex<double>;
    extern int    length                          (const class_state_infinite & state);
    extern int    bond_dimension                  (const class_state_infinite & state);
    extern double truncation_error                (const class_state_infinite & state);
    extern double norm                            (const class_state_infinite & state);
    extern double energy_mpo                      (const class_state_infinite & state);
    extern double energy_mpo                      (const class_state_infinite & state, const Eigen::Tensor<Scalar,4> &theta);
    extern double energy_per_site_mpo             (const class_state_infinite & state);
    extern double energy_per_site_ham             (const class_state_infinite & state);
    extern double energy_per_site_mom             (const class_state_infinite & state);
    extern double energy_variance_mpo             (const class_state_infinite & state, const Eigen::Tensor<Scalar,4> &theta, double &energy_mpo);
    extern double energy_variance_mpo             (const class_state_infinite & state, const Eigen::Tensor<Scalar,4> &theta);
    extern double energy_variance_mpo             (const class_state_infinite & state);
    extern double energy_variance_per_site_mpo    (const class_state_infinite & state);
    extern double energy_variance_per_site_ham    (const class_state_infinite & state);
    extern double energy_variance_per_site_mom    (const class_state_infinite & state);
    extern double entanglement_entropy            (const class_state_infinite & state);
}