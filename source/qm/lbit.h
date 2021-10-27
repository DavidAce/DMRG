#pragma once

#include "gate.h"
#include "qm.h"
#include <complex>
#include <math/tenx/fwd_decl.h>
#include <vector>

namespace qm::lbit {
    /* clang-format off */
    extern Eigen::Tensor<cplx, 2>              get_time_evolution_operator(cplx delta_t, const Eigen::Tensor<cplx, 2> &hamiltonian);
    extern std::vector<qm::SwapGate>           get_time_evolution_swap_gates(cplx delta_t, const std::vector<qm::SwapGate> &hams_nsite);
    extern std::vector<qm::Gate>               get_time_evolution_gates(cplx delta_t, const std::vector<qm::Gate> &hams_nsite);
    extern std::vector<qm::Gate>               get_unitary_2gate_layer(size_t sites, double fmix);
    extern std::vector<Eigen::Tensor<cplx, 2>> get_time_evolution_operators_2site(size_t sites, cplx delta_t, const std::vector<Eigen::Tensor<cplx, 2>> &twosite_hams);
    extern std::vector<Eigen::Tensor<cplx, 2>> get_time_evolution_operators_3site(size_t sites, cplx delta_t, const std::vector<Eigen::Tensor<cplx, 2>> &hams_3site);
    extern std::vector<Eigen::Tensor<cplx, 4>> get_time_evolution_mpos(cplx delta_t, const std::vector<Eigen::Tensor<cplx, 4>> &mpos);
    extern cplx                                get_lbit_exp_value(const std::vector<std::vector<qm::Gate>> &unitary_layers, const Eigen::Matrix2cd &tau, size_t pos_tau, const Eigen::Matrix2cd &sig, size_t pos_sig);
    extern Eigen::Tensor<cplx, 2>              get_lbit_real_overlap(const std::vector<std::vector<qm::Gate>> &unitary_layers, size_t sites);
    extern Eigen::Tensor<cplx, 2>              get_lbit_overlap_averaged(const std::vector<Eigen::Tensor<cplx, 2>> &lbit_overlap_vec);
    extern Eigen::Tensor<cplx, 2>              get_lbit_overlap_permuted(const Eigen::Tensor<cplx, 2> &lbit_overlap);

    extern std::tuple<double, double, std::vector<double>, size_t>
                                               get_characteristic_length_scale(const Eigen::Tensor<cplx, 2> &lbit_overlap_permuted);
    extern std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::Tensor<double, 3>, Eigen::Tensor<double, 4>>
                                               get_lbit_analysis(const std::vector<size_t> &udepth_vec, const std::vector<double> &fmix_vec, size_t sites, size_t reps);
    /* clang-format on */
}