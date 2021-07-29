#pragma once
#include "gate.h"
#include "qm.h"
#include <complex>
#include <general/eigen_tensor_fwd_decl.h>
#include <vector>

namespace qm::time {
    extern std::vector<Eigen::Tensor<cplx, 2>> Suzuki_Trotter_1st_order(cplx delta_t, const Eigen::Tensor<cplx, 2> &h_evn, const Eigen::Tensor<cplx, 2> &h_odd);
    extern std::vector<Eigen::Tensor<cplx, 2>> Suzuki_Trotter_2nd_order(cplx delta_t, const Eigen::Tensor<cplx, 2> &h_evn, const Eigen::Tensor<cplx, 2> &h_odd);
    extern std::vector<Eigen::Tensor<cplx, 2>> Suzuki_Trotter_4th_order(cplx delta_t, const Eigen::Tensor<cplx, 2> &h_evn, const Eigen::Tensor<cplx, 2> &h_odd);
    extern std::vector<Eigen::Tensor<cplx, 2>> get_twosite_time_evolution_operators(cplx delta_t, size_t susuki_trotter_order,
                                                                                    const Eigen::Tensor<cplx, 2> &h_evn, const Eigen::Tensor<cplx, 2> &h_odd);
    extern std::vector<Eigen::Tensor<cplx, 2>> compute_G(cplx a, size_t susuki_trotter_order, const Eigen::Tensor<cplx, 2> &h_evn,
                                                         const Eigen::Tensor<cplx, 2> &h_odd);
    extern std::pair<std::vector<qm::Gate>, std::vector<qm::Gate>> get_time_evolution_gates(cplx delta_t, const std::vector<qm::Gate> &hams_nsite);
}