#include <tensors/edges/class_edges_finite.h>
#include <edges/class_env_ene.h>
#include <edges/class_env_var.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>
#include <iomanip>
#include <math/nmspc_math.h>
#include <tensors/model/class_model_finite.h>
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>



//
//double tools::finite::measure::multisite::energy_variance(const class_state_finite &state, const Eigen::Tensor<Scalar, 3> &multitheta) {
//    // Depending on whether the mpo's are reduced or not we get different formulas.
//    // If mpo's are reduced:
//    //      Var H = <(H-E_red)^2> - <(H-E_red)>^2 = <H^2> - 2<H>E_red + E_red^2 - (<H> - E_red) ^2
//    //                                            = H2    - 2*E*E_red + E_red^2 - E^2 + 2*E*E_red - E_red^2
//    //                                            = H2    - E^2
//    //      so Var H = <(H-E_red)^2> - energy_minus_energy_reduced^2 = H2 - ~0
//    //      where H2 is computed with reduced mpo's. Note that ~0 is not exactly zero
//    //      because E_red != E necessarily (though they are supposed to be very close)
//    // Else:
//    //      Var H = <(H - 0)^2> - <H - 0>^2 = H2 - E^2
//    tools::common::profile::t_var->tic();
//    auto  multimpo = state.get_multimpo();
//    auto &env2L    = state.get_ENV2L(state.active_sites.front()).block;
//    auto &env2R    = state.get_ENV2R(state.active_sites.back()).block;
//
//    auto                     dsizes   = state.active_dimensions();
//    double                   log2chiL = std::log2(dsizes[1]);
//    double                   log2chiR = std::log2(dsizes[2]);
//    double                   log2spin = std::log2(dsizes[0]);
//    Eigen::Tensor<Scalar, 0> H2;
//    OMP                      omp(settings::threading::num_threads);
//    if(log2spin >= std::max(log2chiL, log2chiR)) {
//        if(log2chiL > log2chiR) {
//            //            tools::log->trace("H2 path: log2spin > std::max(log2chiL , log2chiR)  and  log2chiL > log2chiR ");
//            Eigen::Tensor<Scalar, 3> theta = multitheta.shuffle(Textra::array3{1, 0, 2});
//            H2.device(omp.dev)             = theta.contract(env2L, Textra::idx({0}, {0}))
//                                     .contract(multimpo, Textra::idx({0, 3}, {2, 0}))
//                                     .contract(env2R, Textra::idx({0, 3}, {0, 2}))
//                                     .contract(multimpo, Textra::idx({2, 1, 4}, {2, 0, 1}))
//                                     .contract(theta.conjugate(), Textra::idx({2, 0, 1}, {1, 0, 2}));
//        }
//
//        else {
//            //            tools::log->trace("H2 path: log2spin >= std::max(log2chiL , log2chiR) and  log2chiL <= log2chiR ");
//            Eigen::Tensor<Scalar, 3> theta = multitheta.shuffle(Textra::array3{2, 0, 1});
//            H2.device(omp.dev)             = theta.contract(env2R, Textra::idx({0}, {0}))
//                                     .contract(multimpo, Textra::idx({0, 3}, {2, 1}))
//                                     .contract(env2L, Textra::idx({0, 3}, {0, 2}))
//                                     .contract(multimpo, Textra::idx({2, 4, 1}, {2, 0, 1}))
//                                     .contract(theta.conjugate(), Textra::idx({2, 1, 0}, {1, 2, 0}));
//        }
//
//    } else {
//        //        tools::log->trace("H2 path: log2spin < std::max(log2chiL , log2chiR)");
//        Eigen::Tensor<Scalar, 3> theta = multitheta.shuffle(Textra::array3{1, 0, 2});
//        H2.device(omp.dev)             = theta.contract(env2L, Textra::idx({0}, {0}))
//                                 .contract(multimpo, Textra::idx({0, 3}, {2, 0}))
//                                 .contract(multimpo, Textra::idx({4, 2}, {2, 0}))
//                                 .contract(env2R, Textra::idx({0, 2, 3}, {0, 2, 3}))
//                                 .contract(theta.conjugate(), Textra::idx({1, 0, 2}, {1, 0, 2}));
//    }
//
//    //
//    //
//    //
//    //    Eigen::Tensor<Scalar, 0> H2 =
//    //            env2L
//    //            .contract(multitheta                 , idx({0}  ,{1}))
//    //            .contract(multimpo                   , idx({3,1},{2,0}))
//    //            .contract(multimpo                   , idx({4,1},{2,0}))
//    //            .contract(multitheta.conjugate()     , idx({4,0},{0,1}))
//    //            .contract(env2R                      , idx({0,3,1,2},{0,1,2,3}));
//    tools::common::profile::t_var->toc();
//    double energy;
//    if(state.is_reduced()) {
//        energy = multisite::energy_minus_energy_reduced(state, multitheta);
//    } else {
//        energy = multisite::energy(state, multitheta);
//    }
//    double E2  = energy * energy;
//    double var = std::abs(H2(0) - E2);
//    if(std::isnan(var) or std::isinf(var)) throw std::runtime_error(fmt::format("Variance is invalid: {}", var));
//    internal::significant_digits(std::abs(H2(0)), E2);
//    if(var < state.lowest_recorded_variance) {
//        state.lowest_recorded_variance = var;
//    }
//    return var;
//}

//double tools::finite::measure::multisite::energy_variance_per_site(const class_state_finite &state, const Eigen::Tensor<Scalar, 3> &multitheta) {
//    return multisite::energy_variance(state, multitheta) / static_cast<double>(state.get_length());
//}



//
//double tools::finite::measure::multisite::energy(const class_state_finite &state) {
//    if(state.measurements.energy) return state.measurements.energy.value();
//    if(state.active_sites.empty()) return tools::finite::measure::energy(state);
//    tools::common::profile::t_ene->tic();
//    auto theta = state.get_multisite_tensor();
//    tools::common::profile::t_ene->toc();
//    state.measurements.energy = multisite::energy(state, theta);
//    return state.measurements.energy.value();
//}
//
//double tools::finite::measure::multisite::energy_per_site(const class_state_finite &state) {
//    if(state.measurements.energy_per_site) {
//        return state.measurements.energy_per_site.value();
//    } else {
//        if(state.active_sites.empty()) return tools::finite::measure::energy_per_site(state);
//        state.measurements.energy_per_site = multisite::energy(state) / static_cast<double>(state.get_length());
//        return state.measurements.energy_per_site.value();
//    }
//}
//
//double tools::finite::measure::multisite::energy_variance(const class_state_finite &state) {
//    if(state.measurements.energy_variance) {
//        return state.measurements.energy_variance.value();
//    } else {
//        if(state.active_sites.empty()) return tools::finite::measure::energy_variance(state);
//        tools::common::profile::t_var->tic();
//        auto theta = state.get_multisite_tensor();
//        tools::common::profile::t_var->toc();
//        state.measurements.energy_variance = multisite::energy_variance(state, theta);
//        return state.measurements.energy_variance.value();
//    }
//}
//
//double tools::finite::measure::multisite::energy_variance_per_site(const class_state_finite &state) {
//    if(state.measurements.energy_variance_per_site) {
//        return state.measurements.energy_variance_per_site.value();
//    } else {
//        if(state.active_sites.empty()) return tools::finite::measure::energy_variance_per_site(state);
//        state.measurements.energy_variance_per_site = multisite::energy_variance(state) / static_cast<double>(state.get_length());
//        return state.measurements.energy_variance_per_site.value();
//    }
//}
