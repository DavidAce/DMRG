
#include <tools/finite/measure.h>
#include <iomanip>
#include <config/nmspc_settings.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>
#include <algorithms/class_algorithm_status.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>



double tools::finite::measure::twosite::energy_minus_energy_reduced(const class_state_finite &state, const Eigen::Tensor<Scalar,4> & theta){
    // This measures the bare energy as given by the MPO's.
    // On each MPO the site energy *could* be reduced.
    // If they are reduced, then
    //      < H > = E - E_reduced ~ 0
    // Else
    //      < H > = E

    tools::common::profile::t_ene->tic();
    OMP omp(settings::threading::num_threads);
    Eigen::Tensor<Scalar, 0>  E;
    E.device(omp.dev) =
        state.ENV_L.back().block
            .contract(theta,                               idx({0},{1}))
            .contract(state.MPO_L.back()->MPO(),           idx({1,2},{0,2}))
            .contract(state.MPO_R.front()->MPO(),          idx({3,1},{0,2}))
            .contract(theta.conjugate(),                   idx({0,2,4},{1,0,2}))
            .contract(state.ENV_R.front().block,           idx({0,2,1},{0,1,2}));
    if(std::abs(std::imag(E(0))) > 1e-10 ){
        tools::log->critical(fmt::format("Energy has an imaginary part: {:.16f} + i {:.16f}",std::real(E(0)), std::imag(E(0))));
    }
    assert(std::abs(std::imag(E(0))) < 1e-10 and "Energy has an imaginary part!!!");
    double ene = std::real(E(0));
    if (std::isnan(ene) or std::isinf(ene)) throw std::runtime_error(fmt::format("Energy is invalid: {}", ene));
    tools::common::profile::t_ene->toc();
    return  ene;
}


double tools::finite::measure::twosite::energy(const class_state_finite &state, const Eigen::Tensor<Scalar,4> & theta){
    // This measures the actual energy of the system regardless of the reduced/non-reduced state of the MPO's
    // If they are reduced, then
    //      "Actual energy" = (E - E_reduced) + E_reduced = (~0) + E_reduced = E
    // Else
    //      "Actual energy" = (E - E_reduced) + E_reduced = (E)  + 0 = E

    return twosite::energy_minus_energy_reduced(state,theta) + state.get_energy_reduced();
}


double tools::finite::measure::twosite::energy_per_site(const class_state_finite &state, const Eigen::Tensor<Scalar,4> & theta){
    return twosite::energy(state,theta)/static_cast<double>(state.get_length());
}

double tools::finite::measure::twosite::energy_variance(const class_state_finite &state, const Eigen::Tensor<Scalar,4> & theta){
    // Depending on whether the mpo's are reduced or not we get different formulas.
    // If mpo's are reduced:
    //      Var H = <(H-E_red)^2> - <(H-E_red)>^2 = <H^2> - 2<H>E_red + E_red^2 - (<H> - E_red) ^2
    //                                            = H2    - 2*E*E_red + E_red^2 - E^2 + 2*E*E_red - E_red^2
    //                                            = H2    - E^2
    //      so Var H = <(H-E_red)^2> - energy_minus_energy_reduced^2 = H2 - ~0
    //      where H2 is computed with reduced mpo's. Note that ~0 is not exactly zero
    //      because E_red != E necessarily (though they are supposed to be very close)
    // Else:
    //      Var H = <(H - 0)^2> - <H - 0>^2 = H2 - E^2

//    auto dims_theta = theta.dimensions();
//    auto dims_MPO_L = state.MPO_L.back()->MPO().dimensions();
//    auto dims_MPO_R = state.MPO_R.front()->MPO().dimensions();
//    auto dims_ENV2_L = state.ENV2_L.back().block.dimensions();
//    auto dims_ENV2_R = state.ENV2_R.front().block.dimensions();
//
//    tools::log->warn("dims theta {}", dims_theta);
//    tools::log->warn("dims MPO_L {}", dims_MPO_L);
//    tools::log->warn("dims MPO_R {}", dims_MPO_R);
//    tools::log->warn("dims ENV2_L {}", dims_ENV2_L);
//    tools::log->warn("dims ENV2_R {}", dims_ENV2_R);

    tools::common::profile::t_var->tic();
    OMP omp(static_cast<unsigned int>(settings::threading::num_threads));
    Eigen::Tensor<Scalar, 0> H2;
    H2.device(omp.dev) =
        state.ENV2_L.back().block
            .contract(theta                        , idx({0}  ,{1}))
            .contract(state.MPO_L.back()->MPO()    , idx({1,3},{0,2}))
            .contract(state.MPO_R.front()->MPO()   , idx({4,2},{0,2}))
            .contract(state.MPO_L.back()->MPO()    , idx({1,3},{0,2}))
            .contract(state.MPO_R.front()->MPO()   , idx({4,3},{0,2}))
            .contract(theta.conjugate()            , idx({0,3,5},{1,0,2}))
            .contract(state.ENV2_R.front().block   , idx({0,3,1,2},{0,1,2,3}));
    tools::common::profile::t_var->toc();
    double energy;
    if (state.isReduced()){
        energy = tools::finite::measure::twosite::energy_minus_energy_reduced(state,theta);
    }else{
        energy = tools::finite::measure::twosite::energy(state,theta);
    }
    double E2 = energy*energy;
    double var = std::abs(H2(0) - E2);
    if (std::isnan(var) or std::isinf(var)) throw std::runtime_error(fmt::format("Variance is invalid: {}", var));
    if(var < state.lowest_recorded_variance){
        state.lowest_recorded_variance = var;
    }
    return var;
}


double tools::finite::measure::twosite::energy_variance_per_site(const class_state_finite &state, const Eigen::Tensor<Scalar,4> & theta){
    return twosite::energy_variance(state,theta)/static_cast<double>(state.get_length());
}


