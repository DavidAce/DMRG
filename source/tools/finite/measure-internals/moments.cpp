#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>
#include <simulation/nmspc_settings.h>

double tools::finite::measure::multisite::moments::first(const Eigen::Tensor<Scalar,3> &mps,const Eigen::Tensor<Scalar, 4> &mpo, const Eigen::Tensor<Scalar,3> &envL, const Eigen::Tensor<Scalar,3> &envR) {
    // This measures the expectation value <M> of some multisite operator M given multisite mps', mpos and corresponding environments.
    // Note that the environments must contain the correct type of mpos
    // Usually this is the energy E = <H>
//    tools::common::profile::t_ene->tic();
    Eigen::Tensor<Scalar, 0> M1 =
        envL
            .contract(mps, Textra::idx({0}, {1}))
            .contract(mpo, Textra::idx({2, 1}, {2, 0}))
            .contract(mps.conjugate(), Textra::idx({3, 0}, {0, 1}))
            .contract(envR, Textra::idx({0, 2, 1}, {0, 1, 2}));
    if(abs(imag(M1(0))) > 1e-10) {
//        tools::log->critical(fmt::format("Expectation value has an imaginary part: {:.16f} + i {:.16f}", std::real(E(0)), std::imag(E(0))));
        throw std::runtime_error(fmt::format("Expectation value has an imaginary part: {:.16f} + i {:.16f}", std::real(M1(0)), std::imag(M1(0))));
    }
    double moment = std::real(M1(0));
    if(std::isnan(moment) or std::isinf(moment)) throw std::runtime_error(fmt::format("First moment is invalid: {}", moment));
//    tools::common::profile::t_ene->toc();
    return moment;
}

double tools::finite::measure::multisite::moments::second(const Eigen::Tensor<Scalar,3> &mps,const Eigen::Tensor<Scalar, 4> &mpo, const Eigen::Tensor<Scalar,4> &envL, const Eigen::Tensor<Scalar,4> &envR) {
    // This measures the second moment <M²> of some multisite operator M given multisite mps', mpos and corresponding environments.
    // Note that the environments must contain the correct type of mpos
    // Usually this is the second moment of the Hamiltonian = <H²>, which corresponds exactly to the variance if the mpos are energy-reduced, since then
    // the variance is Var H = <(H-E)²> = <H²>

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
    double log2chiL = std::log2(mps.dimension(1));
    double log2chiR = std::log2(mps.dimension(2));
    double log2spin = std::log2(mps.dimension(0));
    Eigen::Tensor<Scalar, 0> M2;
    OMP                      omp(settings::threading::num_threads);
    /* clang-format off */
    if(log2spin >= std::max(log2chiL, log2chiR)) {
        if(log2chiL > log2chiR) {
            //            tools::log->trace("H2 path: log2spin > std::max(log2chiL , log2chiR)  and  log2chiL > log2chiR ");
            Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(Textra::array3{1, 0, 2});
            M2.device(omp.dev) =
                mps_shuffled
                    .contract(envL, Textra::idx({0}, {0}))
                    .contract(mpo,  Textra::idx({0, 3}, {2, 0}))
                    .contract(envR, Textra::idx({0, 3}, {0, 2}))
                    .contract(mpo , Textra::idx({2, 1, 4}, {2, 0, 1}))
                    .contract(mps_shuffled.conjugate(), Textra::idx({2, 0, 1}, {1, 0, 2}));
        }
        else {
            //            tools::log->trace("H2 path: log2spin >= std::max(log2chiL , log2chiR) and  log2chiL <= log2chiR ");
            Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(Textra::array3{2, 0, 1});
            M2.device(omp.dev) =
                mps_shuffled
                    .contract(envR, Textra::idx({0}, {0}))
                    .contract(mpo,  Textra::idx({0, 3}, {2, 1}))
                    .contract(envL, Textra::idx({0, 3}, {0, 2}))
                    .contract(mpo,  Textra::idx({2, 4, 1}, {2, 0, 1}))
                    .contract(mps_shuffled.conjugate(), Textra::idx({2, 1, 0}, {1, 2, 0}));
        }

    } else {
        //        tools::log->trace("H2 path: log2spin < std::max(log2chiL , log2chiR)");
        Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(Textra::array3{1, 0, 2});
        M2.device(omp.dev) =
            mps_shuffled.contract(envL, Textra::idx({0}, {0}))
                .contract(mpo, Textra::idx({0, 3}, {2, 0}))
                .contract(mpo, Textra::idx({4, 2}, {2, 0}))
                .contract(envR, Textra::idx({0, 2, 3}, {0, 2, 3}))
                .contract(mps_shuffled.conjugate(), Textra::idx({1, 0, 2}, {1, 0, 2}));
    }
    /* clang-format on */
    double moment = std::real(M2(0));
    if(std::isnan(moment) or std::isinf(moment)) throw std::runtime_error(fmt::format("Second moment is invalid: {}", moment));
    return moment;
}



double tools::finite::measure::multisite::energy_minus_energy_reduced(const class_state_finite & state, const class_model_finite & model, const class_edges_finite & edges) {
    // This measures the bare energy as given by the MPO's.
    // On each MPO the site energy *could* be reduced.
    // If they are reduced, then
    //      < H > = E - E_reduced ~ 0
    // Else
    //      < H > = E
    if(not math::all_equal(state.active_sites,model.active_sites,edges.active_sites))
        throw std::runtime_error(fmt::format("Could not compute energy: active sites are not equal: state {} | model {} | edges {}",state.active_sites,model.active_sites,edges.active_sites));
    if(state.active_sites.empty())
        throw std::runtime_error("Could not compute energy: active sites are empty");
    const auto &mps = state.get_multisite_mps();
    const auto &mpo = model.get_multisite_mpo();
    const auto &env = edges.get_multisite_ene();
    tools::common::profile::t_ene->tic();
    double e_minus_ered = tools::finite::measure::multisite::first_moment(mps, mpo, env.first.get(), env.second.get());
    tools::common::profile::t_ene->toc();
    return e_minus_ered;
}