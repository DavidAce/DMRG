#include <math/eig.h>
#include <physics/nmspc_quantum_mechanics.h>
#include <tensors/class_tensors_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/common/views.h>
#include <tools/infinite/measure.h>
#include <general/nmspc_tensor_extra.h>

class_state_infinite::Scalar moment_generating_function(const class_state_infinite & state_original,
                                                        std::vector<Eigen::Tensor<class_state_infinite::Scalar, 2>> &Op_vec
                                                        ) {
    using Scalar                       = class_state_infinite::Scalar;
    class_state_infinite state_evolved = state_original;

    long chi_lim = 5 * state_evolved.chiC();
    for(auto &Op : Op_vec) {
        // Evolve
        Eigen::Tensor<Scalar, 3> mps_evo = Op.contract(state_evolved.get_2site_mps(), Textra::idx({0}, {0}));
        state_evolved.set_mps(mps_evo, chi_lim);
        if(&Op != &Op_vec.back()) { state_evolved.swap_AB(); }
    }

    long sizeLB = state_evolved.chiB() * state_evolved.chiB();
    // Normalize
    Eigen::Tensor<Scalar, 2> transfer_matrix_theta_evn =
        tools::common::views::get_transfer_matrix_theta_evn(state_evolved).reshape(Textra::array2{sizeLB, sizeLB});
    using namespace settings::precision;
    eig::solver solver;
    auto        nev = static_cast<eig::size_type>(1);
    auto        ncv = static_cast<eig::size_type>(eig_max_ncv);

    solver.eigs(transfer_matrix_theta_evn.data(), sizeLB, nev, ncv, eig::Ritz::LM, eig::Form::NSYM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::OFF,
                eig::Dephase::OFF);
    auto eigval                   = eig::view::get_eigval<eig::cplx>(solver.result, 0);
    auto new_theta_evn_normalized = tools::common::views::get_theta_evn(state_evolved, sqrt(eigval));
    auto old_theta_evn_normalized = tools::common::views::get_theta_evn(state_original);
    long sizeL                    = new_theta_evn_normalized.dimension(1) * state_original.chiA();
    long sizeR                    = new_theta_evn_normalized.dimension(3) * state_original.chiB();

    Eigen::Tensor<Scalar, 2> transfer_matrix_G = new_theta_evn_normalized.contract(old_theta_evn_normalized.conjugate(), Textra::idx({0, 2}, {0, 2}))
                                                     .shuffle(Textra::array4{0, 2, 1, 3})
                                                     .reshape(Textra::array2{sizeL, sizeR});
    // Compute the characteristic function G(a).
    solver.eigs(transfer_matrix_G.data(), transfer_matrix_G.dimension(0), nev, ncv, eig::Ritz::LM, eig::Form::NSYM, eig::Side::R, std::nullopt, eig::Shinv::OFF,
                eig::Vecs::OFF, eig::Dephase::OFF);
    //    solver.eig(transfer_matrix_G.data(),(int)transfer_matrix_G.dimension(0), 1, eig_max_ncv, Ritz::LM, Side::R, false);
    auto lambdaG = eig::view::get_eigval<eig::cplx>(solver.result, 0);
    //    t_temp1->toc();
    return lambdaG;
}

double tools::infinite::measure::energy_per_site_mom(const class_tensors_infinite &tensors) {
    if(tensors.measurements.energy_per_site_mom) return tensors.measurements.energy_per_site_mom.value();
    const auto &state = *tensors.state;
    const auto &model = *tensors.model;
    if(state.chiC() <= 2) {
        tensors.measurements.energy_per_site_mom          = std::numeric_limits<double>::quiet_NaN();
        tensors.measurements.energy_variance_per_site_mom = std::numeric_limits<double>::quiet_NaN();
        return tensors.measurements.energy_per_site_mom.value();
    }
    tools::log->trace("Measuring energy mom");
    tools::common::profile::get_default_prof()["t_ene_mom"]->tic();

    Scalar a      = Scalar(0.0, 1.0) * 5e-3;
    auto   SX     = qm::gen_manybody_spin(qm::spinOneHalf::sx, 2);
    auto   SY     = qm::gen_manybody_spin(qm::spinOneHalf::sy, 2);
    auto   SZ     = qm::gen_manybody_spin(qm::spinOneHalf::sz, 2);
    auto   h_evn  = model.get_mpo_siteA().single_site_hamiltonian(0, 2, SX, SY, SZ);
    auto   h_odd  = model.get_mpo_siteB().single_site_hamiltonian(1, 2, SX, SY, SZ);
    auto   Op_vec = qm::timeEvolution::compute_G(a, 4, h_evn, h_odd);

    // The following only works if state.MPS has been normalized! I.e, you have to have run MPS->compute_mps_components() prior.
    Scalar lambdaG                                    = moment_generating_function(state, Op_vec);
    Scalar l                                          = 2.0; // Number of sites in unit cell
    Scalar G                                          = pow(lambdaG, 1.0 / l);
    Scalar logG                                       = std::log(lambdaG) * 1.0 / l;
    Scalar logGc                                      = std::log(conj(lambdaG)) * 1.0 / l;
    Scalar O                                          = (logG - logGc) / (2.0 * a);
    Scalar VarO                                       = 2.0 * std::log(abs(G)) / (a * a);
    tensors.measurements.energy_per_site_mom          = std::real(O);
    tensors.measurements.energy_variance_per_site_mom = std::real(VarO);
    tools::common::profile::get_default_prof()["t_ene_mom"]->toc();
    return tensors.measurements.energy_per_site_mom.value();
}

double tools::infinite::measure::energy_variance_per_site_mom(const class_tensors_infinite &tensors) {
    if(tensors.measurements.energy_variance_per_site_mom) return tensors.measurements.energy_variance_per_site_mom.value();
    tensors.measurements.energy_per_site_mom = tools::infinite::measure::energy_per_site_mom(tensors);
    return tensors.measurements.energy_variance_per_site_mom.value();
}
