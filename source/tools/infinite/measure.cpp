//
// Created by david on 2019-02-01.
//

#include <general/nmspc_quantum_mechanics.h>
#include <math/class_eigsolver.h>
#include <math/class_svd_wrapper.h>
#include <tensors/class_tensors_infinite.h>
#include <tensors/edges/class_edges_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/log.h>
#include <tools/common/moments.h>
#include <tools/common/prof.h>
#include <tools/common/views.h>
#include <tools/infinite/measure.h>

void tools::infinite::measure::do_all_measurements(const class_tensors_infinite &tensors) {
    tensors.measurements.length                       = tools::infinite::measure::length(tensors);
    tensors.measurements.energy_mpo                   = tools::infinite::measure::energy_mpo(tensors);
    tensors.measurements.energy_per_site_mpo          = tools::infinite::measure::energy_per_site_mpo(tensors);
    tensors.measurements.energy_variance_mpo          = tools::infinite::measure::energy_variance_mpo(tensors);
    tensors.measurements.energy_per_site_ham          = tools::infinite::measure::energy_per_site_ham(tensors);
    tensors.measurements.energy_per_site_mom          = tools::infinite::measure::energy_per_site_mom(tensors);
    tensors.measurements.energy_variance_per_site_mpo = tools::infinite::measure::energy_variance_per_site_mpo(tensors);
    tensors.measurements.energy_variance_per_site_ham = tools::infinite::measure::energy_variance_per_site_ham(tensors);
    tensors.measurements.energy_variance_per_site_mom = tools::infinite::measure::energy_variance_per_site_mom(tensors);
}

void tools::infinite::measure::do_all_measurements(const class_state_infinite &state) {
    state.measurements.norm                 = tools::infinite::measure::norm(state);
    state.measurements.bond_dimension       = tools::infinite::measure::bond_dimension(state);
    state.measurements.entanglement_entropy = tools::infinite::measure::entanglement_entropy(state);
    state.measurements.truncation_error     = tools::infinite::measure::truncation_error(state);
}

//
//
//
//
//
// Scalar moment_generating_function(const class_mps_2site &MPS_original,
//                                                     std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec){
////    t_temp1->tic();
//    std::unique_ptr<class_mps_2site> MPS_evolved = std::make_unique<class_mps_2site>(MPS_original);
//
//    class_SVD<Scalar> SVD;
//    SVD.setThreshold(settings::precision::svd_threshold);

//    long chi_lim = 5*MPS_evolved->chiC();
////    t_temp2->tic();
//    for (auto &Op: Op_vec) {
//        //Evolve
//        Eigen::Tensor<Scalar, 4> theta_evo = Op.contract(MPS_evolved->get_theta(),Textra::idx({0, 1}, {0, 2})).shuffle(array4{0, 2, 1, 3});
//        auto[U, S, V] = SVD.schmidt(theta_evo,chi_lim);
//        MPS_evolved->LC = S;
//        Eigen::Tensor<Scalar,3> L_U =  asDiagonalInversed(MPS_evolved->MPS_A->get_L()).contract(U,Textra::idx({1}, {1})).shuffle(array3{1, 0, 2});
//        Eigen::Tensor<Scalar,3> V_L =  V.contract(asDiagonalInversed(MPS_evolved->MPS_B->get_L()),Textra::idx({2}, {0}));
//        MPS_evolved->MPS_A->set_G(L_U);
//        MPS_evolved->MPS_B->set_G(V_L);
//
//        if (&Op != &Op_vec.back()) {
//            MPS_evolved->enlarge();
//        }
//    }
////    t_temp2->toc();
//
//    long sizeLB = MPS_evolved->chiB() * MPS_evolved->chiB();
//    //Normalize
////    t_temp3->tic();
//    Eigen::Tensor<Scalar,2> transfer_matrix_theta_evn       =
//    tools::common::views::get_transfer_matrix_theta_evn(*MPS_evolved).reshape(array2{sizeLB,sizeLB});
////    t_temp3->toc();
//    using namespace settings::precision;
//
////    t_temp4->tic();
//    class_eigsolver_arpack<Scalar, Form::GENERAL> solver;
//    solver.eig(transfer_matrix_theta_evn.data(),(int)sizeLB, 1, eig_max_ncv, eigsolver_properties::Ritz::LM, eigsolver_properties::Side::R, false);
//    auto new_theta_evn_normalized        = tools::common::views::get_theta_evn(*MPS_evolved, sqrt(solver.ref_eigvals()[0]));
////    t_temp4->toc();
//    long sizeL = new_theta_evn_normalized.dimension(1) * MPS_original.chiA();// theta_evn_normalized.dimension(1);
//    long sizeR = new_theta_evn_normalized.dimension(3) * MPS_original.chiB();// theta_evn_normalized.dimension(3);
//
//    Eigen::Tensor<Scalar,2> transfer_matrix_G =   new_theta_evn_normalized
//            .contract(tools::common::views::theta_evn_normalized.conjugate(),Textra::idx({0,2},{0,2}))
//            .shuffle(array4{0,2,1,3})
//            .reshape(array2{sizeL,sizeR});
//    //Compute the characteristic function G(a).
//    solver.eig(transfer_matrix_G.data(),(int)transfer_matrix_G.dimension(0), 1, eig_max_ncv, Ritz::LM, Side::R, false);
//    Scalar lambdaG = solver.ref_eigvals()[0];
////    t_temp1->toc();
//    return lambdaG;
//}
//
//
//
//
//
//
//
//
//
//
//
// double tools::infinite::measure::energy(const class_state_infinite & state){
////    state.t_ene_mpo->tic();
//    double energy = tools::finite::measure::energy(state);
//    double L      = tools::finite::measure::length(state);
////    state.t_ene_mpo->toc();
//    return energy / L;
//
//}
//
//
//
// double tools::infinite::measure::energy_per_site_ham(const class_state_infinite & state){
//    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
//    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
//    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
//    auto h_evn = state.HA->single_site_hamiltonian(0,2,SX,SY, SZ);
//    auto h_odd = state.HB->single_site_hamiltonian(1,2,SX,SY, SZ);
//    tools::common::views::compute_mps_components(state);
//    using namespace tools::common::views;
//
//    Eigen::Tensor<Scalar,0>
//    E_evn = theta_evn_normalized
//            .contract(Matrix_to_Tensor(h_evn,2,2,2,2), Textra::idx({0, 2}, {0, 1}))
//            .contract(theta_evn_normalized.conjugate(),Textra::idx({2, 3}, {0, 2}))
//            .contract(l_evn,                           Textra::idx({0, 2}, {0, 1}))
//            .contract(r_evn,                           Textra::idx({0, 1}, {0, 1}));
//
//    Eigen::Tensor<Scalar,0>
//    E_odd  = theta_odd_normalized
//            .contract(Matrix_to_Tensor(h_odd,2,2,2,2) ,idx({0, 2}, {0, 1}))
//            .contract(theta_odd_normalized.conjugate(),idx({2, 3}, {0, 2}))
//            .contract(l_odd,                          Textra::idx({0, 2}, {0, 1}))
//            .contract(r_odd,                          Textra::idx({0, 1}, {0, 1}));
//    assert(abs(imag(E_evn(0)+ E_odd(0))) < 1e-10 and "Energy has an imaginary part!!!" );
////    t_ene_ham->toc();
//    return 0.5*std::real(E_evn(0) + E_odd(0));
//
//}
//
//
// double tools::infinite::measure::energy_per_site_mom(const class_state_infinite & state){
////    t_var_gen->tic();
//    Scalar a  = (0.0 + 1.0i) *5e-3;
//    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
//    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
//    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
//    auto h_evn = state.HA->single_site_hamiltonian(0,2,SX,SY, SZ);
//    auto h_odd = state.HB->single_site_hamiltonian(1,2,SX,SY, SZ);
//    auto Op_vec = qm::timeEvolution::compute_G(a,4, h_evn, h_odd);
//
//
//    //The following only works if state.MPS has been normalized! I.e, you have to have run MPS->compute_mps_components() prior.
//    Scalar lambdaG  = moment_generating_function(*state.MPS, Op_vec);
//    Scalar l        = 2.0; //Number of sites in unit cell
//    Scalar G        = pow(lambdaG,1.0/l);
//    Scalar logG     = log(lambdaG) * 1.0/l;
//    Scalar logGc    = log(conj(lambdaG) ) * 1.0/l;
//    Scalar O        = (logG - logGc)/(2.0*a);
////    Scalar VarO     = 2.0*log(abs(G))/ (a*a);
//    return std::real(O);
////    variance_mom    = real(VarO);
////    t_var_gen->toc();
//
//}
//
//
// double tools::infinite::measure::energy_variance_per_site(const class_state_infinite &state) {
//
//    double VarE  = tools::finite::measure::energy_variance(state);
//    double L     = tools::finite::measure::length(state);
//    return VarE/L;
//}
//
// double tools::infinite::measure::energy_variance_per_site_ham(const class_state_infinite &state) {
////    t_var_ham->tic();
//    using namespace tools::common::views;
//
//    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
//    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
//    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
//    auto h_evn = state.HA->single_site_hamiltonian(0,2,SX,SY, SZ);
//    auto h_odd = state.HB->single_site_hamiltonian(1,2,SX,SY, SZ);
//    tools::common::views::compute_mps_components(state);
//
//    Eigen::Tensor<Scalar,0>
//            E_evn = theta_evn_normalized
//            .contract(Matrix_to_Tensor(h_evn,2,2,2,2), Textra::idx({0, 2}, {0, 1}))
//            .contract(theta_evn_normalized.conjugate(),Textra::idx({2, 3}, {0, 2}))
//            .contract(l_evn,                           Textra::idx({0, 2}, {0, 1}))
//            .contract(r_evn,                           Textra::idx({0, 1}, {0, 1}));
//
//    Eigen::Tensor<Scalar,0>
//            E_odd  = theta_odd_normalized
//            .contract(Matrix_to_Tensor(h_odd,2,2,2,2) ,idx({0, 2}, {0, 1}))
//            .contract(theta_odd_normalized.conjugate(),idx({2, 3}, {0, 2}))
//            .contract(l_odd,                          Textra::idx({0, 2}, {0, 1}))
//            .contract(r_odd,                          Textra::idx({0, 1}, {0, 1}));
//
//    Eigen::Tensor<Scalar,4> h0 =  Matrix_to_Tensor((h_evn - E_evn(0)*MatrixType<Scalar>::Identity(4,4)).eval(), 2,2,2,2);
//    Eigen::Tensor<Scalar,4> h1 =  Matrix_to_Tensor((h_odd - E_odd(0)*MatrixType<Scalar>::Identity(4,4)).eval(), 2,2,2,2);
//
//    Eigen::Tensor<Scalar,0> E2AB =
//            theta_evn_normalized
//                    .contract(h0                                , Textra::idx({0, 2}, {0, 1}))
//                    .contract(h0                                , Textra::idx({2, 3}, {0, 1}))
//                    .contract(theta_evn_normalized.conjugate()  , Textra::idx({2, 3}, {0, 2}))
//                    .contract(l_evn                             , Textra::idx({0, 2}, {0, 1}))
//                    .contract(r_evn                             , Textra::idx({0, 1}, {0, 1}));
//
//
//    Eigen::Tensor<Scalar, 0> E2BA =
//            theta_odd_normalized
//                    .contract(h1                              ,Textra::idx({0, 2}, {0, 1}))
//                    .contract(h1                              ,Textra::idx({2, 3}, {0, 1}))
//                    .contract(theta_odd_normalized.conjugate(),Textra::idx({2, 3}, {0, 2}))
//                    .contract(l_odd                           ,Textra::idx({0, 2}, {0, 1}))
//                    .contract(r_odd                           ,Textra::idx({0, 1}, {0, 1}));
//
//
//
//    Eigen::Tensor<Scalar,5> thetaABA = theta_evn_normalized.contract(LAGA,Textra::idx({3},{1}));
//    Eigen::Tensor<Scalar,5> thetaBAB = theta_odd_normalized.contract(LCGB,Textra::idx({3},{1}));
//
//    Eigen::Tensor<Scalar,0> E2ABA_1  =
//            thetaABA
//                    .contract(h1,                  Textra::idx({2,3},{0,1}))
//                    .contract(h0,                  Textra::idx({0,3},{0,1}))
//                    .contract(thetaABA.conjugate(),Textra::idx({3,4,2},{0,2,3}))
//                    .contract(l_evn,               Textra::idx({0,2},{0,1}))
//                    .contract(r_odd,               Textra::idx({0,1},{0,1})) ;
//
//    Eigen::Tensor<Scalar,0> E2BAB_1  =
//            thetaBAB
//                    .contract(h1,                  Textra::idx({0,2},{0,1}))
//                    .contract(h0,                  Textra::idx({4,1},{0,1}))
//                    .contract(thetaBAB.conjugate(),Textra::idx({2,3,4},{0,2,3}))
//                    .contract(l_odd,               Textra::idx({0,2},{0,1}))
//                    .contract(r_evn,               Textra::idx({0,1},{0,1})) ;
//
//    Eigen::Tensor<Scalar,0> E2ABA_2  =
//            thetaABA
//                    .contract(h0,                  Textra::idx({0,2},{0,1}))
//                    .contract(h1,                  Textra::idx({4,1},{0,1}))
//                    .contract(thetaABA.conjugate(),Textra::idx({2,3,4},{0,2,3}))
//                    .contract(l_evn,               Textra::idx({0,2},{0,1}))
//                    .contract(r_odd,               Textra::idx({0,1},{0,1})) ;
//
//    Eigen::Tensor<Scalar,0> E2BAB_2  =
//            thetaBAB
//                    .contract(h0                  ,Textra::idx({2,3},{0,1}))
//                    .contract(h1                  ,Textra::idx({0,3},{0,1}))
//                    .contract(thetaBAB.conjugate(),Textra::idx({3,4,2},{0,2,3}))
//                    .contract(l_odd               ,Textra::idx({0,2},{0,1}))
//                    .contract(r_evn               ,Textra::idx({0,1},{0,1})) ;
//
//
//    Eigen::Tensor<Scalar,2> E2d_L_evn =
//            theta_evn_normalized
//                    .contract(h0                              ,Textra::idx({0, 2}, {0, 1}))
//                    .contract(theta_evn_normalized.conjugate(),Textra::idx({2, 3}, {0, 2}))
//                    .contract(l_evn                           ,Textra::idx({0, 2}, {0, 1}));
//
//    Eigen::Tensor<Scalar,2> E2d_R_evn =
//            theta_evn_normalized
//                    .contract(h0                              ,Textra::idx({0, 2}, {0, 1}))
//                    .contract(theta_evn_normalized.conjugate(),Textra::idx({2, 3}, {0, 2}))
//                    .contract(r_evn                           ,Textra::idx({1, 3}, {0, 1}));
//
//    Eigen::Tensor<Scalar,2> E2d_L_odd  =
//            theta_odd_normalized
//                    .contract(h1                              , Textra::idx({0, 2}, {0, 1}))
//                    .contract(theta_odd_normalized.conjugate(), Textra::idx({2, 3}, {0, 2}))
//                    .contract(l_odd                           , Textra::idx({0, 2}, {0, 1}));
//
//
//    Eigen::Tensor<Scalar,2> E2d_R_odd =
//            theta_odd_normalized
//                    .contract(h1                              , Textra::idx({0, 2}, {0, 1}))
//                    .contract(theta_odd_normalized.conjugate(), Textra::idx({2, 3}, {0, 2}))
//                    .contract(r_odd                           , Textra::idx({1, 3}, {0, 1}));
//
//    Eigen::array<Eigen::IndexPair<long>,0> pair = {};
//    Eigen::Tensor<Scalar,4> fixpoint_evn = r_evn.contract(l_evn, pair);
//    Eigen::Tensor<Scalar,4> fixpoint_odd = r_odd.contract(l_odd, pair);
//
//    long sizeLA = state.MPS->chiC();
//    long sizeLB = state.MPS->chiB();
//    Eigen::Tensor<Scalar,2> one_minus_transfer_matrix_evn = Matrix_to_Tensor2(MatrixType<Scalar>::Identity(sizeLB*sizeLB, sizeLA*sizeLA).eval()) -
//    (transfer_matrix_evn-fixpoint_evn).reshape(array2{sizeLB*sizeLB, sizeLA*sizeLA}); Eigen::Tensor<Scalar,2> one_minus_transfer_matrix_odd =
//    Matrix_to_Tensor2(MatrixType<Scalar>::Identity(sizeLA*sizeLA, sizeLB*sizeLB).eval()) - (transfer_matrix_odd-fixpoint_odd).reshape(array2{sizeLA*sizeLA,
//    sizeLB*sizeLB}); class_SVD<Scalar> SVD; SVD.setThreshold(settings::precision::svd_threshold);

//    Eigen::Tensor<Scalar,4> E_evn_pinv  = SVD.pseudo_inverse(one_minus_transfer_matrix_evn).reshape(array4{sizeLB,sizeLB,sizeLA,sizeLA});
//    Eigen::Tensor<Scalar,4> E_odd_pinv  = SVD.pseudo_inverse(one_minus_transfer_matrix_odd).reshape(array4{sizeLA,sizeLA,sizeLB,sizeLB});
//    Eigen::Tensor<Scalar,0> E2LRP_ABAB  = E2d_L_evn.contract(E_evn_pinv,idx({0,1},{0,1})).contract(E2d_R_evn,idx({0,1},{0,1}));
//    Eigen::Tensor<Scalar,0> E2LRP_ABBA  = E2d_L_evn.contract(transfer_matrix_LBGA,
//   Textra::idx({0,1},{0,1})).contract(E_odd_pinv,idx({0,1},{0,1})).contract(E2d_R_odd,idx({0,1},{0,1})); Eigen::Tensor<Scalar,0> E2LRP_BABA  =
//    E2d_L_odd.contract(E_odd_pinv,idx({0,1},{0,1})).contract(E2d_R_odd,idx({0,1},{0,1})); Eigen::Tensor<Scalar,0> E2LRP_BAAB  =
//    E2d_L_odd.contract(transfer_matrix_LAGB,Textra::idx({0,1},{0,1})).contract(E_evn_pinv,idx({0,1},{0,1})).contract(E2d_R_evn,idx({0,1},{0,1}));
//
//
//    Scalar e2ab           = E2AB(0);
//    Scalar e2ba           = E2BA(0);
//    Scalar e2aba_1        = E2ABA_1(0);
//    Scalar e2bab_1        = E2BAB_1(0);
//    Scalar e2aba_2        = E2ABA_2(0);
//    Scalar e2bab_2        = E2BAB_2(0);
//    Scalar e2lrpabab      = E2LRP_ABAB(0);
//    Scalar e2lrpabba      = E2LRP_ABBA(0);
//    Scalar e2lrpbaba      = E2LRP_BABA(0);
//    Scalar e2lrpbaab      = E2LRP_BAAB(0);
//
//    return std::real(0.5*(e2ab + e2ba) + 0.5*(e2aba_1  + e2bab_1  + e2aba_2  + e2bab_2 )  + e2lrpabab + e2lrpabba + e2lrpbaba  + e2lrpbaab) ;
////    t_var_ham->toc();
//}
//
//
//
// double tools::infinite::measure::energy_variance_per_site_mom(const class_state_infinite &state){
////    t_var_gen->tic();
//    Scalar a  = (0.0 + 1.0i) *5e-3;
//    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
//    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
//    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
//    auto h_evn = state.HA->single_site_hamiltonian(0,2,SX,SY, SZ);
//    auto h_odd = state.HB->single_site_hamiltonian(1,2,SX,SY, SZ);
//    auto Op_vec = qm::timeEvolution::compute_G(a,4, h_evn, h_odd);
//
//
//    //The following only works if state.MPS has been normalized! I.e, you have to have run MPS->compute_mps_components() prior.
//    Scalar lambdaG  = moment_generating_function(*state.MPS, Op_vec);
//    Scalar l        = 2.0; //Number of sites in unit cell
//    Scalar G        = pow(lambdaG,1.0/l);
//    Scalar logG     = log(lambdaG) * 1.0/l;
//    Scalar logGc    = log(conj(lambdaG) ) * 1.0/l;
//    Scalar O        = (logG - logGc)/(2.0*a);
//    Scalar VarO     = 2.0*log(abs(G))/ (a*a);
////    return std::real(O);
//    return  real(VarO);
////    t_var_gen->toc();
//
//}
//
//

class_state_infinite::Scalar moment_generating_function(const class_state_infinite &state_original, std::vector<Eigen::Tensor<class_state_infinite::Scalar, 2>> &Op_vec) {
    using Scalar                                 = class_state_infinite::Scalar;
    class_state_infinite state_evolved           = state_original;

    class_SVD SVD;
    SVD.setThreshold(settings::precision::svd_threshold);

//    long cfg_chi_lim_max = 5 * state_evolved.chiC();
    for(auto &Op : Op_vec) {
        // Evolve
        Eigen::Tensor<Scalar, 3> mps_evo   = Op.contract(state_evolved.get_2site_tensor(), Textra::idx({0}, {0}));
        state_evolved.set_mps(mps_evo);
        if(&Op != &Op_vec.back()) {
            state_evolved.swap_AB();
        }
    }

    long sizeLB = state_evolved.chiB() * state_evolved.chiB();
    // Normalize
    Eigen::Tensor<Scalar, 2> transfer_matrix_theta_evn = tools::common::views::get_transfer_matrix_theta_evn(state_evolved).reshape(Textra::array2{sizeLB, sizeLB});
    using namespace settings::precision;
    using namespace eigutils::eigSetting;
    class_eigsolver solver;
    solver.eigs<Storage::DENSE>(transfer_matrix_theta_evn.data(), (int) sizeLB, 1, (int) eig_max_ncv, NAN, Form::NONSYMMETRIC, Ritz::LM, Side::R, false);

    auto new_theta_evn_normalized = tools::common::views::get_theta_evn(state_evolved, sqrt(solver.solution.get_eigvals<Form::NONSYMMETRIC>()[0]));
    auto old_theta_evn_normalized = tools::common::views::get_theta_evn(state_original);
    long sizeL = new_theta_evn_normalized.dimension(1) * state_original.chiA();
    long sizeR = new_theta_evn_normalized.dimension(3) * state_original.chiB();

    Eigen::Tensor<Scalar, 2> transfer_matrix_G = new_theta_evn_normalized.contract(old_theta_evn_normalized.conjugate(), Textra::idx({0, 2}, {0, 2}))
                                                     .shuffle(Textra::array4{0, 2, 1, 3})
                                                     .reshape(Textra::array2{sizeL, sizeR});
    // Compute the characteristic function G(a).
    solver.eigs<Storage::DENSE>(transfer_matrix_G.data(), (int) transfer_matrix_G.dimension(0), 1, (int) eig_max_ncv, NAN, Form::NONSYMMETRIC, Ritz::LM,
                                Side::R, false);
    //    solver.eig(transfer_matrix_G.data(),(int)transfer_matrix_G.dimension(0), 1, eig_max_ncv, Ritz::LM, Side::R, false);
    Scalar lambdaG = solver.solution.get_eigvals<Form::NONSYMMETRIC>()[0];
    //    t_temp1->toc();
    return lambdaG;
}

size_t tools::infinite::measure::length(const class_tensors_infinite &tensors) { return tensors.edges->get_length(); }
size_t tools::infinite::measure::length(const class_edges_infinite &edges) { return edges.get_length(); }

double tools::infinite::measure::norm(const class_state_infinite &state) {
    if(state.measurements.norm) return state.measurements.norm.value();
    Eigen::Tensor<Scalar, 0> norm = state.get_2site_tensor().contract(state.get_2site_tensor().conjugate(), Textra::idx({0, 1, 2}, {0, 1, 2}));
    return std::abs(norm(0));
}

long tools::infinite::measure::bond_dimension(const class_state_infinite &state) {
    if(state.measurements.bond_dimension) return state.measurements.bond_dimension.value();
    state.measurements.bond_dimension = state.chiC();
    return state.measurements.bond_dimension.value();
}

double tools::infinite::measure::truncation_error(const class_state_infinite &state) {
    if(state.measurements.truncation_error) return state.measurements.truncation_error.value();
    state.measurements.truncation_error = state.get_truncation_error();
    return state.measurements.truncation_error.value();
}

double tools::infinite::measure::entanglement_entropy(const class_state_infinite &state) {
    tools::common::profile::t_ent->tic();
    if(state.measurements.entanglement_entropy) return state.measurements.entanglement_entropy.value();
    const auto &                   LC = state.LC();
    Eigen::Tensor<Scalar, 0> SA = -LC.square().contract(LC.square().log().eval(), Textra::idx({0}, {0}));
    tools::common::profile::t_ent->toc();
    state.measurements.entanglement_entropy = std::real(SA(0));
    return state.measurements.entanglement_entropy.value();
}

template<typename state_or_mps_type>
double tools::infinite::measure::energy_minus_energy_reduced(const state_or_mps_type &state, const class_model_infinite &model,
                                                             const class_edges_infinite &edges) {
    if constexpr(std::is_same_v<state_or_mps_type, class_state_infinite>) {
        return tools::infinite::measure::energy_minus_energy_reduced(state.get_2site_tensor(), model, edges);
    } else {
        tools::log->trace("Measuring energy mpo");
        const auto &mpo = model.get_2site_tensor();
        const auto &env = edges.get_ene_blk();
        tools::common::profile::t_ene->tic();
        double e_minus_ered = tools::common::moments::first(state, mpo, env.L, env.R);
        tools::common::profile::t_ene->toc();
        return e_minus_ered;
    }
}

template double tools::infinite::measure::energy_minus_energy_reduced(const class_state_infinite &, const class_model_infinite &model,
                                                                      const class_edges_infinite &edges);
template double tools::infinite::measure::energy_minus_energy_reduced(const Eigen::Tensor<Scalar, 3> &, const class_model_infinite &model,
                                                                      const class_edges_infinite &edges);

template<typename state_or_mps_type>
double tools::infinite::measure::energy_mpo(const state_or_mps_type &state, const class_model_infinite &model, const class_edges_infinite &edges) {
    if constexpr(std::is_same_v<state_or_mps_type, class_state_infinite>)
        return tools::infinite::measure::energy_mpo(state.get_2site_tensor(), model, edges);
    else
        return tools::infinite::measure::energy_minus_energy_reduced(state, model, edges) + model.get_energy_per_site_reduced() * static_cast<double>(edges.get_length());
}

template double tools::infinite::measure::energy_mpo(const class_state_infinite &, const class_model_infinite &model,
                                                                      const class_edges_infinite &edges);
template double tools::infinite::measure::energy_mpo(const Eigen::Tensor<Scalar, 3> &, const class_model_infinite &model,
                                                                      const class_edges_infinite &edges);


template<typename state_or_mps_type>
double tools::infinite::measure::energy_per_site_mpo(const state_or_mps_type &state, const class_model_infinite &model, const class_edges_infinite &edges) {
    std::cerr << "energy_per_site_mpo: CHECK DIVISION" << std::endl;
    return tools::infinite::measure::energy_mpo(state,model,edges) /  static_cast<double>(2);
}

template double tools::infinite::measure::energy_per_site_mpo(const class_state_infinite &, const class_model_infinite &model,
                                                     const class_edges_infinite &edges);
template double tools::infinite::measure::energy_per_site_mpo(const Eigen::Tensor<Scalar, 3> &, const class_model_infinite &model,
                                                     const class_edges_infinite &edges);



template<typename state_or_mps_type>
double tools::infinite::measure::energy_variance_mpo(const state_or_mps_type &state, const class_model_infinite &model, const class_edges_infinite &edges) {
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
    if constexpr(std::is_same_v<state_or_mps_type, class_state_infinite>) {
        auto var = tools::infinite::measure::energy_variance_mpo(state.get_2site_tensor(), model, edges);
        if(var < state.lowest_recorded_variance) state.lowest_recorded_variance = var;
        return var;
    }else{
        tools::log->trace("Measuring energy variance mpo");
        double energy = 0;
        if(model.is_reduced())
            energy = tools::infinite::measure::energy_minus_energy_reduced(state, model, edges);
        else
            energy = tools::infinite::measure::energy_mpo(state, model, edges);
        double E2 = energy * energy;
        const auto &mpo = model.get_2site_tensor();
        const auto &env = edges.get_var_blk();
        tools::log->trace("Measuring energy variance mpo");
        tools::common::profile::t_ene->tic();
        double H2 = tools::common::moments::second(state, mpo, env.L, env.R);
        tools::common::profile::t_ene->toc();
        double var = std::abs(H2 - E2);
        return var;
    }
}

template double tools::infinite::measure::energy_variance_mpo(const class_state_infinite &, const class_model_infinite &model,
                                                     const class_edges_infinite &edges);
template double tools::infinite::measure::energy_variance_mpo(const Eigen::Tensor<Scalar, 3> &, const class_model_infinite &model,
                                                     const class_edges_infinite &edges);

template<typename state_or_mps_type>
double tools::infinite::measure::energy_variance_per_site_mpo(const state_or_mps_type &state, const class_model_infinite &model, const class_edges_infinite &edges) {
    std::cerr << "energy_per_site_mpo: CHECK DIVISION" << std::endl;
    return tools::infinite::measure::energy_variance_mpo(state,model,edges) /  static_cast<double>(2);
}

template double tools::infinite::measure::energy_variance_per_site_mpo(const class_state_infinite &, const class_model_infinite &model,
                                                              const class_edges_infinite &edges);
template double tools::infinite::measure::energy_variance_per_site_mpo(const Eigen::Tensor<Scalar, 3> &, const class_model_infinite &model,
                                                              const class_edges_infinite &edges);

double tools::infinite::measure::energy_mpo(const class_tensors_infinite &tensors) {
    if(tensors.measurements.energy_mpo) return tensors.measurements.energy_mpo.value();
    tensors.measurements.energy_mpo = tools::infinite::measure::energy_mpo(*tensors.state, *tensors.model, *tensors.edges);
    return tensors.measurements.energy_mpo.value();
}

double tools::infinite::measure::energy_per_site_mpo(const class_tensors_infinite &tensors) {
    if(tensors.measurements.energy_per_site_mpo) return tensors.measurements.energy_per_site_mpo.value();
    std::cerr << "energy_per_site_mpo: CHECK DIVISION" << std::endl;
    auto L                                   = tools::infinite::measure::length(tensors);
    tensors.measurements.energy_per_site_mpo = tools::infinite::measure::energy_mpo(tensors) / static_cast<double>(L);
    return tensors.measurements.energy_per_site_mpo.value();
}


double tools::infinite::measure::energy_variance_mpo(const class_tensors_infinite &tensors) {
    if(tensors.measurements.energy_variance_mpo) return tensors.measurements.energy_variance_mpo.value();
    tensors.measurements.energy_variance_mpo = tools::infinite::measure::energy_variance_mpo(*tensors.state, *tensors.model, *tensors.edges);
    return tensors.measurements.energy_variance_mpo.value();
}

double tools::infinite::measure::energy_variance_per_site_mpo(const class_tensors_infinite &tensors) {
    if(tensors.measurements.energy_variance_per_site_mpo) return tensors.measurements.energy_variance_per_site_mpo.value();
    auto L                                   = tools::infinite::measure::length(tensors);
    tensors.measurements.energy_variance_per_site_mpo = tools::infinite::measure::energy_variance_mpo(tensors) / static_cast<double>(L);
    return tensors.measurements.energy_variance_per_site_mpo.value();
}

double tools::infinite::measure::energy_mpo(const Eigen::Tensor<Scalar, 3> &mps, const class_tensors_infinite &tensors) {
    return tools::infinite::measure::energy_mpo(mps, *tensors.model, *tensors.edges);
}
double tools::infinite::measure::energy_per_site_mpo(const Eigen::Tensor<Scalar, 3> &mps, const class_tensors_infinite &tensors) {
    return tools::infinite::measure::energy_per_site_mpo(mps, *tensors.model, *tensors.edges);
}
double tools::infinite::measure::energy_variance_mpo(const Eigen::Tensor<Scalar, 3> &mps, const class_tensors_infinite &tensors) {
    return tools::infinite::measure::energy_variance_mpo(mps, *tensors.model, *tensors.edges);
}
double tools::infinite::measure::energy_variance_per_site_mpo(const Eigen::Tensor<Scalar, 3> &mps, const class_tensors_infinite &tensors) {
    return tools::infinite::measure::energy_variance_per_site_mpo(mps, *tensors.model, *tensors.edges);
}



double tools::infinite::measure::energy_per_site_ham(const class_tensors_infinite &tensors) {
    const auto & state = *tensors.state;
    const auto & model = *tensors.model;

    if(tensors.measurements.energy_per_site_ham) return tensors.measurements.energy_per_site_ham.value();
    if(state.measurements.bond_dimension <= 2) return std::numeric_limits<double>::quiet_NaN();
    if(state.chiA() != state.chiB()) return std::numeric_limits<double>::quiet_NaN();
    if(state.chiA() != state.chiC()) return std::numeric_limits<double>::quiet_NaN();
    if(state.chiB() != state.chiC()) return std::numeric_limits<double>::quiet_NaN();
    tools::log->trace("Measuring energy ham");
    tools::common::profile::t_ene_ham->tic();
    auto SX    = qm::gen_manybody_spin(qm::spinOneHalf::sx, 2);
    auto SY    = qm::gen_manybody_spin(qm::spinOneHalf::sy, 2);
    auto SZ    = qm::gen_manybody_spin(qm::spinOneHalf::sz, 2);
    auto h_evn = model.get_mpo_siteA().single_site_hamiltonian(0, 2, SX, SY, SZ);
    auto h_odd = model.get_mpo_siteB().single_site_hamiltonian(1, 2, SX, SY, SZ);
    tools::common::views::compute_mps_components(state);
    using namespace tools::common::views;

    Eigen::Tensor<Scalar, 0> E_evn = theta_evn_normalized.contract(Textra::MatrixTensorMap(h_evn, 2, 2, 2, 2), Textra::idx({0, 2}, {0, 1}))
        .contract(theta_evn_normalized.conjugate(), Textra::idx({2, 3}, {0, 2}))
        .contract(l_evn, Textra::idx({0, 2}, {0, 1}))
        .contract(r_evn, Textra::idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar, 0> E_odd = theta_odd_normalized.contract(Textra::MatrixTensorMap(h_odd, 2, 2, 2, 2), Textra::idx({0, 2}, {0, 1}))
        .contract(theta_odd_normalized.conjugate(), Textra::idx({2, 3}, {0, 2}))
        .contract(l_odd, Textra::idx({0, 2}, {0, 1}))
        .contract(r_odd, Textra::idx({0, 1}, {0, 1}));
    assert(abs(imag(E_evn(0) + E_odd(0))) < 1e-10 and "Energy has an imaginary part!!!");
    tools::common::profile::t_ene_ham->toc();
    tensors.measurements.energy_per_site_ham = 0.5 * std::real(E_evn(0) + E_odd(0));
    return tensors.measurements.energy_per_site_ham.value();
}

double tools::infinite::measure::energy_per_site_mom(const class_tensors_infinite &tensors) {
    if(tensors.measurements.energy_per_site_mom) return tensors.measurements.energy_per_site_mom.value();
    const auto & state = *tensors.state;
    const auto & model = *tensors.model;
    if(state.chiC() <= 2) {
        tensors.measurements.energy_per_site_mom          = std::numeric_limits<double>::quiet_NaN();
        tensors.measurements.energy_variance_per_site_mom = std::numeric_limits<double>::quiet_NaN();
        return tensors.measurements.energy_per_site_mom.value();
    }
    tools::log->trace("Measuring energy mom");

    tools::common::profile::t_ene_mom->tic();
    Scalar a      = Scalar(0.0, 1.0) * 5e-3;
    auto   SX     = qm::gen_manybody_spin(qm::spinOneHalf::sx, 2);
    auto   SY     = qm::gen_manybody_spin(qm::spinOneHalf::sy, 2);
    auto   SZ     = qm::gen_manybody_spin(qm::spinOneHalf::sz, 2);
    auto   h_evn  = model.get_mpo_siteA().single_site_hamiltonian(0, 2, SX, SY, SZ);
    auto   h_odd  = model.get_mpo_siteB().single_site_hamiltonian(1, 2, SX, SY, SZ);
    auto   Op_vec = qm::timeEvolution::compute_G(a, 4, h_evn, h_odd);

    // The following only works if state.MPS has been normalized! I.e, you have to have run MPS->compute_mps_components() prior.
    Scalar lambdaG                                  = moment_generating_function(state, Op_vec);
    Scalar l                                        = 2.0; // Number of sites in unit cell
    Scalar G                                        = pow(lambdaG, 1.0 / l);
    Scalar logG                                     = std::log(lambdaG) * 1.0 / l;
    Scalar logGc                                    = std::log(conj(lambdaG)) * 1.0 / l;
    Scalar O                                        = (logG - logGc) / (2.0 * a);
    Scalar VarO                                     = 2.0 * std::log(abs(G)) / (a * a);
    tensors.measurements.energy_per_site_mom          = std::real(O);
    tensors.measurements.energy_variance_per_site_mom = std::real(VarO);
    tools::common::profile::t_ene_mom->toc();
    return tensors.measurements.energy_per_site_mom.value();
}



double tools::infinite::measure::energy_variance_per_site_ham(const class_tensors_infinite &tensors) {
    if(tensors.measurements.energy_variance_per_site_ham) return tensors.measurements.energy_variance_per_site_ham.value();
    //    if(tensors.MPS->chiA() != tensors.MPS->chiB()) return std::numeric_limits<double>::quiet_NaN();
    //    if(tensors.MPS->chiA() != tensors.MPS->chiC()) return std::numeric_limits<double>::quiet_NaN();
    //    if(tensors.MPS->chiB() != tensors.MPS->chiC()) return std::numeric_limits<double>::quiet_NaN();
    if(tensors.state->chiC() <= 2) return std::numeric_limits<double>::quiet_NaN();

    const auto &state = *tensors.state;
    const auto &model = *tensors.model;

    tools::log->trace("Measuring energy variance ham from tensors");

    tools::common::profile::t_var_ham->tic();
    using namespace tools::common::views;

    auto SX    = qm::gen_manybody_spin(qm::spinOneHalf::sx, 2);
    auto SY    = qm::gen_manybody_spin(qm::spinOneHalf::sy, 2);
    auto SZ    = qm::gen_manybody_spin(qm::spinOneHalf::sz, 2);
    auto h_evn = model.get_mpo_siteA().single_site_hamiltonian(0, 2, SX, SY, SZ);
    auto h_odd = model.get_mpo_siteB().single_site_hamiltonian(1, 2, SX, SY, SZ);
    tools::common::views::compute_mps_components(state);

    Eigen::Tensor<Scalar, 0> E_evn = theta_evn_normalized.contract(Textra::MatrixTensorMap(h_evn, 2, 2, 2, 2), Textra::idx({0, 2}, {0, 1}))
                                         .contract(theta_evn_normalized.conjugate(), Textra::idx({2, 3}, {0, 2}))
                                         .contract(l_evn, Textra::idx({0, 2}, {0, 1}))
                                         .contract(r_evn, Textra::idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar, 0> E_odd = theta_odd_normalized.contract(Textra::MatrixTensorMap(h_odd, 2, 2, 2, 2), Textra::idx({0, 2}, {0, 1}))
                                         .contract(theta_odd_normalized.conjugate(), Textra::idx({2, 3}, {0, 2}))
                                         .contract(l_odd, Textra::idx({0, 2}, {0, 1}))
                                         .contract(r_odd, Textra::idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar, 4> h0 = Textra::MatrixTensorMap((h_evn - E_evn(0) * Textra::MatrixType<Scalar>::Identity(4, 4)).eval(), 2, 2, 2, 2);
    Eigen::Tensor<Scalar, 4> h1 = Textra::MatrixTensorMap((h_odd - E_odd(0) * Textra::MatrixType<Scalar>::Identity(4, 4)).eval(), 2, 2, 2, 2);

    Eigen::Tensor<Scalar, 0> E2AB = theta_evn_normalized.contract(h0, Textra::idx({0, 2}, {0, 1}))
                                        .contract(h0, Textra::idx({2, 3}, {0, 1}))
                                        .contract(theta_evn_normalized.conjugate(), Textra::idx({2, 3}, {0, 2}))
                                        .contract(l_evn, Textra::idx({0, 2}, {0, 1}))
                                        .contract(r_evn, Textra::idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar, 0> E2BA = theta_odd_normalized.contract(h1, Textra::idx({0, 2}, {0, 1}))
                                        .contract(h1, Textra::idx({2, 3}, {0, 1}))
                                        .contract(theta_odd_normalized.conjugate(), Textra::idx({2, 3}, {0, 2}))
                                        .contract(l_odd, Textra::idx({0, 2}, {0, 1}))
                                        .contract(r_odd, Textra::idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar, 5> thetaABA = theta_evn_normalized.contract(LAGA, Textra::idx({3}, {1}));
    Eigen::Tensor<Scalar, 5> thetaBAB = theta_odd_normalized.contract(LCGB, Textra::idx({3}, {1}));

    Eigen::Tensor<Scalar, 0> E2ABA_1 = thetaABA.contract(h1, Textra::idx({2, 3}, {0, 1}))
                                           .contract(h0, Textra::idx({0, 3}, {0, 1}))
                                           .contract(thetaABA.conjugate(), Textra::idx({3, 4, 2}, {0, 2, 3}))
                                           .contract(l_evn, Textra::idx({0, 2}, {0, 1}))
                                           .contract(r_odd, Textra::idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar, 0> E2BAB_1 = thetaBAB.contract(h1, Textra::idx({0, 2}, {0, 1}))
                                           .contract(h0, Textra::idx({4, 1}, {0, 1}))
                                           .contract(thetaBAB.conjugate(), Textra::idx({2, 3, 4}, {0, 2, 3}))
                                           .contract(l_odd, Textra::idx({0, 2}, {0, 1}))
                                           .contract(r_evn, Textra::idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar, 0> E2ABA_2 = thetaABA.contract(h0, Textra::idx({0, 2}, {0, 1}))
                                           .contract(h1, Textra::idx({4, 1}, {0, 1}))
                                           .contract(thetaABA.conjugate(), Textra::idx({2, 3, 4}, {0, 2, 3}))
                                           .contract(l_evn, Textra::idx({0, 2}, {0, 1}))
                                           .contract(r_odd, Textra::idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar, 0> E2BAB_2 = thetaBAB.contract(h0, Textra::idx({2, 3}, {0, 1}))
                                           .contract(h1, Textra::idx({0, 3}, {0, 1}))
                                           .contract(thetaBAB.conjugate(), Textra::idx({3, 4, 2}, {0, 2, 3}))
                                           .contract(l_odd, Textra::idx({0, 2}, {0, 1}))
                                           .contract(r_evn, Textra::idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar, 2> E2d_L_evn = theta_evn_normalized.contract(h0, Textra::idx({0, 2}, {0, 1}))
                                             .contract(theta_evn_normalized.conjugate(), Textra::idx({2, 3}, {0, 2}))
                                             .contract(l_evn, Textra::idx({0, 2}, {0, 1}));

    Eigen::Tensor<Scalar, 2> E2d_R_evn = theta_evn_normalized.contract(h0, Textra::idx({0, 2}, {0, 1}))
                                             .contract(theta_evn_normalized.conjugate(), Textra::idx({2, 3}, {0, 2}))
                                             .contract(r_evn, Textra::idx({1, 3}, {0, 1}));

    Eigen::Tensor<Scalar, 2> E2d_L_odd = theta_odd_normalized.contract(h1, Textra::idx({0, 2}, {0, 1}))
                                             .contract(theta_odd_normalized.conjugate(), Textra::idx({2, 3}, {0, 2}))
                                             .contract(l_odd, Textra::idx({0, 2}, {0, 1}));

    Eigen::Tensor<Scalar, 2> E2d_R_odd = theta_odd_normalized.contract(h1, Textra::idx({0, 2}, {0, 1}))
                                             .contract(theta_odd_normalized.conjugate(), Textra::idx({2, 3}, {0, 2}))
                                             .contract(r_odd, Textra::idx({1, 3}, {0, 1}));

    Eigen::array<Eigen::IndexPair<long>, 0> pair         = {};
    Eigen::Tensor<Scalar, 4>                fixpoint_evn = r_evn.contract(l_evn, pair);
    Eigen::Tensor<Scalar, 4>                fixpoint_odd = r_odd.contract(l_odd, pair);

    long                     sizeLA = state.chiC();
    long                     sizeLB = state.chiB();
    Eigen::Tensor<Scalar, 2> one_minus_transfer_matrix_evn =
        Textra::MatrixTensorMap(Textra::MatrixType<Scalar>::Identity(sizeLB * sizeLB, sizeLA * sizeLA).eval()) -
        (transfer_matrix_evn - fixpoint_evn).reshape(Textra::array2{sizeLB * sizeLB, sizeLA * sizeLA});
    Eigen::Tensor<Scalar, 2> one_minus_transfer_matrix_odd =
        Textra::MatrixTensorMap(Textra::MatrixType<Scalar>::Identity(sizeLA * sizeLA, sizeLB * sizeLB).eval()) -
        (transfer_matrix_odd - fixpoint_odd).reshape(Textra::array2{sizeLA * sizeLA, sizeLB * sizeLB});
    class_SVD SVD;
    SVD.setThreshold(settings::precision::svd_threshold);
    Eigen::Tensor<Scalar, 4> E_evn_pinv = SVD.pseudo_inverse(one_minus_transfer_matrix_evn).reshape(Textra::array4{sizeLB, sizeLB, sizeLA, sizeLA});
    Eigen::Tensor<Scalar, 4> E_odd_pinv = SVD.pseudo_inverse(one_minus_transfer_matrix_odd).reshape(Textra::array4{sizeLA, sizeLA, sizeLB, sizeLB});
    Eigen::Tensor<Scalar, 0> E2LRP_ABAB = E2d_L_evn.contract(E_evn_pinv, Textra::idx({0, 1}, {0, 1})).contract(E2d_R_evn, Textra::idx({0, 1}, {0, 1}));
    Eigen::Tensor<Scalar, 0> E2LRP_ABBA = E2d_L_evn.contract(transfer_matrix_LAGA, Textra::idx({0, 1}, {0, 1}))
                                              .contract(E_odd_pinv, Textra::idx({0, 1}, {0, 1}))
                                              .contract(E2d_R_odd, Textra::idx({0, 1}, {0, 1}));
    Eigen::Tensor<Scalar, 0> E2LRP_BABA = E2d_L_odd.contract(E_odd_pinv, Textra::idx({0, 1}, {0, 1})).contract(E2d_R_odd, Textra::idx({0, 1}, {0, 1}));
    Eigen::Tensor<Scalar, 0> E2LRP_BAAB = E2d_L_odd.contract(transfer_matrix_LCGB, Textra::idx({0, 1}, {0, 1}))
                                              .contract(E_evn_pinv, Textra::idx({0, 1}, {0, 1}))
                                              .contract(E2d_R_evn, Textra::idx({0, 1}, {0, 1}));

    Scalar e2ab      = E2AB(0);
    Scalar e2ba      = E2BA(0);
    Scalar e2aba_1   = E2ABA_1(0);
    Scalar e2bab_1   = E2BAB_1(0);
    Scalar e2aba_2   = E2ABA_2(0);
    Scalar e2bab_2   = E2BAB_2(0);
    Scalar e2lrpabab = E2LRP_ABAB(0);
    Scalar e2lrpabba = E2LRP_ABBA(0);
    Scalar e2lrpbaba = E2LRP_BABA(0);
    Scalar e2lrpbaab = E2LRP_BAAB(0);
    tools::common::profile::t_var_ham->toc();
    tensors.measurements.energy_variance_per_site_ham =
        std::real(0.5 * (e2ab + e2ba) + 0.5 * (e2aba_1 + e2bab_1 + e2aba_2 + e2bab_2) + e2lrpabab + e2lrpabba + e2lrpbaba + e2lrpbaab);
    return tensors.measurements.energy_variance_per_site_ham.value();
}

double tools::infinite::measure::energy_variance_per_site_mom(const class_tensors_infinite &tensors) {
    if(tensors.measurements.energy_variance_per_site_mom) return tensors.measurements.energy_variance_per_site_mom.value();
    tensors.measurements.energy_per_site_mom = tools::infinite::measure::energy_per_site_mom(tensors);
    return tensors.measurements.energy_variance_per_site_mom.value();
}
