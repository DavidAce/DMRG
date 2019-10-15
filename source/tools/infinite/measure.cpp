//
// Created by david on 2019-02-01.
//

#include <tools/nmspc_tools.h>
#include <state/class_infinite_state.h>
#include <state/class_mps_site.h>
#include <state/class_mps_2site.h>
#include <state/class_environment.h>
#include <model/class_model_base.h>
#include <general/nmspc_quantum_mechanics.h>
#include <math/class_eigsolver.h>
#include <math/class_svd_wrapper.h>
using Scalar = std::complex<double>;
using namespace Textra;


//
//
//void tools::infinite::measure::do_all_measurements(const class_infinite_state & state)
//{
//
//
//    using namespace tools::infinite::measure;
//    results::length          = length(state);
//    results::energy          = compute_energy_mpo(state);  //This number is needed for variance calculation!
//    results::energy_per_site = compute_energy_mpo(state);
//    results::energy_variance = compute_energy_variance_mpo(state, results::energy);
//}
//
//
//
//
//
//Scalar moment_generating_function(const class_mps_2site &MPS_original,
//                                                     std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec){
////    t_temp1.tic();
//    std::unique_ptr<class_mps_2site> MPS_evolved = std::make_unique<class_mps_2site>(MPS_original);
//
//    class_SVD<Scalar> SVD;
//    SVD.setThreshold(settings::precision::SVDThreshold);

//    long chi_lim = 5*MPS_evolved->chiC();
////    t_temp2.tic();
//    for (auto &Op: Op_vec) {
//        //Evolve
//        Eigen::Tensor<Scalar, 4> theta_evo = Op.contract(MPS_evolved->get_theta(), idx({0, 1}, {0, 2})).shuffle(array4{0, 2, 1, 3});
//        auto[U, S, V] = SVD.schmidt(theta_evo,chi_lim);
//        MPS_evolved->LC = S;
//        Eigen::Tensor<Scalar,3> L_U =  asDiagonalInversed(MPS_evolved->MPS_A->get_L()).contract(U, idx({1}, {1})).shuffle(array3{1, 0, 2});
//        Eigen::Tensor<Scalar,3> V_L =  V.contract(asDiagonalInversed(MPS_evolved->MPS_B->get_L()), idx({2}, {0}));
//        MPS_evolved->MPS_A->set_G(L_U);
//        MPS_evolved->MPS_B->set_G(V_L);
//
//        if (&Op != &Op_vec.back()) {
//            MPS_evolved->swap_AB();
//        }
//    }
////    t_temp2.toc();
//
//    long sizeLB = MPS_evolved->chiB() * MPS_evolved->chiB();
//    //Normalize
////    t_temp3.tic();
//    Eigen::Tensor<Scalar,2> transfer_matrix_theta_evn       = tools::common::views::get_transfer_matrix_theta_evn(*MPS_evolved).reshape(array2{sizeLB,sizeLB});
////    t_temp3.toc();
//    using namespace settings::precision;
//
////    t_temp4.tic();
//    class_eigsolver_arpack<Scalar, Form::GENERAL> solver;
//    solver.eig(transfer_matrix_theta_evn.data(),(int)sizeLB, 1, eigMaxNcv, eigsolver_properties::Ritz::LM, eigsolver_properties::Side::R, false);
//    auto new_theta_evn_normalized        = tools::common::views::get_theta_evn(*MPS_evolved, sqrt(solver.ref_eigvals()[0]));
////    t_temp4.toc();
//    long sizeL = new_theta_evn_normalized.dimension(1) * MPS_original.chiA();// theta_evn_normalized.dimension(1);
//    long sizeR = new_theta_evn_normalized.dimension(3) * MPS_original.chiB();// theta_evn_normalized.dimension(3);
//
//    Eigen::Tensor<Scalar,2> transfer_matrix_G =   new_theta_evn_normalized
//            .contract(tools::common::views::theta_evn_normalized.conjugate(), idx({0,2},{0,2}))
//            .shuffle(array4{0,2,1,3})
//            .reshape(array2{sizeL,sizeR});
//    //Compute the characteristic function G(a).
//    solver.eig(transfer_matrix_G.data(),(int)transfer_matrix_G.dimension(0), 1, eigMaxNcv, Ritz::LM, Side::R, false);
//    Scalar lambdaG = solver.ref_eigvals()[0];
////    t_temp1.toc();
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
//double tools::infinite::measure::energy(const class_infinite_state & state){
////    state.t_ene_mpo.tic();
//    double energy = tools::finite::measure::energy(state);
//    double L      = tools::finite::measure::length(state);
////    state.t_ene_mpo.toc();
//    return energy / L;
//
//}
//
//
//
//double tools::infinite::measure::energy_per_site_ham(const class_infinite_state & state){
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
//            .contract(Matrix_to_Tensor(h_evn,2,2,2,2),  idx({0, 2}, {0, 1}))
//            .contract(theta_evn_normalized.conjugate(), idx({2, 3}, {0, 2}))
//            .contract(l_evn,                            idx({0, 2}, {0, 1}))
//            .contract(r_evn,                            idx({0, 1}, {0, 1}));
//
//    Eigen::Tensor<Scalar,0>
//    E_odd  = theta_odd_normalized
//            .contract(Matrix_to_Tensor(h_odd,2,2,2,2) ,idx({0, 2}, {0, 1}))
//            .contract(theta_odd_normalized.conjugate(),idx({2, 3}, {0, 2}))
//            .contract(l_odd,                           idx({0, 2}, {0, 1}))
//            .contract(r_odd,                           idx({0, 1}, {0, 1}));
//    assert(abs(imag(E_evn(0)+ E_odd(0))) < 1e-10 and "Energy has an imaginary part!!!" );
////    t_ene_ham.toc();
//    return 0.5*std::real(E_evn(0) + E_odd(0));
//
//}
//
//
//double tools::infinite::measure::energy_per_site_mom(const class_infinite_state & state){
////    t_var_gen.tic();
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
////    t_var_gen.toc();
//
//}
//
//
//double tools::infinite::measure::energy_variance_per_site(const class_infinite_state &state) {
//
//    double VarE  = tools::finite::measure::energy_variance(state);
//    double L     = tools::finite::measure::length(state);
//    return VarE/L;
//}
//
//double tools::infinite::measure::energy_variance_per_site_ham(const class_infinite_state &state) {
////    t_var_ham.tic();
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
//            .contract(Matrix_to_Tensor(h_evn,2,2,2,2),  idx({0, 2}, {0, 1}))
//            .contract(theta_evn_normalized.conjugate(), idx({2, 3}, {0, 2}))
//            .contract(l_evn,                            idx({0, 2}, {0, 1}))
//            .contract(r_evn,                            idx({0, 1}, {0, 1}));
//
//    Eigen::Tensor<Scalar,0>
//            E_odd  = theta_odd_normalized
//            .contract(Matrix_to_Tensor(h_odd,2,2,2,2) ,idx({0, 2}, {0, 1}))
//            .contract(theta_odd_normalized.conjugate(),idx({2, 3}, {0, 2}))
//            .contract(l_odd,                           idx({0, 2}, {0, 1}))
//            .contract(r_odd,                           idx({0, 1}, {0, 1}));
//
//    Eigen::Tensor<Scalar,4> h0 =  Matrix_to_Tensor((h_evn - E_evn(0)*MatrixType<Scalar>::Identity(4,4)).eval(), 2,2,2,2);
//    Eigen::Tensor<Scalar,4> h1 =  Matrix_to_Tensor((h_odd - E_odd(0)*MatrixType<Scalar>::Identity(4,4)).eval(), 2,2,2,2);
//
//    Eigen::Tensor<Scalar,0> E2AB =
//            theta_evn_normalized
//                    .contract(h0                                ,  idx({0, 2}, {0, 1}))
//                    .contract(h0                                ,  idx({2, 3}, {0, 1}))
//                    .contract(theta_evn_normalized.conjugate()  ,  idx({2, 3}, {0, 2}))
//                    .contract(l_evn                             ,  idx({0, 2}, {0, 1}))
//                    .contract(r_evn                             ,  idx({0, 1}, {0, 1}));
//
//
//    Eigen::Tensor<Scalar, 0> E2BA =
//            theta_odd_normalized
//                    .contract(h1                              , idx({0, 2}, {0, 1}))
//                    .contract(h1                              , idx({2, 3}, {0, 1}))
//                    .contract(theta_odd_normalized.conjugate(), idx({2, 3}, {0, 2}))
//                    .contract(l_odd                           , idx({0, 2}, {0, 1}))
//                    .contract(r_odd                           , idx({0, 1}, {0, 1}));
//
//
//
//    Eigen::Tensor<Scalar,5> thetaABA = theta_evn_normalized.contract(LBGA, idx({3},{1}));
//    Eigen::Tensor<Scalar,5> thetaBAB = theta_odd_normalized.contract(LAGB, idx({3},{1}));
//
//    Eigen::Tensor<Scalar,0> E2ABA_1  =
//            thetaABA
//                    .contract(h1,                   idx({2,3},{0,1}))
//                    .contract(h0,                   idx({0,3},{0,1}))
//                    .contract(thetaABA.conjugate(), idx({3,4,2},{0,2,3}))
//                    .contract(l_evn,                idx({0,2},{0,1}))
//                    .contract(r_odd,                idx({0,1},{0,1})) ;
//
//    Eigen::Tensor<Scalar,0> E2BAB_1  =
//            thetaBAB
//                    .contract(h1,                   idx({0,2},{0,1}))
//                    .contract(h0,                   idx({4,1},{0,1}))
//                    .contract(thetaBAB.conjugate(), idx({2,3,4},{0,2,3}))
//                    .contract(l_odd,                idx({0,2},{0,1}))
//                    .contract(r_evn,                idx({0,1},{0,1})) ;
//
//    Eigen::Tensor<Scalar,0> E2ABA_2  =
//            thetaABA
//                    .contract(h0,                   idx({0,2},{0,1}))
//                    .contract(h1,                   idx({4,1},{0,1}))
//                    .contract(thetaABA.conjugate(), idx({2,3,4},{0,2,3}))
//                    .contract(l_evn,                idx({0,2},{0,1}))
//                    .contract(r_odd,                idx({0,1},{0,1})) ;
//
//    Eigen::Tensor<Scalar,0> E2BAB_2  =
//            thetaBAB
//                    .contract(h0                  , idx({2,3},{0,1}))
//                    .contract(h1                  , idx({0,3},{0,1}))
//                    .contract(thetaBAB.conjugate(), idx({3,4,2},{0,2,3}))
//                    .contract(l_odd               , idx({0,2},{0,1}))
//                    .contract(r_evn               , idx({0,1},{0,1})) ;
//
//
//    Eigen::Tensor<Scalar,2> E2d_L_evn =
//            theta_evn_normalized
//                    .contract(h0                              , idx({0, 2}, {0, 1}))
//                    .contract(theta_evn_normalized.conjugate(), idx({2, 3}, {0, 2}))
//                    .contract(l_evn                           , idx({0, 2}, {0, 1}));
//
//    Eigen::Tensor<Scalar,2> E2d_R_evn =
//            theta_evn_normalized
//                    .contract(h0                              , idx({0, 2}, {0, 1}))
//                    .contract(theta_evn_normalized.conjugate(), idx({2, 3}, {0, 2}))
//                    .contract(r_evn                           , idx({1, 3}, {0, 1}));
//
//    Eigen::Tensor<Scalar,2> E2d_L_odd  =
//            theta_odd_normalized
//                    .contract(h1                              ,  idx({0, 2}, {0, 1}))
//                    .contract(theta_odd_normalized.conjugate(),  idx({2, 3}, {0, 2}))
//                    .contract(l_odd                           ,  idx({0, 2}, {0, 1}));
//
//
//    Eigen::Tensor<Scalar,2> E2d_R_odd =
//            theta_odd_normalized
//                    .contract(h1                              ,  idx({0, 2}, {0, 1}))
//                    .contract(theta_odd_normalized.conjugate(),  idx({2, 3}, {0, 2}))
//                    .contract(r_odd                           ,  idx({1, 3}, {0, 1}));
//
//    Eigen::array<Eigen::IndexPair<long>,0> pair = {};
//    Eigen::Tensor<Scalar,4> fixpoint_evn = r_evn.contract(l_evn, pair);
//    Eigen::Tensor<Scalar,4> fixpoint_odd = r_odd.contract(l_odd, pair);
//
//    long sizeLA = state.MPS->chiC();
//    long sizeLB = state.MPS->chiB();
//    Eigen::Tensor<Scalar,2> one_minus_transfer_matrix_evn = Matrix_to_Tensor2(MatrixType<Scalar>::Identity(sizeLB*sizeLB, sizeLA*sizeLA).eval()) - (transfer_matrix_evn-fixpoint_evn).reshape(array2{sizeLB*sizeLB, sizeLA*sizeLA});
//    Eigen::Tensor<Scalar,2> one_minus_transfer_matrix_odd = Matrix_to_Tensor2(MatrixType<Scalar>::Identity(sizeLA*sizeLA, sizeLB*sizeLB).eval()) - (transfer_matrix_odd-fixpoint_odd).reshape(array2{sizeLA*sizeLA, sizeLB*sizeLB});
//    class_SVD<Scalar> SVD;
//    SVD.setThreshold(settings::precision::SVDThreshold);

//    Eigen::Tensor<Scalar,4> E_evn_pinv  = SVD.pseudo_inverse(one_minus_transfer_matrix_evn).reshape(array4{sizeLB,sizeLB,sizeLA,sizeLA});
//    Eigen::Tensor<Scalar,4> E_odd_pinv  = SVD.pseudo_inverse(one_minus_transfer_matrix_odd).reshape(array4{sizeLA,sizeLA,sizeLB,sizeLB});
//    Eigen::Tensor<Scalar,0> E2LRP_ABAB  = E2d_L_evn.contract(E_evn_pinv,idx({0,1},{0,1})).contract(E2d_R_evn,idx({0,1},{0,1}));
//    Eigen::Tensor<Scalar,0> E2LRP_ABBA  = E2d_L_evn.contract(transfer_matrix_LBGA, idx({0,1},{0,1})).contract(E_odd_pinv,idx({0,1},{0,1})).contract(E2d_R_odd,idx({0,1},{0,1}));
//    Eigen::Tensor<Scalar,0> E2LRP_BABA  = E2d_L_odd.contract(E_odd_pinv,idx({0,1},{0,1})).contract(E2d_R_odd,idx({0,1},{0,1}));
//    Eigen::Tensor<Scalar,0> E2LRP_BAAB  = E2d_L_odd.contract(transfer_matrix_LAGB, idx({0,1},{0,1})).contract(E_evn_pinv,idx({0,1},{0,1})).contract(E2d_R_evn,idx({0,1},{0,1}));
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
////    t_var_ham.toc();
//}
//
//
//
//double tools::infinite::measure::energy_variance_per_site_mom(const class_infinite_state &state){
////    t_var_gen.tic();
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
////    t_var_gen.toc();
//
//}
//
//




Scalar moment_generating_function(const class_mps_2site &MPS_original,
                                  std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec){
//    t_temp1.tic();
    std::unique_ptr<class_mps_2site> MPS_evolved = std::make_unique<class_mps_2site>(MPS_original);

    class_SVD SVD;
    SVD.setThreshold(settings::precision::SVDThreshold);

    long chi_max = 5*MPS_evolved->chiC();
//    t_temp2.tic();
    for (auto &Op: Op_vec) {
        //Evolve
        Eigen::Tensor<Scalar, 4> theta_evo = Op.contract(MPS_evolved->get_theta(), idx({0, 1}, {0, 2})).shuffle(array4{0, 2, 1, 3});
        auto[U, S, V] = SVD.schmidt(theta_evo,chi_max);
        MPS_evolved->MPS_A->set_LC(S);
        MPS_evolved->MPS_A->set_M(U);
        MPS_evolved->MPS_B->set_M(V);

        if (&Op != &Op_vec.back()) {
            MPS_evolved->swap_AB();
        }
    }
//    t_temp2.toc();

    long sizeLB = MPS_evolved->chiB() * MPS_evolved->chiB();
    //Normalize
//    t_temp3.tic();
    Eigen::Tensor<Scalar,2> transfer_matrix_theta_evn       = tools::common::views::get_transfer_matrix_theta_evn(*MPS_evolved).reshape(array2{sizeLB,sizeLB});
//    t_temp3.toc();
    using namespace settings::precision;
    using namespace eigutils::eigSetting;
//    t_temp4.tic();
//    class_eigsolver_arpack<Scalar, Form::GENERAL> solver;
    class_eigsolver solver;
    solver.eigs<Storage::DENSE>(transfer_matrix_theta_evn.data(),(int)sizeLB, 1, eigMaxNcv,NAN,Form::NONSYMMETRIC,Ritz::LM,Side::R,false);

//    solver.eig(transfer_matrix_theta_evn.data(),(int)sizeLB, 1, eigMaxNcv, eigsolver_properties::Ritz::LM, eigsolver_properties::Side::R, false);
    auto new_theta_evn_normalized        = tools::common::views::get_theta_evn(*MPS_evolved, sqrt(solver.solution.get_eigvals<Form::NONSYMMETRIC>()[0]));
//    t_temp4.toc();
    long sizeL = new_theta_evn_normalized.dimension(1) * MPS_original.chiA();// theta_evn_normalized.dimension(1);
    long sizeR = new_theta_evn_normalized.dimension(3) * MPS_original.chiB();// theta_evn_normalized.dimension(3);

    Eigen::Tensor<Scalar,2> transfer_matrix_G =   new_theta_evn_normalized
            .contract(tools::common::views::theta_evn_normalized.conjugate(), idx({0,2},{0,2}))
            .shuffle(array4{0,2,1,3})
            .reshape(array2{sizeL,sizeR});
    //Compute the characteristic function G(a).
    solver.eigs<Storage::DENSE>(transfer_matrix_G.data(),(int)transfer_matrix_G.dimension(0), 1, eigMaxNcv,NAN,Form::NONSYMMETRIC,Ritz::LM,Side::R,false);
//    solver.eig(transfer_matrix_G.data(),(int)transfer_matrix_G.dimension(0), 1, eigMaxNcv, Ritz::LM, Side::R, false);
    Scalar lambdaG = solver.solution.get_eigvals<Form::NONSYMMETRIC>()[0];
//    t_temp1.toc();
    return lambdaG;
}


int tools::infinite::measure::length(const class_infinite_state & state){
    return state.get_length();
}


double tools::infinite::measure::norm(const class_infinite_state & state){
    if(state.measurements.norm) {return state.measurements.norm.value();}
    auto theta = state.get_theta();
    Eigen::Tensor<Scalar, 0> norm =
            theta.contract(theta.conjugate(), idx({1, 3, 0, 2}, {1, 3, 0, 2}));
    return std::abs(norm(0));
}


int tools::infinite::measure::bond_dimension(const class_infinite_state & state){
    if(state.measurements.bond_dimension){return state.measurements.bond_dimension.value();}
    return (int) state.MPS->MPS_A->get_chiR();
}

double tools::infinite::measure::truncation_error(const class_infinite_state & state){
    if(state.measurements.truncation_error){return state.measurements.truncation_error.value();}
    return state.MPS->truncation_error;
}



double tools::infinite::measure::current_entanglement_entropy(const class_infinite_state & state){
    tools::log->trace("Measuring entanglement entropy from state");
    tools::common::profile::t_ent.tic();
    if(state.measurements.current_entanglement_entropy){return state.measurements.current_entanglement_entropy.value();}
    auto & LC = state.MPS->MPS_A->get_LC();
    Eigen::Tensor<Scalar,0> SA  = -LC.square()
            .contract(LC.square().log().eval(), idx({0},{0}));
    tools::common::profile::t_ent.toc();
    return std::real(SA(0));
}


double tools::infinite::measure::energy_mpo(const class_infinite_state & state, const Eigen::Tensor<Scalar,4> &theta){
    tools::log->trace("Measuring energy mpo from state");
    tools::common::profile::t_ene_mpo.tic();
    Eigen::Tensor<Scalar, 0>  E =
            state.Lblock->block
                    .contract(theta,                                     idx({0},{1}))
                    .contract(state.HA->MPO(),                      idx({1,2},{0,2}))
                    .contract(state.HB->MPO(),                      idx({3,1},{0,2}))
                    .contract(theta.conjugate(),                         idx({0,2,4},{1,0,2}))
                    .contract(state.Rblock->block,                  idx({0,2,1},{0,1,2}));
    if(abs(imag(E(0))) > 1e-10 ){
        tools::log->critical(fmt::format("Energy has an imaginary part: {:.16f} + i {:.16f}",std::real(E(0)), std::imag(E(0))));
//        throw std::runtime_error("Energy has an imaginary part: " + std::to_string(std::real(E(0))) + " + i " + std::to_string(std::imag(E(0))));
    }
    assert(abs(imag(E(0))) < 1e-10 and "Energy has an imaginary part");
    tools::common::profile::t_ene_mpo.toc();
    return std::real(E(0)) ;
}


double tools::infinite::measure::energy_mpo(const class_infinite_state & state){
    if(state.measurements.energy_mpo){return state.measurements.energy_mpo.value();}
    if(state.sim_type == SimulationType::iTEBD){return std::numeric_limits<double>::quiet_NaN();}
    auto theta    = tools::common::views::get_theta(state);
    double result = tools::infinite::measure::energy_mpo(state,theta);
    return result ;
}


double tools::infinite::measure::energy_per_site_mpo(const class_infinite_state & state){
    if(state.measurements.energy_per_site_mpo){return state.measurements.energy_per_site_mpo.value();}
    auto L     = tools::infinite::measure::length(state);
    return tools::infinite::measure::energy_mpo(state) / L;
}


double tools::infinite::measure::energy_per_site_ham(const class_infinite_state & state){
    if(state.measurements.energy_per_site_ham){return state.measurements.energy_per_site_ham.value();}
    if (state.sim_type == SimulationType::fDMRG){return std::numeric_limits<double>::quiet_NaN();}
    if (state.sim_type == SimulationType::xDMRG){return std::numeric_limits<double>::quiet_NaN();}
    if (state.measurements.bond_dimension <= 2 ){return std::numeric_limits<double>::quiet_NaN();}

    tools::common::profile::t_ene_ham.tic();
    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
    auto h_evn = state.HA->single_site_hamiltonian(0,2,SX,SY, SZ);
    auto h_odd = state.HB->single_site_hamiltonian(1,2,SX,SY, SZ);
    tools::common::views::compute_mps_components(state);
    using namespace tools::common::views;

    Eigen::Tensor<Scalar,0>
            E_evn = theta_evn_normalized
            .contract(MatrixTensorMap(h_evn,2,2,2,2)  , idx({0, 2}, {0, 1}))
            .contract(theta_evn_normalized.conjugate(), idx({2, 3}, {0, 2}))
            .contract(l_evn,                            idx({0, 2}, {0, 1}))
            .contract(r_evn,                            idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar,0>
            E_odd  = theta_odd_normalized
            .contract(MatrixTensorMap(h_odd,2,2,2,2)  ,idx({0, 2}, {0, 1}))
            .contract(theta_odd_normalized.conjugate(),idx({2, 3}, {0, 2}))
            .contract(l_odd,                           idx({0, 2}, {0, 1}))
            .contract(r_odd,                           idx({0, 1}, {0, 1}));
    assert(abs(imag(E_evn(0)+ E_odd(0))) < 1e-10 and "Energy has an imaginary part!!!" );
    tools::common::profile::t_ene_ham.toc();
    return 0.5*std::real(E_evn(0) + E_odd(0));

}


double tools::infinite::measure::energy_per_site_mom(const class_infinite_state & state){
    if(state.measurements.energy_per_site_mom){return state.measurements.energy_per_site_mom.value();}
    if (state.sim_type == SimulationType::fDMRG){return std::numeric_limits<double>::quiet_NaN();}
    if (state.sim_type == SimulationType::xDMRG){return std::numeric_limits<double>::quiet_NaN();}
    if (state.measurements.bond_dimension <= 2 ){return std::numeric_limits<double>::quiet_NaN();}
    tools::common::profile::t_ene_mom.tic();
    Scalar a  = Scalar(0.0 , 1.0) * 5e-3;
    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
    auto h_evn = state.HA->single_site_hamiltonian(0,2,SX,SY, SZ);
    auto h_odd = state.HB->single_site_hamiltonian(1,2,SX,SY, SZ);
    auto Op_vec = qm::timeEvolution::compute_G(a,4, h_evn, h_odd);


    //The following only works if state.MPS has been normalized! I.e, you have to have run MPS->compute_mps_components() prior.
    Scalar lambdaG  = moment_generating_function(*state.MPS, Op_vec);
    Scalar l        = 2.0; //Number of sites in unit cell
    Scalar G        = pow(lambdaG,1.0/l);
    Scalar logG     = std::log(lambdaG) * 1.0/l;
    Scalar logGc    = std::log(conj(lambdaG) ) * 1.0/l;
    Scalar O        = (logG - logGc)/(2.0*a);
    Scalar VarO     = 2.0*std::log(abs(G))/ (a*a);
    state.measurements.energy_per_site_mom           = std::real(O);
    state.measurements.energy_variance_per_site_mom  = std::real(VarO);
    tools::common::profile::t_ene_mom.toc();
    return std::real(O);
}


double tools::infinite::measure::energy_variance_mpo(const class_infinite_state & state,const Eigen::Tensor<std::complex<double>,4> &theta , double &energy_mpo) {
    if (state.sim_type == SimulationType::iTEBD){return std::numeric_limits<double>::quiet_NaN();}
    tools::log->trace("Measuring energy variance mpo from state");
    tools::common::profile::t_var_mpo.tic();
    Eigen::Tensor<Scalar, 0> H2 =
            state.Lblock2->block
                    .contract(theta              ,               idx({0}  ,{1}))
                    .contract(state.HA->MPO(),              idx({1,3},{0,2}))
                    .contract(state.HB->MPO(),              idx({4,2},{0,2}))
                    .contract(state.HA->MPO(),              idx({1,3},{0,2}))
                    .contract(state.HB->MPO(),              idx({4,3},{0,2}))
                    .contract(theta.conjugate()  ,               idx({0,3,5},{1,0,2}))
                    .contract(state.Rblock2->block,         idx({0,3,1,2},{0,1,2,3}));
    tools::common::profile::t_var_mpo.toc();
    if(abs(imag(H2(0))) > 1e-10 ){
        throw std::runtime_error("H2 has an imaginary part: " + std::to_string(std::real(H2(0))) + " + i " + std::to_string(std::imag(H2(0))));
    }
    return std::abs(H2(0) - energy_mpo*energy_mpo);
}



double tools::infinite::measure::energy_variance_mpo(const class_infinite_state & state,const Eigen::Tensor<std::complex<double>,4> &theta) {
    if (state.sim_type == SimulationType::iTEBD){return std::numeric_limits<double>::quiet_NaN();}
    auto energy_mpo = tools::infinite::measure::energy_mpo(state,theta);
    double result = tools::infinite::measure::energy_variance_mpo(state,theta,energy_mpo);
    return result;
}

double tools::infinite::measure::energy_variance_mpo(const class_infinite_state & state) {
    if(state.measurements.energy_variance_mpo){return state.measurements.energy_variance_mpo.value();}
    if (state.sim_type == SimulationType::iTEBD){return std::numeric_limits<double>::quiet_NaN();}
    auto energy_mpo = tools::infinite::measure::energy_mpo(state);
    auto theta      = tools::common::views::get_theta(state);
    return            tools::infinite::measure::energy_variance_mpo(state,theta,energy_mpo);
}


double tools::infinite::measure::energy_variance_per_site_mpo(const class_infinite_state & state) {
    if(state.measurements.energy_variance_per_site_mpo){return state.measurements.energy_variance_per_site_mpo.value();}
    auto L = tools::infinite::measure::length(state);
    return tools::infinite::measure::energy_variance_mpo(state)/L;
}




double tools::infinite::measure::energy_variance_per_site_ham(const class_infinite_state & state) {
    if(state.measurements.energy_variance_per_site_ham){return state.measurements.energy_variance_per_site_ham.value();}
    if (state.MPS->chiA() != state.MPS->chiB()){return std::numeric_limits<double>::quiet_NaN();}
    if (state.MPS->chiA() != state.MPS->chiC()){return std::numeric_limits<double>::quiet_NaN();}
    if (state.sim_type == SimulationType::fDMRG)    {return std::numeric_limits<double>::quiet_NaN();}
    if (state.sim_type == SimulationType::xDMRG)    {return std::numeric_limits<double>::quiet_NaN();}
    if (state.measurements.bond_dimension <= 2 )    {return std::numeric_limits<double>::quiet_NaN();}


    tools::log->trace("Measuring energy variance ham from state");

    tools::common::profile::t_var_ham.tic();
    using namespace tools::common::views;

    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
    auto h_evn = state.HA->single_site_hamiltonian(0,2,SX,SY, SZ);
    auto h_odd = state.HB->single_site_hamiltonian(1,2,SX,SY, SZ);
    tools::common::views::compute_mps_components(state);

    Eigen::Tensor<Scalar,0>
            E_evn = theta_evn_normalized
            .contract(MatrixTensorMap(h_evn,2,2,2,2),   idx({0, 2}, {0, 1}))
            .contract(theta_evn_normalized.conjugate(), idx({2, 3}, {0, 2}))
            .contract(l_evn,                            idx({0, 2}, {0, 1}))
            .contract(r_evn,                            idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar,0>
            E_odd  = theta_odd_normalized
            .contract(MatrixTensorMap(h_odd,2,2,2,2)  ,idx({0, 2}, {0, 1}))
            .contract(theta_odd_normalized.conjugate(),idx({2, 3}, {0, 2}))
            .contract(l_odd,                           idx({0, 2}, {0, 1}))
            .contract(r_odd,                           idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar,4> h0 =  MatrixTensorMap((h_evn - E_evn(0)*MatrixType<Scalar>::Identity(4,4)).eval(), 2,2,2,2);
    Eigen::Tensor<Scalar,4> h1 =  MatrixTensorMap((h_odd - E_odd(0)*MatrixType<Scalar>::Identity(4,4)).eval(), 2,2,2,2);

    Eigen::Tensor<Scalar,0> E2AB =
            theta_evn_normalized
                    .contract(h0                                ,  idx({0, 2}, {0, 1}))
                    .contract(h0                                ,  idx({2, 3}, {0, 1}))
                    .contract(theta_evn_normalized.conjugate()  ,  idx({2, 3}, {0, 2}))
                    .contract(l_evn                             ,  idx({0, 2}, {0, 1}))
                    .contract(r_evn                             ,  idx({0, 1}, {0, 1}));


    Eigen::Tensor<Scalar, 0> E2BA =
            theta_odd_normalized
                    .contract(h1                              , idx({0, 2}, {0, 1}))
                    .contract(h1                              , idx({2, 3}, {0, 1}))
                    .contract(theta_odd_normalized.conjugate(), idx({2, 3}, {0, 2}))
                    .contract(l_odd                           , idx({0, 2}, {0, 1}))
                    .contract(r_odd                           , idx({0, 1}, {0, 1}));



    Eigen::Tensor<Scalar,5> thetaABA = theta_evn_normalized.contract(LBGA, idx({3},{1}));
    Eigen::Tensor<Scalar,5> thetaBAB = theta_odd_normalized.contract(LAGB, idx({3},{1}));

    Eigen::Tensor<Scalar,0> E2ABA_1  =
            thetaABA
                    .contract(h1,                   idx({2,3},{0,1}))
                    .contract(h0,                   idx({0,3},{0,1}))
                    .contract(thetaABA.conjugate(), idx({3,4,2},{0,2,3}))
                    .contract(l_evn,                idx({0,2},{0,1}))
                    .contract(r_odd,                idx({0,1},{0,1})) ;

    Eigen::Tensor<Scalar,0> E2BAB_1  =
            thetaBAB
                    .contract(h1,                   idx({0,2},{0,1}))
                    .contract(h0,                   idx({4,1},{0,1}))
                    .contract(thetaBAB.conjugate(), idx({2,3,4},{0,2,3}))
                    .contract(l_odd,                idx({0,2},{0,1}))
                    .contract(r_evn,                idx({0,1},{0,1})) ;

    Eigen::Tensor<Scalar,0> E2ABA_2  =
            thetaABA
                    .contract(h0,                   idx({0,2},{0,1}))
                    .contract(h1,                   idx({4,1},{0,1}))
                    .contract(thetaABA.conjugate(), idx({2,3,4},{0,2,3}))
                    .contract(l_evn,                idx({0,2},{0,1}))
                    .contract(r_odd,                idx({0,1},{0,1})) ;

    Eigen::Tensor<Scalar,0> E2BAB_2  =
            thetaBAB
                    .contract(h0                  , idx({2,3},{0,1}))
                    .contract(h1                  , idx({0,3},{0,1}))
                    .contract(thetaBAB.conjugate(), idx({3,4,2},{0,2,3}))
                    .contract(l_odd               , idx({0,2},{0,1}))
                    .contract(r_evn               , idx({0,1},{0,1})) ;


    Eigen::Tensor<Scalar,2> E2d_L_evn =
            theta_evn_normalized
                    .contract(h0                              , idx({0, 2}, {0, 1}))
                    .contract(theta_evn_normalized.conjugate(), idx({2, 3}, {0, 2}))
                    .contract(l_evn                           , idx({0, 2}, {0, 1}));

    Eigen::Tensor<Scalar,2> E2d_R_evn =
            theta_evn_normalized
                    .contract(h0                              , idx({0, 2}, {0, 1}))
                    .contract(theta_evn_normalized.conjugate(), idx({2, 3}, {0, 2}))
                    .contract(r_evn                           , idx({1, 3}, {0, 1}));

    Eigen::Tensor<Scalar,2> E2d_L_odd  =
            theta_odd_normalized
                    .contract(h1                              ,  idx({0, 2}, {0, 1}))
                    .contract(theta_odd_normalized.conjugate(),  idx({2, 3}, {0, 2}))
                    .contract(l_odd                           ,  idx({0, 2}, {0, 1}));


    Eigen::Tensor<Scalar,2> E2d_R_odd =
            theta_odd_normalized
                    .contract(h1                              ,  idx({0, 2}, {0, 1}))
                    .contract(theta_odd_normalized.conjugate(),  idx({2, 3}, {0, 2}))
                    .contract(r_odd                           ,  idx({1, 3}, {0, 1}));

    Eigen::array<Eigen::IndexPair<long>,0> pair = {};
    Eigen::Tensor<Scalar,4> fixpoint_evn = r_evn.contract(l_evn, pair);
    Eigen::Tensor<Scalar,4> fixpoint_odd = r_odd.contract(l_odd, pair);

    long sizeLA = state.MPS->chiC();
    long sizeLB = state.MPS->chiB();
    Eigen::Tensor<Scalar,2> one_minus_transfer_matrix_evn = MatrixTensorMap(MatrixType<Scalar>::Identity(sizeLB*sizeLB, sizeLA*sizeLA).eval()) - (transfer_matrix_evn-fixpoint_evn).reshape(array2{sizeLB*sizeLB, sizeLA*sizeLA});
    Eigen::Tensor<Scalar,2> one_minus_transfer_matrix_odd = MatrixTensorMap(MatrixType<Scalar>::Identity(sizeLA*sizeLA, sizeLB*sizeLB).eval()) - (transfer_matrix_odd-fixpoint_odd).reshape(array2{sizeLA*sizeLA, sizeLB*sizeLB});
    class_SVD SVD;
    SVD.setThreshold(settings::precision::SVDThreshold);
    Eigen::Tensor<Scalar,4> E_evn_pinv  = SVD.pseudo_inverse(one_minus_transfer_matrix_evn).reshape(array4{sizeLB,sizeLB,sizeLA,sizeLA});
    Eigen::Tensor<Scalar,4> E_odd_pinv  = SVD.pseudo_inverse(one_minus_transfer_matrix_odd).reshape(array4{sizeLA,sizeLA,sizeLB,sizeLB});
    Eigen::Tensor<Scalar,0> E2LRP_ABAB  = E2d_L_evn.contract(E_evn_pinv,idx({0,1},{0,1})).contract(E2d_R_evn,idx({0,1},{0,1}));
    Eigen::Tensor<Scalar,0> E2LRP_ABBA  = E2d_L_evn.contract(transfer_matrix_LBGA, idx({0,1},{0,1})).contract(E_odd_pinv,idx({0,1},{0,1})).contract(E2d_R_odd,idx({0,1},{0,1}));
    Eigen::Tensor<Scalar,0> E2LRP_BABA  = E2d_L_odd.contract(E_odd_pinv,idx({0,1},{0,1})).contract(E2d_R_odd,idx({0,1},{0,1}));
    Eigen::Tensor<Scalar,0> E2LRP_BAAB  = E2d_L_odd.contract(transfer_matrix_LAGB, idx({0,1},{0,1})).contract(E_evn_pinv,idx({0,1},{0,1})).contract(E2d_R_evn,idx({0,1},{0,1}));


    Scalar e2ab           = E2AB(0);
    Scalar e2ba           = E2BA(0);
    Scalar e2aba_1        = E2ABA_1(0);
    Scalar e2bab_1        = E2BAB_1(0);
    Scalar e2aba_2        = E2ABA_2(0);
    Scalar e2bab_2        = E2BAB_2(0);
    Scalar e2lrpabab      = E2LRP_ABAB(0);
    Scalar e2lrpabba      = E2LRP_ABBA(0);
    Scalar e2lrpbaba      = E2LRP_BABA(0);
    Scalar e2lrpbaab      = E2LRP_BAAB(0);
    tools::common::profile::t_var_ham.toc();

    return std::real(0.5*(e2ab + e2ba) + 0.5*(e2aba_1  + e2bab_1  + e2aba_2  + e2bab_2 )  + e2lrpabab + e2lrpabba + e2lrpbaba  + e2lrpbaab) ;
}


double tools::infinite::measure::energy_variance_per_site_mom(const class_infinite_state & state){
    if(state.measurements.energy_variance_per_site_mom){return state.measurements.energy_variance_per_site_mom.value();}
    if (state.sim_type == SimulationType::fDMRG)    {return std::numeric_limits<double>::quiet_NaN();}
    if (state.sim_type == SimulationType::xDMRG)    {return std::numeric_limits<double>::quiet_NaN();}
    [[maybe_unused]] auto dummy = energy_per_site_mom(state);
    return state.measurements.energy_variance_per_site_mom.value();
}
