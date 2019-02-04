//
// Created by david on 2019-02-04.
//



#include <mps_routines/nmspc_mps_tools.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mps_2site.h>
#include <mps_routines/class_environment.h>
#include <model/class_hamiltonian_base.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/class_svd_wrapper.h>
#include <general/class_eigsolver_arpack.h>

using Scalar = std::complex<double>;
using namespace eigsolver_properties;
using namespace Textra;


namespace MPS_Tools::Common::Measure::Results {
    size_t length                           = 0;
    size_t bond_dimension                   = 0;
    double norm                             = 0;
    double truncation_error                 = 0;

    double energy_mpo                       = 0;
    double energy_per_site_mpo              = 0;
    double energy_variance_mpo              = 0;
    double energy_variance_per_site_mpo     = 0;

    double energy_ham                       = 0;
    double energy_variance_ham              = 0;

    double energy_mom                       = 0;
    double energy_variance_mom              = 0;

    double midchain_entanglement_entropy    = 0;
    bool   superblock_measured              = false;
}




void MPS_Tools::Common::Measure::do_all_measurements(const class_superblock & superblock)
{
    using namespace MPS_Tools::Common::Measure;
    if (Results::superblock_measured){return;}

    Results::length                         = length(superblock);
    Results::norm                           = norm(superblock);

    Results::energy_mpo                     = energy_mpo(superblock);  //This number is needed for variance calculation!
    Results::energy_per_site_mpo            = Results::energy_mpo/Results::length;
    Results::energy_variance_mpo            = energy_variance_mpo(superblock, Results::energy_mpo);
    Results::energy_variance_per_site_mpo   = Results::energy_variance_mpo / Results::length;

    Results::energy_ham                     = energy_ham(superblock);  //This number is needed for variance calculation!
    Results::energy_variance_ham            = energy_variance_ham(superblock);

//    Results::energy_mom                     = energy_mom(superblock);  //This number is needed for variance calculation!
    Results::energy_variance_mom            = energy_variance_mom(superblock, Results::energy_mom); // This function computes the energy also

    Results::midchain_entanglement_entropy  = midchain_entanglement_entropy(superblock);
    Results::superblock_measured = true;
}

void MPS_Tools::Common::Measure::set_not_measured(){
    MPS_Tools::Common::Measure::Results::superblock_measured = false;
    MPS_Tools::Finite::Measure::Results::state_measured = false;
}



Scalar moment_generating_function(const class_mps_2site &MPS_original,
                                  std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec){
//    t_temp1.tic();
    std::unique_ptr<class_mps_2site> MPS_evolved = std::make_unique<class_mps_2site>(MPS_original);

    class_SVD<Scalar> SVD;
    long chi_max = 5*MPS_evolved->chiC();
//    t_temp2.tic();
    for (auto &Op: Op_vec) {
        //Evolve
        Eigen::Tensor<Scalar, 4> theta_evo = Op.contract(MPS_evolved->get_theta(), idx({0, 1}, {0, 2})).shuffle(array4{0, 2, 1, 3});
        auto[U, S, V] = SVD.schmidt(theta_evo,chi_max);
        MPS_evolved->LC = S;
        Eigen::Tensor<Scalar,3> L_U =  asDiagonalInversed(MPS_evolved->MPS_A->get_L()).contract(U, idx({1}, {1})).shuffle(array3{1, 0, 2});
        Eigen::Tensor<Scalar,3> V_L =  V.contract(asDiagonalInversed(MPS_evolved->MPS_B->get_L()), idx({2}, {0}));
        MPS_evolved->MPS_A->set_G(L_U);
        MPS_evolved->MPS_B->set_G(V_L);

        if (&Op != &Op_vec.back()) {
            MPS_evolved->swap_AB();
        }
    }
//    t_temp2.toc();

    long sizeLB = MPS_evolved->chiB() * MPS_evolved->chiB();
    //Normalize
//    t_temp3.tic();
    Eigen::Tensor<Scalar,2> transfer_matrix_theta_evn       = MPS_Tools::Common::Views::get_transfer_matrix_theta_evn(*MPS_evolved).reshape(array2{sizeLB,sizeLB});
//    t_temp3.toc();
    using namespace settings::precision;

//    t_temp4.tic();
    class_eigsolver_arpack<Scalar, Form::GENERAL> solver;
    solver.eig(transfer_matrix_theta_evn.data(),(int)sizeLB, 1, eigMaxNcv, eigsolver_properties::Ritz::LM, eigsolver_properties::Side::R, false);
    auto new_theta_evn_normalized        = MPS_Tools::Common::Views::get_theta_evn(*MPS_evolved, sqrt(solver.ref_eigvals()[0]));
//    t_temp4.toc();
    long sizeL = new_theta_evn_normalized.dimension(1) * MPS_original.chiA();// theta_evn_normalized.dimension(1);
    long sizeR = new_theta_evn_normalized.dimension(3) * MPS_original.chiB();// theta_evn_normalized.dimension(3);

    Eigen::Tensor<Scalar,2> transfer_matrix_G =   new_theta_evn_normalized
            .contract(MPS_Tools::Common::Views::theta_evn_normalized.conjugate(), idx({0,2},{0,2}))
            .shuffle(array4{0,2,1,3})
            .reshape(array2{sizeL,sizeR});
    //Compute the characteristic function G(a).
    solver.eig(transfer_matrix_G.data(),(int)transfer_matrix_G.dimension(0), 1, eigMaxNcv, Ritz::LM, Side::R, false);
    Scalar lambdaG = solver.ref_eigvals()[0];
//    t_temp1.toc();
    return lambdaG;
}


size_t MPS_Tools::Common::Measure::length(const class_superblock & superblock){
    return superblock.get_length();
}


double MPS_Tools::Common::Measure::norm(const class_superblock & superblock){
    auto theta = MPS_Tools::Common::Views::get_theta(superblock);
    Eigen::Tensor<Scalar, 0> norm =
            theta.contract(theta.conjugate(), idx({1, 3, 0, 2}, {1, 3, 0, 2}));
    return std::real(norm(0));
}


size_t MPS_Tools::Common::Measure::bond_dimension(const class_superblock & superblock){
    return (size_t) superblock.MPS->LC.dimension(0);
}


double MPS_Tools::Common::Measure::midchain_entanglement_entropy(const class_superblock & superblock){
//    t_entropy.tic();
    auto & LC = superblock.MPS->LC;

    Eigen::Tensor<Scalar,0> SA  = -LC.square()
            .contract(LC.square().log().eval(), idx({0},{0}));
    return std::real(SA(0));
//    t_entropy.toc();
}

double MPS_Tools::Common::Measure::energy_mpo(const class_superblock & superblock){
//    superblock.t_ene_mpo.tic();

    auto theta = MPS_Tools::Common::Views::get_theta(superblock);
    Eigen::Tensor<Scalar, 0>  E =
            superblock.Lblock->block
                    .contract(theta,                                     idx({0},{1}))
                    .contract(superblock.HA->MPO,                        idx({1,2},{0,2}))
                    .contract(superblock.HB->MPO,                        idx({3,1},{0,2}))
                    .contract(theta.conjugate(),                         idx({0,2,4},{1,0,2}))
                    .contract(superblock.Rblock->block,                  idx({0,2,1},{0,1,2}));
    assert(abs(imag(E(0))) < 1e-10 and "Energy has an imaginary part!!!");

    return std::real(E(0)) ;
    //    superblock.t_ene_mpo.toc();

}


double MPS_Tools::Common::Measure::energy_ham(const class_superblock & superblock){
    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
    auto h_evn = superblock.HA->single_site_hamiltonian(0,2,SX,SY, SZ);
    auto h_odd = superblock.HB->single_site_hamiltonian(1,2,SX,SY, SZ);
    MPS_Tools::Common::Views::compute_mps_components(superblock);
    using namespace MPS_Tools::Common::Views;

    Eigen::Tensor<Scalar,0>
            E_evn = theta_evn_normalized
            .contract(Matrix_to_Tensor(h_evn,2,2,2,2),  idx({0, 2}, {0, 1}))
            .contract(theta_evn_normalized.conjugate(), idx({2, 3}, {0, 2}))
            .contract(l_evn,                            idx({0, 2}, {0, 1}))
            .contract(r_evn,                            idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar,0>
            E_odd  = theta_odd_normalized
            .contract(Matrix_to_Tensor(h_odd,2,2,2,2) ,idx({0, 2}, {0, 1}))
            .contract(theta_odd_normalized.conjugate(),idx({2, 3}, {0, 2}))
            .contract(l_odd,                           idx({0, 2}, {0, 1}))
            .contract(r_odd,                           idx({0, 1}, {0, 1}));
    assert(abs(imag(E_evn(0)+ E_odd(0))) < 1e-10 and "Energy has an imaginary part!!!" );
//    t_ene_ham.toc();
    return 0.5*std::real(E_evn(0) + E_odd(0));

}


double MPS_Tools::Common::Measure::energy_mom(const class_superblock & superblock){
//    t_var_gen.tic();
    Scalar a  = (0.0 + 1.0i) *5e-3;
    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
    auto h_evn = superblock.HA->single_site_hamiltonian(0,2,SX,SY, SZ);
    auto h_odd = superblock.HB->single_site_hamiltonian(1,2,SX,SY, SZ);
    auto Op_vec = qm::timeEvolution::compute_G(a,4, h_evn, h_odd);


    //The following only works if superblock.MPS has been normalized! I.e, you have to have run MPS->compute_mps_components() prior.
    Scalar lambdaG  = moment_generating_function(*superblock.MPS, Op_vec);
    Scalar l        = 2.0; //Number of sites in unit cell
//    Scalar G        = pow(lambdaG,1.0/l);
    Scalar logG     = log(lambdaG) * 1.0/l;
    Scalar logGc    = log(conj(lambdaG) ) * 1.0/l;
    Scalar O        = (logG - logGc)/(2.0*a);
//    Scalar VarO     = 2.0*log(abs(G))/ (a*a);
    return std::real(O);
//    variance_mom    = real(VarO);
//    t_var_gen.toc();

}


double MPS_Tools::Common::Measure::energy_variance_mpo(const class_superblock &superblock) {
    double energy = MPS_Tools::Common::Measure::energy_mpo(superblock);
    auto theta = MPS_Tools::Common::Views::get_theta(superblock);
    Eigen::Tensor<Scalar, 0> H2 =
            superblock.Lblock2->block
                    .contract(theta              ,               idx({0}  ,{1}))
                    .contract(superblock.HA->MPO,                idx({1,3},{0,2}))
                    .contract(superblock.HB->MPO,                idx({4,2},{0,2}))
                    .contract(superblock.HA->MPO,                idx({1,3},{0,2}))
                    .contract(superblock.HB->MPO,                idx({4,3},{0,2}))
                    .contract(theta.conjugate()  ,               idx({0,3,5},{1,0,2}))
                    .contract(superblock.Rblock2->block,         idx({0,3,1,2},{0,1,2,3}));

    return std::abs(H2(0) - energy*energy);

}


double MPS_Tools::Common::Measure::energy_variance_mpo(const class_superblock &superblock, double & energy_mpo) {
    auto theta = MPS_Tools::Common::Views::get_theta(superblock);
    Eigen::Tensor<Scalar, 0> H2 =
            superblock.Lblock2->block
                    .contract(theta              ,               idx({0}  ,{1}))
                    .contract(superblock.HA->MPO,                idx({1,3},{0,2}))
                    .contract(superblock.HB->MPO,                idx({4,2},{0,2}))
                    .contract(superblock.HA->MPO,                idx({1,3},{0,2}))
                    .contract(superblock.HB->MPO,                idx({4,3},{0,2}))
                    .contract(theta.conjugate()  ,               idx({0,3,5},{1,0,2}))
                    .contract(superblock.Rblock2->block,         idx({0,3,1,2},{0,1,2,3}));

    return std::abs(H2(0) - energy_mpo*energy_mpo);

}


double MPS_Tools::Common::Measure::energy_variance_ham(const class_superblock &superblock) {
//    t_var_ham.tic();
    using namespace MPS_Tools::Common::Views;

    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
    auto h_evn = superblock.HA->single_site_hamiltonian(0,2,SX,SY, SZ);
    auto h_odd = superblock.HB->single_site_hamiltonian(1,2,SX,SY, SZ);
    MPS_Tools::Common::Views::compute_mps_components(superblock);

    Eigen::Tensor<Scalar,0>
            E_evn = theta_evn_normalized
            .contract(Matrix_to_Tensor(h_evn,2,2,2,2),  idx({0, 2}, {0, 1}))
            .contract(theta_evn_normalized.conjugate(), idx({2, 3}, {0, 2}))
            .contract(l_evn,                            idx({0, 2}, {0, 1}))
            .contract(r_evn,                            idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar,0>
            E_odd  = theta_odd_normalized
            .contract(Matrix_to_Tensor(h_odd,2,2,2,2) ,idx({0, 2}, {0, 1}))
            .contract(theta_odd_normalized.conjugate(),idx({2, 3}, {0, 2}))
            .contract(l_odd,                           idx({0, 2}, {0, 1}))
            .contract(r_odd,                           idx({0, 1}, {0, 1}));

    Eigen::Tensor<Scalar,4> h0 =  Matrix_to_Tensor((h_evn - E_evn(0)*MatrixType<Scalar>::Identity(4,4)).eval(), 2,2,2,2);
    Eigen::Tensor<Scalar,4> h1 =  Matrix_to_Tensor((h_odd - E_odd(0)*MatrixType<Scalar>::Identity(4,4)).eval(), 2,2,2,2);

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

    long sizeLA = superblock.MPS->chiC();
    long sizeLB = superblock.MPS->chiB();
    Eigen::Tensor<Scalar,2> one_minus_transfer_matrix_evn = Matrix_to_Tensor2(MatrixType<Scalar>::Identity(sizeLB*sizeLB, sizeLA*sizeLA).eval()) - (transfer_matrix_evn-fixpoint_evn).reshape(array2{sizeLB*sizeLB, sizeLA*sizeLA});
    Eigen::Tensor<Scalar,2> one_minus_transfer_matrix_odd = Matrix_to_Tensor2(MatrixType<Scalar>::Identity(sizeLA*sizeLA, sizeLB*sizeLB).eval()) - (transfer_matrix_odd-fixpoint_odd).reshape(array2{sizeLA*sizeLA, sizeLB*sizeLB});
    class_SVD<Scalar> SVD;
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

    return std::real(0.5*(e2ab + e2ba) + 0.5*(e2aba_1  + e2bab_1  + e2aba_2  + e2bab_2 )  + e2lrpabab + e2lrpabba + e2lrpbaba  + e2lrpbaab) ;
//    t_var_ham.toc();
}


double MPS_Tools::Common::Measure::energy_variance_mom(const class_superblock &superblock){
//    t_var_gen.tic();
    Scalar a  = (0.0 + 1.0i) *5e-3;
    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
    auto h_evn = superblock.HA->single_site_hamiltonian(0,2,SX,SY, SZ);
    auto h_odd = superblock.HB->single_site_hamiltonian(1,2,SX,SY, SZ);
    auto Op_vec = qm::timeEvolution::compute_G(a,4, h_evn, h_odd);


    //The following only works if superblock.MPS has been normalized! I.e, you have to have run MPS->compute_mps_components() prior.
    Scalar lambdaG  = moment_generating_function(*superblock.MPS, Op_vec);
    Scalar l        = 2.0; //Number of sites in unit cell
    Scalar G        = pow(lambdaG,1.0/l);
    Scalar logG     = log(lambdaG) * 1.0/l;
    Scalar logGc    = log(conj(lambdaG) ) * 1.0/l;
    Scalar O        = (logG - logGc)/(2.0*a);
    Scalar VarO     = 2.0*log(abs(G))/ (a*a);
//    return std::real(O);
    Results::energy_mom =  std::real(O);

    return  real(VarO);
//    t_var_gen.toc();

}


double MPS_Tools::Common::Measure::energy_variance_mom(const class_superblock &superblock, double & energy_mom){
//    t_var_gen.tic();
    Scalar a  = (0.0 + 1.0i) *5e-3;
    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
    auto h_evn = superblock.HA->single_site_hamiltonian(0,2,SX,SY, SZ);
    auto h_odd = superblock.HB->single_site_hamiltonian(1,2,SX,SY, SZ);
    auto Op_vec = qm::timeEvolution::compute_G(a,4, h_evn, h_odd);


    //The following only works if superblock.MPS has been normalized! I.e, you have to have run MPS->compute_mps_components() prior.
    Scalar lambdaG  = moment_generating_function(*superblock.MPS, Op_vec);
    Scalar l        = 2.0; //Number of sites in unit cell
    Scalar G        = pow(lambdaG,1.0/l);
    Scalar logG     = log(lambdaG) * 1.0/l;
    Scalar logGc    = log(conj(lambdaG) ) * 1.0/l;
    Scalar O        = (logG - logGc)/(2.0*a);
    Scalar VarO     = 2.0*log(abs(G))/ (a*a);
    energy_mom      =  std::real(O);
    return  real(VarO);
//    t_var_gen.toc();

}

