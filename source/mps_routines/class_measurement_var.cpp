//
// Created by david on 2018-02-21.
//
#include <iomanip>
#include <complex>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_mps_2site.h>
#include <general/nmspc_tensor_extra.h>
#include <general/class_svd_wrapper.h>
#include <general/class_eigsolver_arpack.h>
#include <mps_routines/class_finite_chain_sweeper.h>
#include <mps_routines/class_mps_util.h>
using namespace std;
using namespace Textra;
using Scalar = class_measurement::Scalar;
using namespace std::complex_literals;
using namespace eigsolver_properties;

Scalar class_measurement::moment_generating_function(const std::unique_ptr<class_mps_2site> &MPS_original,
                                                     std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec){
    t_temp1.tic();
    std::unique_ptr<class_mps_2site> MPS_evolved = std::make_unique<class_mps_2site>(*MPS_original);

    class_SVD<Scalar> SVD;
    long chi_max = 5*MPS_evolved->chiC();
    t_temp2.tic();
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
    t_temp2.toc();

    long sizeLB = MPS_evolved->chiB() * MPS_evolved->chiB();
    //Normalize
    t_temp3.tic();
    Eigen::Tensor<Scalar,2> transfer_matrix_theta_evn       = mps_util->get_transfer_matrix_theta_evn(MPS_evolved).reshape(array2{sizeLB,sizeLB});
    t_temp3.toc();
    using namespace settings::precision;

    t_temp4.tic();
    class_eigsolver_arpack<Scalar, Form::GENERAL> solver;
    solver.eig(transfer_matrix_theta_evn.data(),(int)sizeLB, 1, eigMaxNcv, eigsolver_properties::Ritz::LM, eigsolver_properties::Side::R, false);
    auto new_theta_evn_normalized        = mps_util->get_theta_evn(MPS_evolved, sqrt(solver.ref_eigvals()[0]));
    t_temp4.toc();
    long sizeL = new_theta_evn_normalized.dimension(1) * MPS_original->chiA();// theta_evn_normalized.dimension(1);
    long sizeR = new_theta_evn_normalized.dimension(3) * MPS_original->chiB();// theta_evn_normalized.dimension(3);

    Eigen::Tensor<Scalar,2> transfer_matrix_G =   new_theta_evn_normalized
            .contract(mps_util->theta_evn_normalized.conjugate(), idx({0,2},{0,2}))
            .shuffle(array4{0,2,1,3})
            .reshape(array2{sizeL,sizeR});
    //Compute the characteristic function G(a).
    solver.eig(transfer_matrix_G.data(),(int)transfer_matrix_G.dimension(0), 1, eigMaxNcv, Ritz::LM, Side::R, false);
    Scalar lambdaG = solver.ref_eigvals()[0];
    t_temp1.toc();
    return lambdaG;
}


void class_measurement::compute_energy_and_variance_mom(Scalar a,std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec){
    t_var_gen.tic();
    using T = Scalar;
    //The following only works if superblock->MPS has been normalized! I.e, you have to have run MPS->compute_mps_components() prior.
    T lambdaG  = moment_generating_function(superblock->MPS, Op_vec);
    T l        = 2.0; //Number of sites in unit cell
    T G        = pow(lambdaG,1.0/l);
    T logG     = log(lambdaG) * 1.0/l;
    T logGc    = log(conj(lambdaG) ) * 1.0/l;
    T O        = (logG - logGc)/(2.0*a);
    T VarO     = 2.0*log(abs(G))/ (a*a);
    energy_mom = real(O);
    variance_mom = real(VarO);
    t_var_gen.toc();
}


void class_measurement::compute_energy_variance_mpo(){
    t_var_mpo.tic();
    double L = superblock->Lblock2->size + superblock->Rblock2->size + 2.0;
    Eigen::Tensor<Scalar, 0> H2 =
            superblock->Lblock2->block
                    .contract(mps_util->theta,             idx({0}  ,{1}))
                    .contract(superblock->HA->MPO,                idx({1,3},{0,2}))
                    .contract(superblock->HB->MPO,                idx({4,2},{0,2}))
                    .contract(superblock->HA->MPO,                idx({1,3},{0,2}))
                    .contract(superblock->HB->MPO,                idx({4,3},{0,2}))
                    .contract(mps_util->theta.conjugate(), idx({0,3,5},{1,0,2}))
                    .contract(superblock->Rblock2->block,         idx({0,3,1,2},{0,1,2,3}));

    if (sim_type == SimulationType::iDMRG) {
        t_var_mpo.toc();
        variance_mpo_all_sites = std::real(H2(0))/2.0;
        variance_mpo           = std::real(H2(0))/2.0/L;
    }else{
        variance_mpo_all_sites = std::abs(H2(0) - energy_mpo_all_sites*energy_mpo_all_sites);
        variance_mpo           = variance_mpo_all_sites/L;
    }
    t_var_mpo.toc();
}


template<typename T>
double class_measurement::compute_energy_variance_mpo(const T * theta_ptr, Eigen::DSizes<long,4> dsizes, double energy_all_sites){
    Eigen::Tensor<std::complex<double>,4> theta = Eigen::TensorMap<const Eigen::Tensor<const T,4>> (theta_ptr, dsizes).template cast<std::complex<double>>();
//    t_var_mpo.tic();
    double L = superblock->Lblock2->size + superblock->Rblock2->size + 2.0;
    Eigen::Tensor<Scalar, 0> H2 =
            superblock->Lblock2->block
                    .contract(theta,                              idx({0}  ,{1}))
                    .contract(superblock->HA->MPO,                idx({1,3},{0,2}))
                    .contract(superblock->HB->MPO,                idx({4,2},{0,2}))
                    .contract(superblock->HA->MPO,                idx({1,3},{0,2}))
                    .contract(superblock->HB->MPO,                idx({4,3},{0,2}))
                    .contract(theta.conjugate(),                  idx({0,3,5},{1,0,2}))
                    .contract(superblock->Rblock2->block,         idx({0,3,1,2},{0,1,2,3}));
    if (sim_type == SimulationType::iDMRG) {
//        t_var_mpo.toc();
        return std::real(H2(0))/2.0/L;
    }else{
//        t_var_mpo.toc();
        return std::abs(H2(0) - energy_all_sites*energy_all_sites)/L;
    }
}

template double class_measurement::compute_energy_variance_mpo(const double *,Eigen::DSizes<long,4>,double);
template double class_measurement::compute_energy_variance_mpo(const std::complex<double> *, Eigen::DSizes<long,4>, double);


void class_measurement::compute_energy_variance_ham(){
    t_var_ham.tic();
    const Eigen::Tensor<Scalar,4> & theta_evn                  = mps_util->theta_evn_normalized;
    const Eigen::Tensor<Scalar,4> & theta_odd                  = mps_util->theta_odd_normalized;
    const Eigen::Tensor<Scalar,3> & LBGA                       = mps_util->LBGA;
    const Eigen::Tensor<Scalar,3> & LAGB                       = mps_util->LAGB;
    const Eigen::Tensor<Scalar,2> & l_evn                      = mps_util->l_evn;
    const Eigen::Tensor<Scalar,2> & r_evn                      = mps_util->r_evn;
    const Eigen::Tensor<Scalar,2> & l_odd                      = mps_util->l_odd;
    const Eigen::Tensor<Scalar,2> & r_odd                      = mps_util->r_odd;
    const Eigen::Tensor<Scalar,4> & transfer_matrix_evn        = mps_util->transfer_matrix_evn;
    const Eigen::Tensor<Scalar,4> & transfer_matrix_odd        = mps_util->transfer_matrix_odd;
    const Eigen::Tensor<Scalar,4> & transfer_matrix_LBGA       = mps_util->transfer_matrix_LBGA;
    const Eigen::Tensor<Scalar,4> & transfer_matrix_LAGB       = mps_util->transfer_matrix_LAGB;

    Eigen::Tensor<Scalar,4> h0 =  Matrix_to_Tensor((h_evn - E_evn(0)*MatrixType<Scalar>::Identity(4,4)).eval(), 2,2,2,2);
    Eigen::Tensor<Scalar,4> h1 =  Matrix_to_Tensor((h_odd - E_odd(0)*MatrixType<Scalar>::Identity(4,4)).eval(), 2,2,2,2);

    Eigen::Tensor<Scalar,0> E2AB =
            theta_evn
                    .contract(h0,                     idx({0, 2}, {0, 1}))
                    .contract(h0,                     idx({2, 3}, {0, 1}))
                    .contract(theta_evn.conjugate(),  idx({2, 3}, {0, 2}))
                    .contract(l_evn,                  idx({0, 2}, {0, 1}))
                    .contract(r_evn,                  idx({0, 1}, {0, 1}));


    Eigen::Tensor<Scalar, 0> E2BA =
            theta_odd
                    .contract(h1,                    idx({0, 2}, {0, 1}))
                    .contract(h1,                    idx({2, 3}, {0, 1}))
                    .contract(theta_odd.conjugate(), idx({2, 3}, {0, 2}))
                    .contract(l_odd,                 idx({0, 2}, {0, 1}))
                    .contract(r_odd,                 idx({0, 1}, {0, 1}));



    Eigen::Tensor<Scalar,5> thetaABA = theta_evn.contract(LBGA, idx({3},{1}));
    Eigen::Tensor<Scalar,5> thetaBAB = theta_odd.contract(LAGB, idx({3},{1}));

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
                    .contract(h0,                   idx({2,3},{0,1}))
                    .contract(h1,                   idx({0,3},{0,1}))
                    .contract(thetaBAB.conjugate(), idx({3,4,2},{0,2,3}))
                    .contract(l_odd,                idx({0,2},{0,1}))
                    .contract(r_evn,                idx({0,1},{0,1})) ;


    Eigen::Tensor<Scalar,2> E2d_L_evn =
            theta_evn
                    .contract(h0,                    idx({0, 2}, {0, 1}))
                    .contract(theta_evn.conjugate(), idx({2, 3}, {0, 2}))
                    .contract(l_evn,                 idx({0, 2}, {0, 1}));

    Eigen::Tensor<Scalar,2> E2d_R_evn =
            theta_evn
                    .contract(h0,                    idx({0, 2}, {0, 1}))
                    .contract(theta_evn.conjugate(), idx({2, 3}, {0, 2}))
                    .contract(r_evn,                 idx({1, 3}, {0, 1}));

    Eigen::Tensor<Scalar,2> E2d_L_odd  =
            theta_odd
                    .contract(h1,                     idx({0, 2}, {0, 1}))
                    .contract(theta_odd.conjugate(),  idx({2, 3}, {0, 2}))
                    .contract(l_odd,                  idx({0, 2}, {0, 1}));


    Eigen::Tensor<Scalar,2> E2d_R_odd =
            theta_odd
                    .contract(h1,                     idx({0, 2}, {0, 1}))
                    .contract(theta_odd.conjugate(),  idx({2, 3}, {0, 2}))
                    .contract(r_odd,                  idx({1, 3}, {0, 1}));

    Eigen::array<Eigen::IndexPair<long>,0> pair = {};
    Eigen::Tensor<Scalar,4> fixpoint_evn = r_evn.contract(l_evn, pair);
    Eigen::Tensor<Scalar,4> fixpoint_odd = r_odd.contract(l_odd, pair);

    long sizeLA = superblock->MPS->chiC();
    long sizeLB = superblock->MPS->chiB();
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

    variance_ham = std::real(0.5*(e2ab + e2ba) + 0.5*(e2aba_1  + e2bab_1  + e2aba_2  + e2bab_2 )  + e2lrpabab + e2lrpabba + e2lrpbaba  + e2lrpbaab) ;
    t_var_ham.toc();
}




void class_measurement::compute_finite_chain_energy_variance(){
    Eigen::Tensor<Scalar,4> L = env_storage->get_ENV2_L().front().block;
    Eigen::Tensor<Scalar,4> R = env_storage->get_ENV2_R().back().block;
    Eigen::Tensor<Scalar,4> temp;

    auto mpsL  = env_storage->get_MPS_L().begin();
    auto mpoL  = env_storage->get_MPO_L().begin();
    auto endL  = env_storage->get_MPS_L().end  ();
    while(mpsL != endL){
        const Eigen::Tensor<Scalar,1> &LA = mpsL->get_L(); //std::get<0>(*mpsL);
        const Eigen::Tensor<Scalar,3> &GA = mpsL->get_G(); //std::get<1>(*mpsL);
        assert(LA.dimension(0) == L.dimension(0));
        assert(LA.dimension(0) == GA.dimension(1));

        temp = L.contract(asDiagonal(LA), idx({0},{0}))
                .contract(asDiagonal(LA), idx({0},{0}))
                .contract(mpoL->get()->MPO,           idx({0},{0}))
                .contract(mpoL->get()->MPO,           idx({0,5},{0,2}))
                .contract(GA,                  idx({0,3},{1,0}))
                .contract(GA.conjugate(),      idx({0,3},{1,0}))
                .shuffle(array4{2,3,0,1});
        L = temp;
        mpsL++;
        mpoL++;
    }


    //Contract the center point
    auto &MPS_C = env_storage->get_MPS_C();
    temp = L.contract(asDiagonal(MPS_C) , idx({0},{0}))
            .contract(asDiagonal(MPS_C) , idx({0},{0}))
            .shuffle(array4{2,3,0,1});
    L = temp;

    //Contract the right half of the chain
    auto mpsR  = env_storage->get_MPS_R().begin();
    auto mpoR  = env_storage->get_MPO_R().begin();
    auto endR  = env_storage->get_MPS_R().end  ();
    while(mpsR != endR){
        const Eigen::Tensor<Scalar,3> &GB = mpsR->get_G(); // std::get<0>(*mpsR);
        const Eigen::Tensor<Scalar,1> &LB = mpsR->get_L(); // std::get<1>(*mpsR);
        assert(GB.dimension(1) == L.dimension(0));
        assert(LB.dimension(0) == GB.dimension(2));
        temp = L.contract(GB,            idx({0},{1}))
                .contract(mpoR->get()->MPO,     idx({1,3},{0,2}))
                .contract(mpoR->get()->MPO,     idx({1,4},{0,2}))
                .contract(GB.conjugate(),idx({0,4},{1,0}))
                .contract(asDiagonal(LB),idx({0},{0}))
                .contract(asDiagonal(LB),idx({2},{0}))
                .shuffle(array4{2,3,0,1});
        L = temp;
        mpsR++;
        mpoR++;
    }

    Eigen::Tensor<Scalar,0> H2_all_sites = L.contract(R, idx({0,1,2,3},{0,1,2,3}));
    variance_chain = std::abs(std::real(H2_all_sites(0)) - energy_chain*energy_chain) / env_storage->get_length() ;
}
