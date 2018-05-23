//
// Created by david on 2018-02-21.
//
#include <iomanip>
#include <complex>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_mps.h>
#include <mps_routines/class_mpo.h>
#include <mps_routines/class_hamiltonian.h>
#include <general/nmspc_tensor_extra.h>
#include <general/class_svd_wrapper.h>
#include <general/class_arpack_eigsolver.h>
#include <mps_routines/class_finite_chain_storage.h>

using namespace std;
using namespace Textra;
using Scalar = class_measurement::Scalar;
using namespace std::complex_literals;


Scalar class_measurement::moment_generating_function(std::shared_ptr<class_mps> MPS_original,
                                                     std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec){
    t_temp1.tic();
    std::shared_ptr<class_mps> MPS_evolved = std::make_shared<class_mps>(*MPS_original);
    class_SVD<Scalar> SVD;
    SVD.setThreshold(1e-12);
    long chi_max = 4*MPS_evolved->LA.dimension(0);
    for (auto &Op: Op_vec) {
        //Evolve
        Textra::Tensor<Scalar, 4> theta_evo = Op.contract(MPS_evolved->get_theta(), idx({0, 1}, {0, 2})).shuffle(array4{0, 2, 1, 3});
        auto chiL = theta_evo.dimension(1);
        auto chiR = theta_evo.dimension(3);
        auto[U, S, V] = SVD.schmidt(theta_evo, 2, chiL,chi_max,chiR);
        MPS_evolved->LA= S;
        MPS_evolved->GA = asDiagonalInversed(MPS_evolved->LB).contract(U, idx({1}, {1})).shuffle(array3{1, 0, 2});
        MPS_evolved->GB = V.contract(asDiagonalInversed(MPS_evolved->LB), idx({2}, {0}));
        if (&Op != &Op_vec.back()) {
            MPS_evolved->swap_AB();
        }
    }
    long sizeLB = MPS_evolved->LB.size() * MPS_evolved->LB.size();
    //Normalize
    Tensor<Scalar,2> transfer_matrix_theta_evn       = MPS_evolved->get_transfer_matrix_theta_evn().reshape(array2{sizeLB,sizeLB});
    using namespace settings::precision;
    class_arpack_eigsolver<Scalar, Form::GENERAL> solver;
    solver.eig(transfer_matrix_theta_evn.data(), Ritz::LM, Side::R, (int)sizeLB, 1, eigMaxNcv, eigThreshold,eigMaxIter, false);
    auto new_theta_evn_normalized        = MPS_evolved->get_theta_evn(sqrt(solver.ref_eigvals()[0]));

    long sizeL = new_theta_evn_normalized.dimension(1) * MPS_original->theta_evn_normalized.dimension(1);
    long sizeR = new_theta_evn_normalized.dimension(3) * MPS_original->theta_evn_normalized.dimension(3);

    Tensor<Scalar,2> transfer_matrix_G =   new_theta_evn_normalized
            .contract(MPS_original->theta_evn_normalized.conjugate(), idx({0,2},{0,2}))
            .shuffle(array4{0,2,1,3})
            .reshape(array2{sizeL,sizeR});
    //Compute the characteristic function G(a).
    solver.eig(transfer_matrix_G.data(), Ritz::LM, Side::R, (int)transfer_matrix_G.dimension(0), 1, eigMaxNcv, eigThreshold,eigMaxIter, false);
    Scalar lambdaG = solver.ref_eigvals()[0];
    t_temp1.toc();
    return lambdaG;
}


std::pair<double,double> class_measurement::compute_infinite_moments_G(Scalar a, std::vector<Eigen::Tensor<Scalar, 4>> &Op_vec){
    t_var_gen.tic();
    using T = Scalar;
    //The following only works if superblock->MPS has been normalized! I.e, you have to have run MPS->compute_mps_components() prior.
    T lambdaG  = moment_generating_function(superblock->MPS, Op_vec);
//    LT lambdaG2 = (LT) moment_generating_function_2(superblock->MPS, Op_vec);
    T l        = superblock->H->mps_sites;
    T G        = pow(lambdaG,1.0/l);
    T logG     = log(lambdaG) * 1.0/l;
    T logGc    = log(conj(lambdaG) ) * 1.0/l;
    T O        = (logG - logGc)/(2.0*a);
    T VarO     = 2.0*log(abs(G))/ (a*a);
    t_var_gen.toc();
    return std::make_pair(real(O), real(VarO));
}


double class_measurement::compute_infinite_variance_MPO(){
    t_var_mpo.tic();
    if (sim == SimulationType::iDMRG) {
        Tensor<Scalar, 0> H2_minus_E2_two_sites =
                superblock->Lblock2->block
                .contract(superblock->MPS->theta,             idx({0}  ,{1}))
                .contract(superblock->HA->MPO,                idx({1,3},{0,2}))
                .contract(superblock->HB->MPO,                idx({4,2},{0,2}))
                .contract(superblock->HA->MPO,                idx({1,3},{0,2}))
                .contract(superblock->HB->MPO,                idx({4,3},{0,2}))
                .contract(superblock->MPS->theta.conjugate(), idx({0,3,5},{1,0,2}))
                .contract(superblock->Rblock2->block,         idx({0,3,1,2},{0,1,2,3}));
        double L =  superblock->Lblock2->size + superblock->Rblock2->size + 2;
        t_var_mpo.toc();
        return std::real( H2_minus_E2_two_sites(0)/ 2.0 /  L );
    }else{
        Tensor<Scalar, 0> H2_all_sites =
                superblock->Lblock2->block
                        .contract(superblock->MPS->theta,             idx({0}  ,{1}))
                        .contract(superblock->HA->MPO,                idx({1,3},{0,2}))
                        .contract(superblock->HB->MPO,                idx({4,2},{0,2}))
                        .contract(superblock->HA->MPO,                idx({1,3},{0,2}))
                        .contract(superblock->HB->MPO,                idx({4,3},{0,2}))
                        .contract(superblock->MPS->theta.conjugate(), idx({0,3,5},{1,0,2}))
                        .contract(superblock->Rblock2->block,         idx({0,3,1,2},{0,1,2,3}));

        double L = superblock->Lblock2->size + superblock->Rblock2->size + 2;
        t_var_mpo.toc();
        return std::real(H2_all_sites(0) - energy1*energy1*L*L ) / L;
    }
}

double class_measurement::compute_infinite_variance_H(){
    t_var_ham.tic();
    const Tensor<Scalar,4> & theta_evn                  = superblock->MPS->theta_evn_normalized;
    const Tensor<Scalar,4> & theta_odd                  = superblock->MPS->theta_odd_normalized;
    const Tensor<Scalar,3> & LBGA                       = superblock->MPS->LBGA;
    const Tensor<Scalar,3> & LAGB                       = superblock->MPS->LAGB;
    const Tensor<Scalar,2> & l_evn                      = superblock->MPS->l_evn;
    const Tensor<Scalar,2> & r_evn                      = superblock->MPS->r_evn;
    const Tensor<Scalar,2> & l_odd                      = superblock->MPS->l_odd;
    const Tensor<Scalar,2> & r_odd                      = superblock->MPS->r_odd;
    const Tensor<Scalar,4> & transfer_matrix_evn        = superblock->MPS->transfer_matrix_evn;
    const Tensor<Scalar,4> & transfer_matrix_odd        = superblock->MPS->transfer_matrix_odd;
    const Tensor<Scalar,4> & transfer_matrix_LBGA       = superblock->MPS->transfer_matrix_LBGA;
    const Tensor<Scalar,4> & transfer_matrix_LAGB       = superblock->MPS->transfer_matrix_LAGB;

    Tensor<Scalar,4> h0 =  Matrix_to_Tensor((superblock->H->h[0] - E_evn(0)*MatrixType<Scalar>::Identity(4,4)).eval(), 2,2,2,2);
    Tensor<Scalar,4> h1 =  Matrix_to_Tensor((superblock->H->h[1] - E_odd(0)*MatrixType<Scalar>::Identity(4,4)).eval(), 2,2,2,2);

    Tensor<Scalar,0> E2AB =
             theta_evn
            .contract(h0,                     idx({0, 2}, {0, 1}))
            .contract(h0,                     idx({2, 3}, {0, 1}))
            .contract(theta_evn.conjugate(),  idx({2, 3}, {0, 2}))
            .contract(l_evn,                  idx({0, 2}, {0, 1}))
            .contract(r_evn,                  idx({0, 1}, {0, 1}));


    Tensor<Scalar, 0> E2BA =
             theta_odd
            .contract(h1,                    idx({0, 2}, {0, 1}))
            .contract(h1,                    idx({2, 3}, {0, 1}))
            .contract(theta_odd.conjugate(), idx({2, 3}, {0, 2}))
            .contract(l_odd,                 idx({0, 2}, {0, 1}))
            .contract(r_odd,                 idx({0, 1}, {0, 1}));



    Tensor<Scalar,5> thetaABA = theta_evn.contract(LBGA, idx({3},{1}));
    Tensor<Scalar,5> thetaBAB = theta_odd.contract(LAGB, idx({3},{1}));

    Tensor<Scalar,0> E2ABA_1  =
                      thetaABA
                     .contract(h1,                   idx({2,3},{0,1}))
                     .contract(h0,                   idx({0,3},{0,1}))
                     .contract(thetaABA.conjugate(), idx({3,4,2},{0,2,3}))
                     .contract(l_evn,                idx({0,2},{0,1}))
                     .contract(r_odd,                idx({0,1},{0,1})) ;

    Tensor<Scalar,0> E2BAB_1  =
                    thetaBAB
                    .contract(h1,                   idx({0,2},{0,1}))
                    .contract(h0,                   idx({4,1},{0,1}))
                    .contract(thetaBAB.conjugate(), idx({2,3,4},{0,2,3}))
                    .contract(l_odd,                idx({0,2},{0,1}))
                    .contract(r_evn,                idx({0,1},{0,1})) ;

    Tensor<Scalar,0> E2ABA_2  =
                    thetaABA
                    .contract(h0,                   idx({0,2},{0,1}))
                    .contract(h1,                   idx({4,1},{0,1}))
                    .contract(thetaABA.conjugate(), idx({2,3,4},{0,2,3}))
                    .contract(l_evn,                idx({0,2},{0,1}))
                    .contract(r_odd,                idx({0,1},{0,1})) ;

    Tensor<Scalar,0> E2BAB_2  =
                     thetaBAB
                    .contract(h0,                   idx({2,3},{0,1}))
                    .contract(h1,                   idx({0,3},{0,1}))
                    .contract(thetaBAB.conjugate(), idx({3,4,2},{0,2,3}))
                    .contract(l_odd,                idx({0,2},{0,1}))
                    .contract(r_evn,                idx({0,1},{0,1})) ;


    Tensor<Scalar,2> E2d_L_evn =
            theta_evn
                    .contract(h0,                    idx({0, 2}, {0, 1}))
                    .contract(theta_evn.conjugate(), idx({2, 3}, {0, 2}))
                    .contract(l_evn,                 idx({0, 2}, {0, 1}));

    Tensor<Scalar,2> E2d_R_evn =
                     theta_evn
                    .contract(h0,                    idx({0, 2}, {0, 1}))
                    .contract(theta_evn.conjugate(), idx({2, 3}, {0, 2}))
                    .contract(r_evn,                 idx({1, 3}, {0, 1}));

    Tensor<Scalar,2> E2d_L_odd  =
                     theta_odd
                    .contract(h1,                     idx({0, 2}, {0, 1}))
                    .contract(theta_odd.conjugate(),  idx({2, 3}, {0, 2}))
                    .contract(l_odd,                  idx({0, 2}, {0, 1}));


    Tensor<Scalar,2> E2d_R_odd =
                     theta_odd
                    .contract(h1,                     idx({0, 2}, {0, 1}))
                    .contract(theta_odd.conjugate(),  idx({2, 3}, {0, 2}))
                    .contract(r_odd,                  idx({1, 3}, {0, 1}));

    Eigen::array<Eigen::IndexPair<long>,0> pair = {};
    Tensor<Scalar,4> fixpoint_evn = r_evn.contract(l_evn, pair);
    Tensor<Scalar,4> fixpoint_odd = r_odd.contract(l_odd, pair);

    long sizeLA = superblock->MPS->LA.size();
    long sizeLB = superblock->MPS->LB.size();
    Tensor<Scalar,2> one_minus_transfer_matrix_evn = Matrix_to_Tensor2(MatrixType<Scalar>::Identity(sizeLB*sizeLB, sizeLA*sizeLA).eval()) - (transfer_matrix_evn-fixpoint_evn).reshape(array2{sizeLB*sizeLB, sizeLA*sizeLA});
    Tensor<Scalar,2> one_minus_transfer_matrix_odd = Matrix_to_Tensor2(MatrixType<Scalar>::Identity(sizeLA*sizeLA, sizeLB*sizeLB).eval()) - (transfer_matrix_odd-fixpoint_odd).reshape(array2{sizeLA*sizeLA, sizeLB*sizeLB});

    class_SVD<Scalar> SVD;
    Tensor<Scalar,4> E_evn_pinv  = SVD.pseudo_inverse(one_minus_transfer_matrix_evn).reshape(array4{sizeLB,sizeLB,sizeLA,sizeLA});
    Tensor<Scalar,4> E_odd_pinv  = SVD.pseudo_inverse(one_minus_transfer_matrix_odd).reshape(array4{sizeLA,sizeLA,sizeLB,sizeLB});

    Tensor<Scalar,0> E2LRP_ABAB  = E2d_L_evn.contract(E_evn_pinv,idx({0,1},{0,1})).contract(E2d_R_evn,idx({0,1},{0,1}));
    Tensor<Scalar,0> E2LRP_ABBA  = E2d_L_evn.contract(transfer_matrix_LBGA, idx({0,1},{0,1})).contract(E_odd_pinv,idx({0,1},{0,1})).contract(E2d_R_odd,idx({0,1},{0,1}));
    Tensor<Scalar,0> E2LRP_BABA  = E2d_L_odd.contract(E_odd_pinv,idx({0,1},{0,1})).contract(E2d_R_odd,idx({0,1},{0,1}));
    Tensor<Scalar,0> E2LRP_BAAB  = E2d_L_odd.contract(transfer_matrix_LAGB, idx({0,1},{0,1})).contract(E_evn_pinv,idx({0,1},{0,1})).contract(E2d_R_evn,idx({0,1},{0,1}));


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

    Scalar VarE  = 0.5*(e2ab + e2ba) + 0.5*(e2aba_1  + e2bab_1  + e2aba_2  + e2bab_2 )  + e2lrpabab + e2lrpabba + e2lrpbaba  + e2lrpbaab ;
    t_var_ham.toc();
    return real(VarE) ;
}



double class_measurement::compute_finite_variance(){
    Eigen::Tensor<Scalar,4> L = std::as_const(env_storage->ENV2_L).front().block;
    Eigen::Tensor<Scalar,4> R = std::as_const(env_storage->ENV2_R).back().block;
    Eigen::Tensor<Scalar,4> temp;

    auto mpsL  = std::as_const(env_storage->MPS_L).begin();
    auto mpoL  = std::as_const(env_storage->MPO_L).begin();
    auto endL  = std::as_const(env_storage->MPS_L).end  ();
    while(mpsL != endL){
        const Eigen::Tensor<Scalar,1> &LB_left = std::get<0>(*mpsL);
        const Eigen::Tensor<Scalar,3> &GA      = std::get<1>(*mpsL);
        assert(LB_left.dimension(0) == L.dimension(0));
        assert(LB_left.dimension(0) == GA.dimension(1));

        temp = L.contract(asDiagonal(LB_left), idx({0},{0}))
                .contract(asDiagonal(LB_left), idx({0},{0}))
                .contract(mpoL->MPO,           idx({0},{0}))
                .contract(mpoL->MPO,           idx({0,5},{0,2}))
                .contract(GA,                  idx({0,3},{1,0}))
                .contract(GA.conjugate(),      idx({0,3},{1,0}))
                .shuffle(array4{2,3,0,1});
        L = temp;
        mpsL++;
        mpoL++;
    }


    //Contract the center point
    auto &LA = std::as_const(env_storage->LA);
    temp = L.contract(asDiagonal(LA) , idx({0},{0}))
            .contract(asDiagonal(LA) , idx({0},{0}))
            .shuffle(array4{2,3,0,1});
    L = temp;

    //Contract the right half of the chain
    auto mpsR  = std::as_const(env_storage->MPS_R).begin();
    auto mpoR  = std::as_const(env_storage->MPO_R).begin();
    auto endR  = std::as_const(env_storage->MPS_R).end  ();
    while(mpsR != endR){
        const Eigen::Tensor<Scalar,3> &GB      = std::get<0>(*mpsR);
        const Eigen::Tensor<Scalar,1> &LB      = std::get<1>(*mpsR);
        assert(GB.dimension(1) == L.dimension(0));
        assert(LB.dimension(0) == GB.dimension(2));
        temp = L.contract(GB,            idx({0},{1}))
                .contract(mpoR->MPO,     idx({1,3},{0,2}))
                .contract(mpoR->MPO,     idx({1,4},{0,2}))
                .contract(GB.conjugate(),idx({0,4},{1,0}))
                .contract(asDiagonal(LB),idx({0},{0}))
                .contract(asDiagonal(LB),idx({2},{0}))
                .shuffle(array4{2,3,0,1});
        L = temp;
        mpsR++;
        mpoR++;
    }

    Eigen::Tensor<Scalar,0> H2_all_sites = L.contract(R, idx({0,1,2,3},{0,1,2,3}));
    double Var = (std::real(H2_all_sites(0)) - energy_chain*energy_chain) / env_storage->current_length ;
    std::cout << setprecision(16) << "H2 all sites: " << H2_all_sites(0) << " Variance: " << Var <<  " log10: " << std::log10(Var) << std::endl;
    return Var;
}
