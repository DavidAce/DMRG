//
// Created by david on 2019-01-21.
//

#include "class_finite_chain_ops.h"
#include "class_finite_chain_state.h"
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mpo.h>
#include <general/nmspc_tensor_extra.h>
#include <general/class_svd_wrapper.h>
#include <general/nmspc_quantum_mechanics.h>

#include <iomanip>
using namespace Textra;

void class_finite_chain_ops::apply_mpo(state_ptr state,const Eigen::Tensor<Scalar,4> mpo,const Eigen::Tensor<Scalar,3> Ledge, const Eigen::Tensor<Scalar,3> Redge){
    // First the left side
    auto mps            = state->ref_MPS_L().begin();
    auto mpsR_it        = state->ref_MPS_R().rbegin();
    long mpoDim = mpo.dimension(0);
    auto I_tensor1 = Eigen::Tensor<Scalar,1>(mpoDim);
    I_tensor1.setConstant(1);


    // Apply MPO's on Gamma matrices and
    // increase the size on all Lambdas except at the edges

    for (auto & mps : state->ref_MPS_L){
        if (mps.get_position() == 0) continue;
        auto [d,chiL,chiR] = mps.get_dims();
        Eigen::Tensor<Scalar,3> G_temp = mps.get_G()
                .contract(mpo, idx({0},{2}))
                .shuffle(array5{4,2,0,3,1})
                .reshape(array3{d, chiL*mpoDim, chiR*mpoDim});
        Eigen::Tensor<Scalar,1> L_temp = mps.get_L().contract(I_tensor1, idx()).reshape(array1{chiL * mpoDim});
        mps.set_G(G_temp);
        mps.set_L(L_temp);
//        std::cout << "Current Lambda: \n" << std::setprecision(16)  << mps.get_L() << std::endl;

    }

    for (auto & mps : state->ref_MPS_R){
        if (mps.get_position() == state->max_length - 1) continue;
        auto [d,chiL,chiR] = mps.get_dims();
        Eigen::Tensor<Scalar,3> G_temp = mps.get_G()
                .contract(mpo, idx({0},{2}))
                .shuffle(array5{4,2,0,3,1})
                .reshape(array3{d, chiL*mpoDim, chiR*mpoDim});
        Eigen::Tensor<Scalar,1> L_temp = mps.get_L().contract(I_tensor1, idx()).reshape(array1{chiR * mpoDim});
        mps.set_G(G_temp);
        mps.set_L(L_temp);
//        std::cout << "Current Lambda: \n" << std::setprecision(16) << mps.get_L() << std::endl;

    }

    // Take care of the edges. It's enough to apply the left and right MPO-edges on Gamma's and
    // leave edge lambda's unchanged.
    {
        auto[d, chiL, chiR] = state->ref_MPS_L().front().get_dims();
        Eigen::Tensor<Scalar, 3> G_temp =
                state->ref_MPS_L().front().get_G()
                        .contract(mpo, idx({0}, {2}))
                        .shuffle(array5{4, 2, 0, 3, 1})
                        .contract(Ledge, idx({1, 2}, {2, 0}))
                        .reshape(array3{d, mpoDim*chiR, Ledge.dimension(1)})
                        .shuffle(array3{0, 2, 1});
        state->ref_MPS_L().front().set_G(G_temp);

    }
    {
        auto[d, chiL, chiR] = state->ref_MPS_R().back().get_dims();
        Eigen::Tensor<Scalar, 3> G_temp =
                state->ref_MPS_R().back().get_G()
                        .contract(mpo, idx({0}, {2}))
                        .shuffle(array5{4, 2, 0, 3, 1})
                        .contract(Redge, idx({3, 4}, {2, 0}))
                        .reshape(array3{d, mpoDim*chiL, Redge.dimension(1)});
        state->ref_MPS_R().back().set_G(G_temp);
    }


    //Take care of the middle
    auto chiC = state->ref_MPS_C().dimension(0);
    Eigen::Tensor<Scalar,1> C_temp = state->ref_MPS_C().contract(I_tensor1, idx()).reshape(array1{chiC * mpoDim});
    state->ref_MPS_C() = C_temp;


}


void class_finite_chain_ops::normalize_chain(state_ptr state){

    // First the left side
//    std::cout << "Starting Normalization" << std::endl;
    class_SVD<Scalar> svd;
    Eigen::Tensor<Scalar,2> C_temp;
    auto mpsL_it            = state->ref_MPS_L().begin();
    Eigen::Tensor<Scalar,3> A_temp = mpsL_it->get_A();

    while(mpsL_it != state->ref_MPS_L().end()) {
//        auto [d,chiL,chiR] = mpsL_it->get_dims();
        auto d    = A_temp.dimension(0);
        auto chiL = A_temp.dimension(1);
        auto chiR = A_temp.dimension(2);

        auto [U,S,V] = svd.decompose(A_temp, d*chiL, chiR);
        long chiC = S.dimension(0);

//        std::cout << "Current Lambda: \n" << mpsL_it->get_L() << std::endl;
//        std::cout << "Current S: \n" << S << std::endl;
//        std::cout << "Current dims L: [ " << d << " " << chiL << " " << chiR << " ]"<< std::endl;
//        std::cout << "  L dims: " << mpsL_it->get_L().dimensions() << std::endl;
//        std::cout << "  U dims: " << U.dimensions() << std::endl;
//        std::cout << "  S dims: " << S.dimensions() << std::endl;
//        std::cout << "  V dims: " << V.dimensions() << std::endl;

        Eigen::Tensor<Scalar,3> G_temp =
                Textra::asDiagonalInversed(mpsL_it->get_L())
                .contract(U.reshape(Textra::array3{d,chiL,chiC}), idx({1},{1}))
                .shuffle(array3{1,0,2});
        mpsL_it->set_G(G_temp);

        Eigen::Tensor<Scalar,2> leftcanonical =
                mpsL_it->get_A()
                .contract(mpsL_it->get_A().conjugate(), idx({0,1},{0,1}));
//        std::cout << "Left Canonical: \n" << std::setprecision(16) << leftcanonical.real() << std::endl;


        mpsL_it++;
        if (mpsL_it == state->ref_MPS_L().end()){
            C_temp = Textra::asDiagonal(S)
                    .contract(V               , idx({1},{0}))
                    .contract(Textra::asDiagonal(state->ref_MPS_C()), idx({1},{0}));
            break;
        }
        A_temp =
                Textra::asDiagonal(S)
                .contract(V               , idx({1},{0}))
                .contract(mpsL_it->get_A(), idx({1},{1}))
                .shuffle(array3{1,0,2});

        mpsL_it->set_L(S);
    }


    auto mpsR_it            = state->ref_MPS_R().rbegin();
    Eigen::Tensor<Scalar,3> B_temp = mpsR_it->get_B().shuffle(array3{1,0,2});
    while(mpsR_it != state->ref_MPS_R().rend()) {
        auto d    = B_temp.dimension(1);
        auto chiL = B_temp.dimension(0);
        auto chiR = B_temp.dimension(2);

        auto [U,S,V] = svd.decompose(B_temp, chiL, d*chiR);
        long chiC = S.dimension(0);

//        std::cout << "Current Lambda: \n" << mpsR_it->get_L() << std::endl;
//        std::cout << "Current S: \n" << S << std::endl;
//        std::cout << "Current dims L: [ " << d << " " << chiL << " " << chiR << " ]"<< std::endl;
//        std::cout << "  L dims: " << mpsR_it->get_L().dimensions() << std::endl;
//        std::cout << "  U dims: " << U.dimensions() << std::endl;
//        std::cout << "  S dims: " << S.dimensions() << std::endl;
//        std::cout << "  V dims: " << V.dimensions() << std::endl;

        Eigen::Tensor<Scalar,3> G_temp =
                V
                .reshape(Textra::array3{chiC,d,chiR})
                .shuffle(Textra::array3{1,0,2})
                .contract(Textra::asDiagonalInversed(mpsR_it->get_L()), idx({2},{0}));

        mpsR_it->set_G(G_temp);

        Eigen::Tensor<Scalar,2> rightcanonical =
                mpsR_it->get_B()
                        .contract(mpsR_it->get_B().conjugate(), idx({0,2},{0,2}));
//        std::cout << "Right Canonical: \n" << std::setprecision(16) << rightcanonical << std::endl;


        mpsR_it++;
        if (mpsR_it == state->ref_MPS_R().rend()){
            Eigen::Tensor<Scalar,2> C_temptemp =
                    C_temp
                    .contract(U, idx({1},{0}))
                    .contract(Textra::asDiagonal(S),idx({1},{0}));
            C_temp = C_temptemp;
            break;
        }
        B_temp =
                mpsR_it->get_B()
                .contract(U, idx({2},{0}))
                .contract(Textra::asDiagonal(S),idx({2},{0}))
                .shuffle(array3{1,0,2});
        mpsR_it->set_L(S);
    }


    // Now take care of the middle

    Eigen::Tensor<Scalar,4> ACB =
            state->ref_MPS_L().back().get_A()
            .contract(C_temp , idx({2},{0}))
            .contract(state->ref_MPS_R().front().get_B(), idx({2},{1}));
    auto [U,S,V] = svd.schmidt(ACB);
    state->ref_MPS_L().back().set_G (Textra::asDiagonalInversed(state->ref_MPS_L().back().get_L()).contract(U, idx({1},{1})).shuffle(array3{1,0,2}));
    state->ref_MPS_R().front().set_G(V.contract(Textra::asDiagonalInversed(state->ref_MPS_R().front().get_L()), idx({2},{0})));
    state->ref_MPS_C() = S;


}


void class_finite_chain_ops::set_parity_projected_mps(state_ptr state,const Eigen::MatrixXcd  paulimatrix) {
//    std::cout << "Dimensions before projecting \n";
//    for (auto & mps : state->ref_MPS_L) std::cout << "state->ref_MPS_L dims: " << mps.get_G().dimensions()<< std::endl;
//    std::cout << "state->ref_MPS_C dims: " << state->ref_MPS_C.dimensions() << std::endl;
//    for (auto & mps : state->ref_MPS_R) std::cout << "state->ref_MPS_R dims: " << mps.get_G().dimensions()<< std::endl;



    const auto [mpo,L,R]    = class_mpo::parity_selector_mpo(paulimatrix);
//    const auto [mpo,L,R]    = class_mpo::pauli_mpo(paulimatrix);

    apply_mpo(mpo,L,R);

    normalize_chain();

    refresh_environments();
    refresh_superblock();


}


void class_finite_chain_ops::refresh_superblock(state_ptr state, superblock_ptr superblock) {
    superblock->set_superblock(
            state->get_ENV2_L.back().block,
            state->get_ENV_L.back().block,
            state->get_MPO_L().back(),
            state->get_MPS_L().back()->get_L(),
            state->get_MPS_L().back()->get_G(),
            state->get_MPS_C(),
            state->get_MPS_R().front()->get_G(),
            state->get_MPS_R().front()->get_L(),
            state->get_MPO_R().front(),
            state->get_ENV_R.front().block,
            state->get_ENV2_R.front().block
            );
}


void class_finite_chain_ops::refresh_environments(state_ptr state){

    // Generate new environments

    auto mpsL_it   = state->ref_MPS_L().begin();
    auto mpoL_it   = state->ref_MPO_L().begin();
    auto envL_it   = state->ref_ENV_L().begin();
    auto env2L_it  = state->ref_ENV2_L().begin();
    int i = 0;
    while(mpsL_it != state->ref_MPS_L().end() and mpoL_it != state->ref_MPO_L().end()) {
        const Eigen::Tensor<Scalar,3> A = mpsL_it->get_A();

        Eigen::Tensor<Scalar,3> temp =
                envL_it->block
                .contract(A                   , idx({0},{1}))
                .contract(mpoL_it->get()->MPO , idx({1,2},{0,2}))
                .contract(A.conjugate()       , idx({0,3},{1,0}))
                .shuffle(array3{0,2,1});


        Eigen::Tensor<Scalar,4> temp2 =
                env2L_it->block
                        .contract(A                   , idx({0},{1}))
                        .contract(mpoL_it->get()->MPO , idx({1,3},{0,2}))
                        .contract(mpoL_it->get()->MPO , idx({1,4},{0,2}))
                        .contract(A.conjugate()       , idx({0,4},{1,0}))
                        .shuffle(array4{0,3,1,2});
        envL_it++;
        env2L_it++;
        mpsL_it++;
        mpoL_it++;
        if (mpsL_it == state->ref_MPS_L().end()) break;
        envL_it->block = temp;
        env2L_it->block = temp2;
        i++;
    }

    auto mpsR_it   = state->ref_MPS_R().rbegin();
    auto mpoR_it   = state->ref_MPO_R().rbegin();
    auto envR_it   = state->ref_ENV_R().rbegin();
    auto env2R_it  = state->ref_ENV2_R().rbegin();
    while(mpsR_it != state->ref_MPS_R().rend() and mpoR_it != state->ref_MPO_R().rend()){
        const Eigen::Tensor<Scalar,3> B = mpsR_it->get_B();
        Eigen::Tensor<Scalar,3> temp =
                envR_it->block
                .contract(B                   , idx({0},{2}))
                .contract(mpoR_it->get()->MPO , idx({1,2},{1,2}))
                .contract(B.conjugate()       , idx({0,3},{2,0}))
                .shuffle(array3{0,2,1});

        Eigen::Tensor<Scalar,4> temp2 =
                env2R_it->block
                .contract(B                   , idx({0},{2}))
                .contract(mpoR_it->get()->MPO , idx({1,3},{1,2}))
                .contract(mpoR_it->get()->MPO , idx({1,4},{1,2}))
                .contract(B.conjugate()       , idx({0,4},{2,0}))
                .shuffle(array4{0,3,1,2});


        envR_it++;
        env2R_it++;
        mpsR_it++;
        mpoR_it++;
        if (mpsR_it == state->ref_MPS_R().rend()) break;
        envR_it->block = temp;
        env2R_it->block = temp2;
    }

}
