//
// Created by david on 2019-01-30.
//


#include <mps_routines/nmspc_mps_tools.h>
#include <mps_routines/class_finite_chain_state.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mpo.h>
#include <model/class_hamiltonian_base.h>
#include <general/nmspc_tensor_extra.h>
#include <general/class_svd_wrapper.h>
#include <general/nmspc_quantum_mechanics.h>
#include <iomanip>
#include <spdlog/spdlog.h>



using Scalar         = std::complex<double>;
using namespace Textra;

std::list<Eigen::Tensor<Scalar,4>> MPS_Tools::Finite::Ops::make_mpo_list (const std::list<std::shared_ptr<class_hamiltonian_base>> &mpos_L, const std::list<std::shared_ptr<class_hamiltonian_base>> &mpos_R){
    std::list<Eigen::Tensor<Scalar,4>> mpos;
    for(auto &mpo_L : mpos_L){
        mpos.push_back(mpo_L->MPO);
    }
    for(auto &mpo_R : mpos_R){
        mpos.push_back(mpo_R->MPO);
    }
    return mpos;
}


void MPS_Tools::Finite::Ops::apply_mpo(class_finite_chain_state &state,const Eigen::Tensor<Scalar,4> mpo,const Eigen::Tensor<Scalar,3> Ledge, const Eigen::Tensor<Scalar,3> Redge){
    std::list<Eigen::Tensor<Scalar,4>> mpos(state.get_length(), mpo);
    apply_mpos(state,mpos,Ledge,Redge);
}



void MPS_Tools::Finite::Ops::apply_mpo2(class_finite_chain_state &state,const Eigen::Tensor<Scalar,4> mpo,const Eigen::Tensor<Scalar,3> Ledge, const Eigen::Tensor<Scalar,3> Redge){
    // First the left side
    long mpoDim = mpo.dimension(0);
    auto I_tensor1 = Eigen::Tensor<Scalar,1>(mpoDim);
    I_tensor1.setConstant(1);

    // Apply MPO's on Gamma matrices and
    // increase the size on all Lambdas

    for (auto & mps : state.get_MPS_L()){
        auto [d,chiL,chiR] = mps.get_dims();
        Eigen::Tensor<Scalar,3> G_temp =
                mps.get_G()
                        .contract(mpo, idx({0},{2}))
                        .shuffle(array5{4,2,0,3,1})
                        .reshape(array3{d, chiL*mpoDim, chiR*mpoDim});
        Eigen::Tensor<Scalar,1> L_temp = I_tensor1.contract(mps.get_L(), idx()).reshape(array1{chiL * mpoDim});
        mps.set_G(G_temp);
        mps.set_L(L_temp);
//        std::cout << "LCurrent Lambda: \n" << std::setprecision(16)  << mps.get_L() << std::endl;

    }

    for (auto & mps : state.get_MPS_R()){
//        if (mps.get_position() == state.get_length() - 1) continue;
        auto [d,chiL,chiR] = mps.get_dims();
        Eigen::Tensor<Scalar,3> G_temp =
                mps.get_G()
                        .contract(mpo, idx({0},{2}))
                        .shuffle(array5{4,2,0,3,1})
                        .reshape(array3{d, chiL*mpoDim, chiR*mpoDim});
        Eigen::Tensor<Scalar,1> L_temp = I_tensor1.contract(mps.get_L(), idx()).reshape(array1{chiR * mpoDim});
        mps.set_G(G_temp);
        mps.set_L(L_temp);
//        std::cout << "RCurrent Lambda: \n" << std::setprecision(16) << mps.get_L() << std::endl;

    }



    //Take care of the middle
    auto chiC = state.get_MPS_C().dimension(0);
    Eigen::Tensor<Scalar,1> C_temp = I_tensor1.contract(state.get_MPS_C(), idx()).reshape(array1{chiC * mpoDim});
    state.get_MPS_C() = C_temp;
//    std::cout << "C_temp: \n" << C_temp << std::endl;

    // Take care of the edges. It's enough to apply the left and right MPO-edges on Gamma's and
    // leave edge lambda's unchanged.
    {
//        auto[d, chiL, chiR] = state.get_MPS_L().front().get_dims();
        auto Ldim = Ledge.dimension(0);
        Eigen::Tensor<Scalar, 3> G_temp =
                Ledge
                        .shuffle(Textra::array3{0,2,1})
                        .reshape(Textra::array2{Ldim*mpoDim,Ldim})
                        .contract(state.get_MPS_L().front().get_A(),Textra::idx({0},{1}))
                        .shuffle(Textra::array3{1,0,2});
        state.get_MPS_L().front().set_G(G_temp);
        state.get_MPS_L().front().set_L(Eigen::Tensor<Scalar,1>(Ldim).constant(1));
//        std::cout << "LEdge Lambda: \n" << std::setprecision(16)  << state.get_MPS_L().front().get_L() << std::endl;

    }
    {
//        auto[d, chiL, chiR] = state.get_MPS_R().back().get_dims();
        auto Rdim = Redge.dimension(0);
        Eigen::Tensor<Scalar, 3> G_temp =
                Redge
                        .shuffle(Textra::array3{0,2,1})
                        .reshape(Textra::array2{Rdim*mpoDim,Rdim})
                        .contract(state.get_MPS_R().back().get_B(),Textra::idx({0},{2}))
                        .shuffle(Textra::array3{1,2,0});
        state.get_MPS_R().back().set_G(G_temp);
        state.get_MPS_R().back().set_L(Eigen::Tensor<Scalar,1>(Rdim).constant(1));
//        std::cout << "REdge Lambda: \n" << std::setprecision(16)  << state.get_MPS_R().back().get_L() << std::endl;

    }


}

void MPS_Tools::Finite::Ops::apply_mpo3(class_finite_chain_state &state,const Eigen::Tensor<Scalar,4> mpo,const Eigen::Tensor<Scalar,3> Ledge, const Eigen::Tensor<Scalar,3> Redge){
    // First the left side
    long mpoDim = mpo.dimension(0);
    auto I_tensor1 = Eigen::Tensor<Scalar,1>(mpoDim);
    I_tensor1.setConstant(1);

    // Apply MPO's on Gamma matrices and
    // increase the size on all Lambdas

    for (auto & mps : state.get_MPS_L()){
        auto [d,chiL,chiR] = mps.get_dims();
        Eigen::Tensor<Scalar,3> G_temp =
                mps.get_A()
                        .contract(mpo, idx({0},{2}))
                        .shuffle(array5{4,2,0,3,1})
                        .reshape(array3{d, chiL*mpoDim, chiR*mpoDim});
//        Eigen::Tensor<Scalar,1> L_temp = mps.get_L().contract(I_tensor1, idx()).reshape(array1{chiL * mpoDim});
        double normalization = 1.0 ;// / std::sqrt(chiL*mpoDim);
        mps.set_L(Eigen::Tensor<Scalar,1>(chiL*mpoDim).constant(normalization));
//        std::cout << "LCurrent Lambda  : \n" << std::setprecision(16)  << mps.get_L() << std::endl;
//        std::cout << "LCurrent G_before: \n" << std::setprecision(16)  << mps.get_G() << std::endl;
//        std::cout << "LCurrent G_temp  : \n" << std::setprecision(16)  << G_temp << std::endl;
//        std::cout << "LCurrent mpo     : \n" << std::setprecision(16)  << mpo << std::endl;
        mps.set_G(G_temp);
//        std::cout << "LCurrent G_after: \n" << std::setprecision(16)  << mps.get_G() << std::endl;


    }

    for (auto & mps : state.get_MPS_R()){
//        if (mps.get_position() == state.get_length() - 1) continue;
        auto [d,chiL,chiR] = mps.get_dims();
        Eigen::Tensor<Scalar,3> G_temp =
                mps.get_B()
                        .contract(mpo, idx({0},{2}))
                        .shuffle(array5{4,2,0,3,1})
                        .reshape(array3{d, chiL*mpoDim, chiR*mpoDim});
//        Eigen::Tensor<Scalar,1> L_temp = mps.get_L().contract(I_tensor1, idx()).reshape(array1{chiR * mpoDim});
        mps.set_G(G_temp);
        double normalization = 1.0 ;/// std::sqrt(chiR*mpoDim);
        mps.set_L(Eigen::Tensor<Scalar,1>(chiR*mpoDim).constant(normalization));
//        std::cout << "RCurrent Lambda: \n" << std::setprecision(16) << mps.get_L() << std::endl;

    }



    //Take care of the middle
    auto chiC = state.get_MPS_C().dimension(0);
    Eigen::Tensor<Scalar,1> C_temp =  I_tensor1.contract(state.get_MPS_C(), idx()).reshape(array1{chiC * mpoDim});
//    Eigen::Tensor<Scalar,1> C_temp = state.get_MPS_C().contract(I_tensor1, idx()).reshape(array1{chiC * mpoDim});
    state.get_MPS_C() = C_temp;
//    double normalization = 1.0;// / std::sqrt(chiC*mpoDim);
//    state.get_MPS_C() = Eigen::Tensor<Scalar,1>(chiC*mpoDim).constant(normalization);

//    std::cout << "C_temp: \n" << C_temp << std::endl;

    // Take care of the edges. It's enough to apply the left and right MPO-edges on Gamma's and
    // leave edge lambda's unchanged.
    {
//        auto[d, chiL, chiR] = state.get_MPS_L().front().get_dims();
        auto Ldim = Ledge.dimension(0);
        Eigen::Tensor<Scalar, 3> G_temp =
                Ledge
                        .shuffle(Textra::array3{0,2,1})
                        .reshape(Textra::array2{Ldim*mpoDim,Ldim})
                        .contract(state.get_MPS_L().front().get_A(),Textra::idx({0},{1}))
                        .shuffle(Textra::array3{1,0,2});
        state.get_MPS_L().front().set_G(G_temp);
        double normalization = 1.0 ;// / std::sqrt(Ldim*mpoDim);
        state.get_MPS_L().front().set_L(Eigen::Tensor<Scalar,1>(Ldim).constant(normalization));
//        std::cout << "LEdge Lambda: \n" << std::setprecision(16)  << state.get_MPS_L().front().get_L() << std::endl;
    }
    {
//        auto[d, chiL, chiR] = state.get_MPS_R().back().get_dims();
        auto Rdim = Redge.dimension(0);
        Eigen::Tensor<Scalar, 3> G_temp =
                Redge
                        .shuffle(Textra::array3{0,2,1})
                        .reshape(Textra::array2{Rdim*mpoDim,Rdim})
                        .contract(state.get_MPS_R().back().get_B(),Textra::idx({0},{2}))
                        .shuffle(Textra::array3{1,2,0});
        state.get_MPS_R().back().set_G(G_temp);
        double normalization = 1.0 ;//  std::sqrt(Rdim*mpoDim);
        state.get_MPS_R().back().set_L(Eigen::Tensor<Scalar,1>(Rdim).constant(normalization));
//        std::cout << "REdge Lambda: \n" << std::setprecision(16)  << state.get_MPS_R().back().get_L() << std::endl;

    }


}


void MPS_Tools::Finite::Ops::apply_mpos(class_finite_chain_state &state,const std::list<Eigen::Tensor<Scalar,4>> &mpos,const Eigen::Tensor<Scalar,3> Ledge, const Eigen::Tensor<Scalar,3> Redge){
    assert(mpos.size() == state.get_length() and "ERROR: number of mpo's doesn't match the number of sites on the system");

    // Apply MPO's on Gamma matrices and
    // increase the size on all Lambdas
    auto mpo = mpos.begin();
    for (auto & mps : state.get_MPS_L()){
        long mpoDimL = mpo->dimension(0);
        long mpoDimR = mpo->dimension(1);
        auto [d,chiL,chiR] = mps.get_dims();
        Eigen::Tensor<Scalar,3> G_temp =
                mps.get_G()
                        .contract(*mpo, idx({0},{2}))
                        .shuffle(array5{4,2,0,3,1})
                        .reshape(array3{d, chiL*mpoDimL, chiR*mpoDimR});

        mps.set_G(G_temp);

        auto I_tensor1 = Eigen::Tensor<Scalar,1>(mpoDimL);
        I_tensor1.setConstant(1);
        Eigen::Tensor<Scalar,1> L_temp = I_tensor1.contract(mps.get_L(), idx()).reshape(array1{chiL * mpoDimL});
        mps.set_L(L_temp);
        mpo++;
    }


    //Take care of the middle
    {
        long mpoDim = mpo->dimension(0);
        auto I_tensor1 = Eigen::Tensor<Scalar, 1>(mpoDim);
        auto chiC = state.get_MPS_C().dimension(0);
        Eigen::Tensor<Scalar, 1> C_temp = I_tensor1.contract(state.get_MPS_C(), idx()).reshape(array1{chiC * mpoDim});
        state.get_MPS_C() = C_temp;
    }

    for (auto & mps : state.get_MPS_R()){
        long mpoDimL = mpo->dimension(0);
        long mpoDimR = mpo->dimension(1);
        auto [d,chiL,chiR] = mps.get_dims();
        Eigen::Tensor<Scalar,3> G_temp =
                mps.get_G()
                        .contract(*mpo, idx({0},{2}))
                        .shuffle(array5{4,2,0,3,1})
                        .reshape(array3{d, chiL*mpoDimL, chiR*mpoDimR});

        mps.set_G(G_temp);

        auto I_tensor1 = Eigen::Tensor<Scalar,1>(mpoDimR);
        Eigen::Tensor<Scalar,1> L_temp = I_tensor1.contract(mps.get_L(), idx()).reshape(array1{chiR * mpoDimR});
        mps.set_L(L_temp);

        mpo++;
    }


    // Take care of the edges. Apply the left and right MPO-edges on A's and B's
    // so the left and right most lambdas become ones
    {
        long mpoDimL = mpos.front().dimension(0);
//        auto[d, chiL, chiR] = state.get_MPS_L().front().get_dims();
        auto Ldim = Ledge.dimension(0);
        Eigen::Tensor<Scalar, 3> G_temp =
                Ledge
                        .shuffle(Textra::array3{0,2,1})
                        .reshape(Textra::array2{Ldim*mpoDimL,Ldim})
                        .contract(state.get_MPS_L().front().get_A(),Textra::idx({0},{1}))
                        .shuffle(Textra::array3{1,0,2});
        state.get_MPS_L().front().set_G(G_temp);
        state.get_MPS_L().front().set_L(Eigen::Tensor<Scalar,1>(Ldim).constant(1.0));
//        std::cout << "LEdge Lambda: \n" << std::setprecision(16)  << state.get_MPS_L().front().get_L() << std::endl;
    }
    {
        long mpoDimR = mpos.back().dimension(1);
//        auto[d, chiL, chiR] = state.get_MPS_R().back().get_dims();
        auto Rdim = Redge.dimension(0);
        Eigen::Tensor<Scalar, 3> G_temp =
                Redge
                        .shuffle(Textra::array3{0,2,1})
                        .reshape(Textra::array2{Rdim*mpoDimR,Rdim})
                        .contract(state.get_MPS_R().back().get_B(),Textra::idx({0},{2}))
                        .shuffle(Textra::array3{1,2,0});
        state.get_MPS_R().back().set_G(G_temp);
        state.get_MPS_R().back().set_L(Eigen::Tensor<Scalar,1>(Rdim).constant(1.0));
//        std::cout << "REdge Lambda: \n" << std::setprecision(16)  << state.get_MPS_R().back().get_L() << std::endl;

    }


}


void MPS_Tools::Finite::Ops::normalize_chain(class_finite_chain_state & state){

    // First the left side
//    std::cout << "Starting Normalization" << std::endl;
    class_SVD<Scalar> svd;
    Eigen::Tensor<Scalar,2> C_temp;
    auto mpsL_it            = state.get_MPS_L().begin();
    Eigen::Tensor<Scalar,3> A_temp = mpsL_it->get_A();

    while(mpsL_it != state.get_MPS_L().end()) {
//        auto [d,chiL,chiR] = mpsL_it->get_dims();
        auto d    = A_temp.dimension(0);
        auto chiL = A_temp.dimension(1);
        auto chiR = A_temp.dimension(2);
        std::setprecision(16);
        auto [U,S,V] = svd.decompose(A_temp, d*chiL, chiR);
//        std::cout << "A_temp: \n" << A_temp << std::endl;
//        std::cout << "U: \n" << U << std::endl;
//        std::cout << "S: \n" << S << std::endl;
//        std::cout << "V: \n" << V << std::endl;

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
        if (mpsL_it == state.get_MPS_L().end()){
            C_temp = Textra::asDiagonal(S)
                    .contract(V               , idx({1},{0}))
                    .contract(Textra::asDiagonal(state.get_MPS_C()), idx({1},{0}));
            break;
        }
//        std::cout << "A_next : \n" << mpsL_it->get_A() << std::endl;
//        std::cout << "G_next : \n" << mpsL_it->get_G() << std::endl;
//        std::cout << "L_next : \n" << mpsL_it->get_L() << std::endl;
        A_temp =
                Textra::asDiagonal(S)
                        .contract(V               , idx({1},{0}))
                        .contract(mpsL_it->get_A(), idx({1},{1}))
                        .shuffle(array3{1,0,2});

        mpsL_it->set_L(S);
    }


    auto mpsR_it            = state.get_MPS_R().rbegin();
    Eigen::Tensor<Scalar,3> B_temp = mpsR_it->get_B().shuffle(array3{1,0,2});
    while(mpsR_it != state.get_MPS_R().rend()) {
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
        if (mpsR_it == state.get_MPS_R().rend()){
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
            state.get_MPS_L().back().get_A()
                    .contract(C_temp , idx({2},{0}))
                    .contract(state.get_MPS_R().front().get_B(), idx({2},{1}));
    auto [U,S,V] = svd.schmidt(ACB);
    state.get_MPS_L().back().set_G (Textra::asDiagonalInversed(state.get_MPS_L().back().get_L()).contract(U, idx({1},{1})).shuffle(array3{1,0,2}));
    state.get_MPS_R().front().set_G(V.contract(Textra::asDiagonalInversed(state.get_MPS_R().front().get_L()), idx({2},{0})));
    state.get_MPS_C() = S;

}


void MPS_Tools::Finite::Ops::apply_energy_mpo_test(class_finite_chain_state &state, class_superblock &superblock) {

//    const auto [mpo,L,R]    = class_mpo::pauli_mpo(paulimatrix);
    auto hamiltonian_mpos = make_mpo_list(state.get_MPO_L(),state.get_MPO_R());

//    const auto mpo = superblock.HA->MPO;
    const auto L   = state.get_ENV_L().front().block;
    const auto R   = state.get_ENV_R().back().block;


    // Make a copy of the state on which to apply the mpo
    std::cout << std::setprecision(16);
    class_finite_chain_state state_copy1 = state;
    class_finite_chain_state state_copy2 = state;
    double overlap0 = MPS_Tools::Finite::Ops::overlap(state,state_copy1);
    double exp_val0 = MPS_Tools::Finite::Ops::expectation_value(state,state_copy1,hamiltonian_mpos,L,R);
    std::cout << "overlap 0: " << overlap0 << " |  per site: " << overlap0/state.get_length() << std::endl;
    std::cout << "energy  0: " << exp_val0 << " |  per site: " << exp_val0/state.get_length() << std::endl;

    apply_mpos(state_copy1,hamiltonian_mpos,L,R);
    apply_mpos(state_copy2,hamiltonian_mpos,L,R);

    double overlap11 = MPS_Tools::Finite::Ops::overlap(state,state_copy1);
    double exp_val11 = MPS_Tools::Finite::Ops::expectation_value(state,state_copy1,hamiltonian_mpos,L,R);
    std::cout << "overlap 11: " << overlap11 << " |  per site: " << overlap11/state.get_length() << std::endl;
    std::cout << "energy  11: " << exp_val11 << " |  per site: " << exp_val11/state.get_length() << std::endl;

    double overlap12 = MPS_Tools::Finite::Ops::overlap(state,state_copy1);
    double exp_val12 = MPS_Tools::Finite::Ops::expectation_value(state,state_copy1,hamiltonian_mpos,L,R);
    std::cout << "overlap 12: " << overlap12 << " |  per site: " << overlap12/state.get_length() << std::endl;
    std::cout << "energy  12: " << exp_val12 << " |  per site: " << exp_val12/state.get_length() << std::endl;



    normalize_chain(state_copy1);

    double overlap21 = MPS_Tools::Finite::Ops::overlap(state,state_copy1);
    double exp_val21 = MPS_Tools::Finite::Ops::expectation_value(state,state_copy1,hamiltonian_mpos,L,R);
    std::cout << "overlap per site 21: " << overlap21 << " |  per site: " << overlap21/state.get_length() << std::endl;
    std::cout << "energy  per site 21: " << exp_val21 << " |  per site: " << exp_val21/state.get_length() << std::endl;


    normalize_chain(state_copy2);

    double overlap22 = MPS_Tools::Finite::Ops::overlap(state,state_copy2);
    double exp_val22 = MPS_Tools::Finite::Ops::expectation_value(state,state_copy2,hamiltonian_mpos,L,R);
    std::cout << "overlap per site 22: " << overlap22 << " |  per site: " << overlap22/state.get_length() << std::endl;
    std::cout << "energy  per site 22: " << exp_val22 << " |  per site: " << exp_val22/state.get_length() << std::endl;
//

    class_finite_chain_state state_copy3 = state;
    normalize_chain(state_copy3);
    double overlap3 = MPS_Tools::Finite::Ops::overlap(state,state_copy3);
    std::cout << "overlap per site 3: " << overlap3 << " |  per site: " << overlap3/state.get_length() << std::endl;
//
//    rebuild_environments(state);
//
//    rebuild_superblock(state, superblock);

}

void MPS_Tools::Finite::Ops::set_parity_projected_mps(class_finite_chain_state &state, class_superblock &superblock,const Eigen::MatrixXcd  paulimatrix) {



}

void MPS_Tools::Finite::Ops::check_parity_properties(class_finite_chain_state &state) {
//    std::cout << "Dimensions before projecting \n";
//    for (auto & mps : state.ref_MPS_L) std::cout << "state.get_MPS_L dims: " << mps.get_G().dimensions()<< std::endl;
//    std::cout << "state.ref_MPS_C dims: " << state.get_MPS_C.dimensions() << std::endl;
//    for (auto & mps : state.ref_MPS_R) std::cout << "state.get_MPS_R dims: " << mps.get_G().dimensions()<< std::endl;



//    const auto [mpo,L,R]                   = class_mpo::parity_selector_mpo(paulimatrix);
//    const auto [mpo,L,R]                   = class_mpo::parity_projector_mpos(paulimatrix,state.get_length(), 1);
    spdlog::trace("Computing spin components");
    const auto sx = MPS_Tools::Finite::Measure::spin_component(state,qm::spinOneHalf::sx);
    const auto sy = MPS_Tools::Finite::Measure::spin_component(state,qm::spinOneHalf::sy);
    const auto sz = MPS_Tools::Finite::Measure::spin_component(state,qm::spinOneHalf::sz);
    std::cout << std::setprecision(16);
    std::cout << "<psi | sx | psi> = " << sx << std::endl;
    std::cout << "<psi | sy | psi> = " << sy << std::endl;
    std::cout << "<psi | sz | psi> = " << sz << std::endl;


    spdlog::trace("Computing parity projection operators");
    const auto [mpo_up_x,L_up_x,R_up_x]    = class_mpo::parity_projector_mpos(qm::spinOneHalf::sx,state.get_length(), 1);
    const auto [mpo_dn_x,L_dn_x,R_dn_x]    = class_mpo::parity_projector_mpos(qm::spinOneHalf::sx,state.get_length(), -1);
    const auto [mpo_up_y,L_up_y,R_up_y]    = class_mpo::parity_projector_mpos(qm::spinOneHalf::sy,state.get_length(), 1);
    const auto [mpo_dn_y,L_dn_y,R_dn_y]    = class_mpo::parity_projector_mpos(qm::spinOneHalf::sy,state.get_length(), -1);
    const auto [mpo_up_z,L_up_z,R_up_z]    = class_mpo::parity_projector_mpos(qm::spinOneHalf::sz,state.get_length(), 1);
    const auto [mpo_dn_z,L_dn_z,R_dn_z]    = class_mpo::parity_projector_mpos(qm::spinOneHalf::sz,state.get_length(), -1);

//    const auto [mpo,L,R]    = class_mpo::pauli_mpo(paulimatrix);


    // Make a copy of the state on which to apply the mpo
//    std::cout << std::setprecision(16);
//    for(const auto & mps: state.get_MPS_L()){
//        std::cout << "G state: \n" << mps.get_G() << std::endl;
//    }
    spdlog::trace("Copying states into temporaries");
    class_finite_chain_state state_up_x = state;
    class_finite_chain_state state_dn_x = state;
    class_finite_chain_state state_up_y = state;
    class_finite_chain_state state_dn_y = state;
    class_finite_chain_state state_up_z = state;
    class_finite_chain_state state_dn_z = state;
//    class_finite_chain_state state_copy = state;
    class_finite_chain_state state_copy2 = state;


//    for(const auto & mps: state_copy.get_MPS_L()){
//        std::cout << "G state_copy: \n" << mps.get_G() << std::endl;
//    }
//    apply_mpos(state_copy,mpo,L,R);
    spdlog::trace("Applying parity projection operatorss");
    apply_mpos(state_up_x,mpo_up_x,L_up_x,R_up_x);
    apply_mpos(state_dn_x,mpo_dn_x,L_dn_x,R_dn_x);
    apply_mpos(state_up_y,mpo_up_y,L_up_y,R_up_y);
    apply_mpos(state_dn_y,mpo_dn_y,L_dn_y,R_dn_y);
    apply_mpos(state_up_z,mpo_up_z,L_up_z,R_up_z);
    apply_mpos(state_dn_z,mpo_dn_z,L_dn_z,R_dn_z);
    spdlog::trace("Normalizing chains");

//    normalize_chain(state_copy);
    normalize_chain(state_up_x);
    normalize_chain(state_dn_x);
    normalize_chain(state_up_y);
    normalize_chain(state_dn_y);
    normalize_chain(state_up_z);
    normalize_chain(state_dn_z);

    spdlog::trace("Computing more spin components");
    std::cout << "<psi_up_x | sx | psi_up_x> = " << MPS_Tools::Finite::Measure::spin_component(state_up_x,qm::spinOneHalf::sx) << std::endl;
    std::cout << "<psi_up_x | sy | psi_up_x> = " << MPS_Tools::Finite::Measure::spin_component(state_up_x,qm::spinOneHalf::sy) << std::endl;
    std::cout << "<psi_up_x | sz | psi_up_x> = " << MPS_Tools::Finite::Measure::spin_component(state_up_x,qm::spinOneHalf::sz) << std::endl;
    std::cout << "<psi_dn_x | sx | psi_dn_x> = " << MPS_Tools::Finite::Measure::spin_component(state_dn_x,qm::spinOneHalf::sx) << std::endl;
    std::cout << "<psi_dn_x | sy | psi_dn_x> = " << MPS_Tools::Finite::Measure::spin_component(state_dn_x,qm::spinOneHalf::sy) << std::endl;
    std::cout << "<psi_dn_x | sz | psi_dn_x> = " << MPS_Tools::Finite::Measure::spin_component(state_dn_x,qm::spinOneHalf::sz) << std::endl;

    std::cout << "<psi_up_y | sx | psi_up_y> = " << MPS_Tools::Finite::Measure::spin_component(state_up_y,qm::spinOneHalf::sx) << std::endl;
    std::cout << "<psi_up_y | sy | psi_up_y> = " << MPS_Tools::Finite::Measure::spin_component(state_up_y,qm::spinOneHalf::sy) << std::endl;
    std::cout << "<psi_up_y | sz | psi_up_y> = " << MPS_Tools::Finite::Measure::spin_component(state_up_y,qm::spinOneHalf::sz) << std::endl;
    std::cout << "<psi_dn_y | sx | psi_dn_y> = " << MPS_Tools::Finite::Measure::spin_component(state_dn_y,qm::spinOneHalf::sx) << std::endl;
    std::cout << "<psi_dn_y | sy | psi_dn_y> = " << MPS_Tools::Finite::Measure::spin_component(state_dn_y,qm::spinOneHalf::sy) << std::endl;
    std::cout << "<psi_dn_y | sz | psi_dn_y> = " << MPS_Tools::Finite::Measure::spin_component(state_dn_y,qm::spinOneHalf::sz) << std::endl;

    std::cout << "<psi_up_z | sx | psi_up_z> = " << MPS_Tools::Finite::Measure::spin_component(state_up_z,qm::spinOneHalf::sx) << std::endl;
    std::cout << "<psi_up_z | sy | psi_up_z> = " << MPS_Tools::Finite::Measure::spin_component(state_up_z,qm::spinOneHalf::sy) << std::endl;
    std::cout << "<psi_up_z | sz | psi_up_z> = " << MPS_Tools::Finite::Measure::spin_component(state_up_z,qm::spinOneHalf::sz) << std::endl;
    std::cout << "<psi_dn_z | sx | psi_dn_z> = " << MPS_Tools::Finite::Measure::spin_component(state_dn_z,qm::spinOneHalf::sx) << std::endl;
    std::cout << "<psi_dn_z | sy | psi_dn_z> = " << MPS_Tools::Finite::Measure::spin_component(state_dn_z,qm::spinOneHalf::sy) << std::endl;
    std::cout << "<psi_dn_z | sz | psi_dn_z> = " << MPS_Tools::Finite::Measure::spin_component(state_dn_z,qm::spinOneHalf::sz) << std::endl;


    std::cout << std::setprecision(16);
    std::cout << "\nNormalization check" << std::endl;
    spdlog::trace("Computing overlaps");
    std::cout << "<psi_up_x|psi_up_x> = "  <<  MPS_Tools::Finite::Ops::overlap(state_up_x,state_up_x) <<std::endl;
    std::cout << "<psi_dn_x|psi_dn_x> = "  <<  MPS_Tools::Finite::Ops::overlap(state_dn_x,state_dn_x) <<std::endl;
    std::cout << "<psi_up_y|psi_up_y> = "  <<  MPS_Tools::Finite::Ops::overlap(state_up_y,state_up_y) <<std::endl;
    std::cout << "<psi_dn_y|psi_dn_y> = "  <<  MPS_Tools::Finite::Ops::overlap(state_dn_y,state_dn_y) <<std::endl;
    std::cout << "<psi_up_z|psi_up_z> = "  <<  MPS_Tools::Finite::Ops::overlap(state_up_z,state_up_z) <<std::endl;
    std::cout << "<psi_dn_z|psi_dn_z> = "  <<  MPS_Tools::Finite::Ops::overlap(state_dn_z,state_dn_z) <<std::endl;

    std::cout << "\nOverlaps with original state" << std::endl;
    std::cout << "<psi|psi_up_x> = "  <<  MPS_Tools::Finite::Ops::overlap(state,state_up_x) <<std::endl;
    std::cout << "<psi|psi_dn_x> = "  <<  MPS_Tools::Finite::Ops::overlap(state,state_dn_x) <<std::endl;
    std::cout << "<psi|psi_up_y> = "  <<  MPS_Tools::Finite::Ops::overlap(state,state_up_y) <<std::endl;
    std::cout << "<psi|psi_dn_y> = "  <<  MPS_Tools::Finite::Ops::overlap(state,state_dn_y) <<std::endl;
    std::cout << "<psi|psi_up_z> = "  <<  MPS_Tools::Finite::Ops::overlap(state,state_up_z) <<std::endl;
    std::cout << "<psi|psi_dn_z> = "  <<  MPS_Tools::Finite::Ops::overlap(state,state_dn_z) <<std::endl;

    std::cout << "\nOverlaps between up/down sectors" << std::endl;
    std::cout << "<psi_dn_x|psi_up_x> = "  <<  MPS_Tools::Finite::Ops::overlap(state_dn_x,state_up_x) <<std::endl;
    std::cout << "<psi_dn_y|psi_up_y> = "  <<  MPS_Tools::Finite::Ops::overlap(state_dn_y,state_up_y) <<std::endl;
    std::cout << "<psi_dn_z|psi_up_z> = "  <<  MPS_Tools::Finite::Ops::overlap(state_dn_z,state_up_z) <<std::endl;
    std::cout << "\nOverlaps between different direction sectors" << std::endl;
    std::cout << "<psi_up_x|psi_up_y> = "  <<  MPS_Tools::Finite::Ops::overlap(state_up_x,state_up_y) <<std::endl;
    std::cout << "<psi_up_x|psi_up_z> = "  <<  MPS_Tools::Finite::Ops::overlap(state_up_x,state_up_z) <<std::endl;
    spdlog::trace("Generating single hamiltonian MPO list");

    auto hamiltonian_mpos = make_mpo_list(state.get_MPO_L(),state.get_MPO_R());
    auto & Ledge = state.get_ENV_L().front().block;
    auto & Redge = state.get_ENV_R().back().block ;
    auto & Ledge2 = state.get_ENV2_L().front().block;
    auto & Redge2 = state.get_ENV2_R().back().block ;


    apply_mpos(state_copy2,hamiltonian_mpos,Ledge,Redge);
    double overlap_H = MPS_Tools::Finite::Ops::overlap(state_copy2,state);

    spdlog::trace("Computing expectation values");

    double energy      = MPS_Tools::Finite::Ops::expectation_value(state     ,state     ,hamiltonian_mpos, Ledge,Redge);
    double energy_up_x = MPS_Tools::Finite::Ops::expectation_value(state_up_x,state_up_x,hamiltonian_mpos, Ledge,Redge);
    double energy_dn_x = MPS_Tools::Finite::Ops::expectation_value(state_dn_x,state_dn_x,hamiltonian_mpos, Ledge,Redge);
    double energy_up_y = MPS_Tools::Finite::Ops::expectation_value(state_up_y,state_up_y,hamiltonian_mpos, Ledge,Redge);
    double energy_dn_y = MPS_Tools::Finite::Ops::expectation_value(state_dn_y,state_dn_y,hamiltonian_mpos, Ledge,Redge);
    double energy_up_z = MPS_Tools::Finite::Ops::expectation_value(state_up_z,state_up_z,hamiltonian_mpos, Ledge,Redge);
    double energy_dn_z = MPS_Tools::Finite::Ops::expectation_value(state_dn_z,state_dn_z,hamiltonian_mpos, Ledge,Redge);
    spdlog::trace("Computing variances");

    double variance      = MPS_Tools::Finite::Ops::exp_sq_value(state     ,state     ,hamiltonian_mpos, Ledge2,Redge2) - energy      * energy ;
    double variance_up_x = MPS_Tools::Finite::Ops::exp_sq_value(state_up_x,state_up_x,hamiltonian_mpos, Ledge2,Redge2) - energy_up_x * energy_up_x;
    double variance_dn_x = MPS_Tools::Finite::Ops::exp_sq_value(state_dn_x,state_dn_x,hamiltonian_mpos, Ledge2,Redge2) - energy_dn_x * energy_dn_x;
    double variance_up_y = MPS_Tools::Finite::Ops::exp_sq_value(state_up_y,state_up_y,hamiltonian_mpos, Ledge2,Redge2) - energy_up_y * energy_up_y;
    double variance_dn_y = MPS_Tools::Finite::Ops::exp_sq_value(state_dn_y,state_dn_y,hamiltonian_mpos, Ledge2,Redge2) - energy_dn_y * energy_dn_y;
    double variance_up_z = MPS_Tools::Finite::Ops::exp_sq_value(state_up_z,state_up_z,hamiltonian_mpos, Ledge2,Redge2) - energy_up_z * energy_up_z;
    double variance_dn_z = MPS_Tools::Finite::Ops::exp_sq_value(state_dn_z,state_dn_z,hamiltonian_mpos, Ledge2,Redge2) - energy_dn_z * energy_dn_z;
    std::cout << "Energy per site" << std::endl;

    std::cout << "<psi     | H/L |psi     >       = "  <<  energy    /state.get_length() << std::endl;
    std::cout << "<psi     | H/L  psi     >       = "  <<  overlap_H /state.get_length() << std::endl;

    std::cout << "<psi_up_x| H/L |psi_up_x>       = "  <<  energy_up_x / state.get_length() <<std::endl;
    std::cout << "<psi_dn_x| H/L |psi_dn_x>       = "  <<  energy_dn_x / state.get_length() <<std::endl;
    std::cout << "<psi_up_y| H/L |psi_up_y>       = "  <<  energy_up_y / state.get_length() <<std::endl;
    std::cout << "<psi_dn_y| H/L |psi_dn_y>       = "  <<  energy_dn_y / state.get_length() <<std::endl;
    std::cout << "<psi_up_z| H/L |psi_up_z>       = "  <<  energy_up_z / state.get_length() <<std::endl;
    std::cout << "<psi_dn_z| H/L |psi_dn_z>       = "  <<  energy_dn_z / state.get_length() <<std::endl;

    std::cout << "\nVariance per site" << std::endl;

    std::cout << "<psi     | (H2-E2)/L |psi     > = "  << variance / state.get_length()      << "  |   log10 =  "<< std::log10(variance / state.get_length())       << std::endl;
    std::cout << "<psi_up_x| (H2-E2)/L |psi_up_x> = "  << variance_up_x / state.get_length() << "  |   log10 =  "<< std::log10(variance_up_x / state.get_length())  << std::endl;
    std::cout << "<psi_dn_x| (H2-E2)/L |psi_dn_x> = "  << variance_dn_x / state.get_length() << "  |   log10 =  "<< std::log10(variance_dn_x / state.get_length())  << std::endl;
    std::cout << "<psi_up_y| (H2-E2)/L |psi_up_y> = "  << variance_up_y / state.get_length() << "  |   log10 =  "<< std::log10(variance_up_y / state.get_length())  << std::endl;
    std::cout << "<psi_dn_y| (H2-E2)/L |psi_dn_y> = "  << variance_dn_y / state.get_length() << "  |   log10 =  "<< std::log10(variance_dn_y / state.get_length())  << std::endl;
    std::cout << "<psi_up_z| (H2-E2)/L |psi_up_z> = "  << variance_up_z / state.get_length() << "  |   log10 =  "<< std::log10(variance_up_z / state.get_length())  << std::endl;
    std::cout << "<psi_dn_z| (H2-E2)/L |psi_dn_z> = "  << variance_dn_z / state.get_length() << "  |   log10 =  "<< std::log10(variance_dn_z / state.get_length())  << std::endl;
    std::cout << "\nMidchain entanglement entropies" << std::endl;
    spdlog::trace("Computing entanglement entropies");
    std::cout << "S(L/2) psi       = "  <<  MPS_Tools::Finite::Measure::midchain_entanglement_entropy(state     ) <<std::endl;
    std::cout << "S(L/2) psi_up_x  = "  <<  MPS_Tools::Finite::Measure::midchain_entanglement_entropy(state_up_x) <<std::endl;
    std::cout << "S(L/2) psi_dn_x  = "  <<  MPS_Tools::Finite::Measure::midchain_entanglement_entropy(state_dn_x) <<std::endl;
    std::cout << "S(L/2) psi_up_y  = "  <<  MPS_Tools::Finite::Measure::midchain_entanglement_entropy(state_up_y) <<std::endl;
    std::cout << "S(L/2) psi_dn_y  = "  <<  MPS_Tools::Finite::Measure::midchain_entanglement_entropy(state_dn_y) <<std::endl;
    std::cout << "S(L/2) psi_up_z  = "  <<  MPS_Tools::Finite::Measure::midchain_entanglement_entropy(state_up_z) <<std::endl;
    std::cout << "S(L/2) psi_dn_z  = "  <<  MPS_Tools::Finite::Measure::midchain_entanglement_entropy(state_dn_z) <<std::endl;



//    rebuild_environments(state_copy);
//
//    rebuild_superblock(state_copy, superblock);


}


void MPS_Tools::Finite::Ops::rebuild_superblock(class_finite_chain_state &state, class_superblock &superblock) {
    superblock.set_superblock(
            state.get_ENV2_L().back().block,
            state.get_ENV_L().back().block,
            state.get_MPO_L().back()->MPO,
            state.get_MPS_L().back().get_L(),
            state.get_MPS_L().back().get_G(),
            state.get_MPS_C(),
            state.get_MPS_R().front().get_G(),
            state.get_MPS_R().front().get_L(),
            state.get_MPO_R().front()->MPO,
            state.get_ENV_R().front().block,
            state.get_ENV2_R().front().block
    );
}


void MPS_Tools::Finite::Ops::rebuild_environments(class_finite_chain_state &state){

    // Generate new environments

    auto mpsL_it   = state.get_MPS_L().begin();
    auto mpoL_it   = state.get_MPO_L().begin();
    auto envL_it   = state.get_ENV_L().begin();
    auto env2L_it  = state.get_ENV2_L().begin();
    int i = 0;
    while(mpsL_it != state.get_MPS_L().end() and mpoL_it != state.get_MPO_L().end()) {
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
        if (mpsL_it == state.get_MPS_L().end()) break;
        envL_it->block = temp;
        env2L_it->block = temp2;
        i++;
    }

    auto mpsR_it   = state.get_MPS_R().rbegin();
    auto mpoR_it   = state.get_MPO_R().rbegin();
    auto envR_it   = state.get_ENV_R().rbegin();
    auto env2R_it  = state.get_ENV2_R().rbegin();
    while(mpsR_it != state.get_MPS_R().rend() and mpoR_it != state.get_MPO_R().rend()){
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
        if (mpsR_it == state.get_MPS_R().rend()) break;
        envR_it->block = temp;
        env2R_it->block = temp2;
    }

}


double MPS_Tools::Finite::Ops::overlap(const class_finite_chain_state &state1, const class_finite_chain_state &state2){

    assert(state1.get_length() == state2.get_length() and "ERROR: States have different lengths! Can't do overlap.");
    assert(state1.get_position() == state2.get_position() and "ERROR: States need to be at the same position! Can't do overlap.");

    Eigen::Tensor<Scalar,2> overlap =
            state1.get_MPS_L().front().get_A()
            .contract(state2.get_MPS_L().front().get_A().conjugate(), Textra::idx({0,1},{0,1}));
    auto mps1_it_L = state1.get_MPS_L().begin();
    auto mps2_it_L = state2.get_MPS_L().begin();
    mps1_it_L++;
    mps2_it_L++;

    while (mps1_it_L != state1.get_MPS_L().end() ){
        Eigen::Tensor<Scalar,2> temp = overlap
                .contract(mps1_it_L->get_A()            , Textra::idx({0},{1}))
                .contract(mps2_it_L->get_A().conjugate(), Textra::idx({0,1},{1,0}));
        overlap = temp;
        mps1_it_L++;
        mps2_it_L++;
    }

    Eigen::Tensor<Scalar,2> temp = overlap
            .contract(Textra::asDiagonal(state1.get_MPS_C()), Textra::idx({0},{0}))
            .contract(Textra::asDiagonal(state2.get_MPS_C()), Textra::idx({0},{0}));
    overlap = temp;



    auto mps1_it_R = state1.get_MPS_R().begin();
    auto mps2_it_R = state2.get_MPS_R().begin();
    while (mps1_it_R != state1.get_MPS_R().end() ){
        Eigen::Tensor<Scalar,2> temp = overlap
                .contract(mps1_it_R->get_B()            , Textra::idx({0},{1}))
                .contract(mps2_it_R->get_B().conjugate(), Textra::idx({0,1},{1,0}));
        overlap = temp;
        mps1_it_R++;
        mps2_it_R++;
    }
    double norm_chain = std::real(Textra::Tensor2_to_Matrix(overlap).trace());
//    std::cout << "Overlap state1 and state2: " << std::setprecision(16) << norm_chain << std::endl;
    return norm_chain;
}


//double MPS_Tools::Finite::Ops::expectation_value(
//        const class_finite_chain_state &state1,
//        const class_finite_chain_state &state2,
//        const Eigen::Tensor<std::complex<double>,4> mpo,
//        const Eigen::Tensor<std::complex<double>,3> Ledge,
//        const Eigen::Tensor<std::complex<double>,3> Redge)
//{
//    std::list<Eigen::Tensor<std::complex<double>,4>> mpos(state1.get_length(),mpo);
//    return expectation_value(state1,state2,mpos,Ledge,Redge);
//}

double MPS_Tools::Finite::Ops::expectation_value(const class_finite_chain_state &state1, const class_finite_chain_state &state2,const std::list<Eigen::Tensor<std::complex<double>,4>>  &mpos, const Eigen::Tensor<std::complex<double>,3> Ledge, const Eigen::Tensor<std::complex<double>,3> Redge){

    assert(state1.get_length() == state2.get_length() and "ERROR: States have different lengths! Can't do overlap.");
    assert(state1.get_position() == state2.get_position() and "ERROR: States need to be at the same position! Can't do overlap.");
    auto mps1_it_L = state1.get_MPS_L().begin();
    auto mps2_it_L = state2.get_MPS_L().begin();
    auto mpo_it    = mpos.begin();
    Eigen::Tensor<Scalar,3> L = Ledge;
    while(mps1_it_L != state1.get_MPS_L().end()){
        Eigen::Tensor<Scalar,3> temp =
                L
                .contract(mps1_it_L->get_A()                   , idx({0},{1}))
                .contract(*mpo_it                              , idx({1,2},{0,2}))
                .contract(mps2_it_L->get_A().conjugate()       , idx({0,3},{1,0}))
                .shuffle(array3{0,2,1});

        L = temp;
        mps1_it_L++;
        mps2_it_L++;
        mpo_it++;
    }

    {
        //Contract the center point
        auto &MPS_C1 = state1.get_MPS_C();
        auto &MPS_C2 = state2.get_MPS_C();
        Eigen::Tensor<Scalar,3> temp =
                L
                .contract(asDiagonal(MPS_C1), idx({0}, {0}))
                .contract(asDiagonal(MPS_C2), idx({0}, {0}))
                .shuffle(array3{1, 2, 0});
        L = temp;
    }
    //Contract the right half of the chain
    Eigen::Tensor<Scalar,3> R = Redge;
    auto mps1_it_R = state1.get_MPS_R().begin();
    auto mps2_it_R = state2.get_MPS_R().begin();
    while(mps1_it_R != state1.get_MPS_R().end()){
        Eigen::Tensor<Scalar,3> temp =
                L
                .contract(mps1_it_R->get_B()                   , idx({0},{1}))
                .contract(*mpo_it                              , idx({1,2},{0,2}))
                .contract(mps2_it_R->get_B().conjugate()       , idx({0,3},{1,0}))
                .shuffle(array3{0,2,1});
        L = temp;
        mps1_it_R++;
        mps2_it_R++;
        mpo_it++;
    }

    assert(L.dimensions() == R.dimensions());
    Eigen::Tensor<Scalar,0> E_all_sites = L.contract(R, idx({0,1,2},{0,1,2}));
    double energy_chain = std::real(E_all_sites(0));
//    std::cout << std::setprecision(16) << "E all sites: " << energy_chain << " | per site: " << energy_chain / state1.get_length() << std::endl;
    return energy_chain;
}

double MPS_Tools::Finite::Ops::exp_sq_value(const class_finite_chain_state &state1, const class_finite_chain_state &state2,const std::list<Eigen::Tensor<std::complex<double>,4>>  &mpos, const Eigen::Tensor<std::complex<double>,4> Ledge, const Eigen::Tensor<std::complex<double>,4> Redge){

    assert(state1.get_length() == state2.get_length() and "ERROR: States have different lengths! Can't do overlap.");
    assert(state1.get_position() == state2.get_position() and "ERROR: States need to be at the same position! Can't do overlap.");
    auto mps1_it_L = state1.get_MPS_L().begin();
    auto mps2_it_L = state2.get_MPS_L().begin();
    auto mpo_it    = mpos.begin();
    Eigen::Tensor<Scalar,4> L = Ledge;
    while(mps1_it_L != state1.get_MPS_L().end()){
        Eigen::Tensor<Scalar,4> temp =
                L
                .contract(mps1_it_L->get_A()                   , idx({0},{1}))
                .contract(*mpo_it                              , idx({1,3},{0,2}))
                .contract(*mpo_it                              , idx({1,4},{0,2}))
                .contract(mps2_it_L->get_A().conjugate()       , idx({0,4},{1,0}))
                .shuffle(array4{0,3,1,2});

        L = temp;
        mps1_it_L++;
        mps2_it_L++;
        mpo_it++;
    }

    {
        //Contract the center point
        auto &MPS_C1 = state1.get_MPS_C();
        auto &MPS_C2 = state2.get_MPS_C();
        Eigen::Tensor<Scalar,4> temp =
                L
                        .contract(asDiagonal(MPS_C1), idx({0}, {0}))
                        .contract(asDiagonal(MPS_C2), idx({0}, {0}))
                        .shuffle(array4{2,3,0,1});
        L = temp;
    }
    //Contract the right half of the chain
    Eigen::Tensor<Scalar,4> R = Redge;
    auto mps1_it_R = state1.get_MPS_R().begin();
    auto mps2_it_R = state2.get_MPS_R().begin();
    while(mps1_it_R != state1.get_MPS_R().end()){
        Eigen::Tensor<Scalar,4> temp =
                L
                .contract(mps1_it_R->get_B()                  , idx({0},{1}))
                .contract(*mpo_it                             , idx({1,3},{0,2}))
                .contract(*mpo_it                             , idx({1,4},{0,2}))
                .contract(mps2_it_R->get_B().conjugate()      , idx({0,4},{1,0}))
                .shuffle(array4{0,3,1,2});
        L = temp;
        mps1_it_R++;
        mps2_it_R++;
        mpo_it++;
    }

    assert(L.dimensions() == R.dimensions());
    Eigen::Tensor<Scalar,0> H2_all_sites = L.contract(R, idx({0,1,2,3},{0,1,2,3}));
    return std::real(H2_all_sites(0));
}
