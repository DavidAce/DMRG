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
    if (mpos.size() != state.get_length()) throw std::runtime_error("Number of mpo's doesn't match the number of sites on the system");
    // Apply MPO's on Gamma matrices and
    // increase the size on all Lambdas by chi*mpoDim
    spdlog::trace("Applying MPO's");
    // Print bond dimensions:
    //    for(auto & mps : state.get_MPS_L()){
    //        std::cout << "start chi : " << mps.get_chiL() << std::endl;
    //    }
    //    std::cout << "start chiC: " << state.get_MPS_C().dimension(0) << std::endl;
    //    for(auto & mps : state.get_MPS_R()){
    //        std::cout << "start chi : " << mps.get_chiR() << std::endl;
    //    }

    Eigen::Tensor<Scalar,1> I_tensor(1);
    I_tensor.setConstant(1);

    auto mpo = mpos.begin();
    {
        for (auto & mps : state.get_MPS_L()){
            long mpoDimL = mpo->dimension(0);
            long mpoDimR = mpo->dimension(1);
            auto [d,chiL,chiR] = mps.get_dims();

            Eigen::Tensor<Scalar,3> G_temp =
                    mps.get_G()
                            .contract(*mpo, idx({0},{2}))
                            .shuffle(array5{4,0,2,1,3})
                            .reshape(array3{d, chiL*mpoDimL, chiR*mpoDimR});

            mps.set_G(G_temp);

            Eigen::Tensor<Scalar,1> L_temp = mps.get_L().broadcast(array1{mpoDimL});
//            Eigen::Tensor<Scalar,1> L_temp = mps.get_L().contract(I_tensor.broadcast(array1{mpoDimL}), idx()).reshape(array1{chiL * mpoDimL});
//            Eigen::Tensor<Scalar,1> L_temp = I_tensor.broadcast(array1{mpoDimL}).contract(mps.get_L(), idx()).reshape(array1{chiL * mpoDimL});

            double Lnorm = Textra::Tensor1_to_Vector(L_temp).cwiseAbs2().sum();
            if (Lnorm < 1e-1){
                std::cout << "Ltemp\n" << L_temp<< std::endl;
                std::cout << "chiL = " << chiL << std::endl;
                std::cout << "mpoDimL = " << mpoDimL << std::endl;
                throw std::runtime_error("Schmidt values became too small when applying MPO");
            }
//            std::cout << "Ltemp\n" << L_temp<< std::endl;
            mps.set_L(L_temp);
//            std::cout << "applied chi  : " << mps.get_chiL() << std::endl;
            mpo++;
        }
    }

    //Take care of the middle
    {
        long mpoDim = mpo->dimension(0);
//        auto chiC = state.get_MPS_C().dimension(0);
//        Eigen::Tensor<Scalar,1> C_temp = state.get_MPS_C().contract(I_tensor.broadcast(array1{mpoDim}), idx()).reshape(array1{chiC * mpoDim});
//        Eigen::Tensor<Scalar,1> C_temp = I_tensor.broadcast(array1{mpoDim}).contract(state.get_MPS_C(), idx()).reshape(array1{chiC * mpoDim});

        Eigen::Tensor<Scalar, 1> C_temp = state.get_MPS_C().broadcast(array1{mpoDim});
        state.get_MPS_C() = C_temp;
//        std::cout << "applied chiC : " << state.get_MPS_C().dimension(0) << std::endl;

    }
    {
        for (auto & mps : state.get_MPS_R()){
            long mpoDimL = mpo->dimension(0);
            long mpoDimR = mpo->dimension(1);
            auto [d,chiL,chiR] = mps.get_dims();
            Eigen::Tensor<Scalar,3> G_temp =
                    mps.get_G()
                            .contract(*mpo, idx({0},{2}))
                            .shuffle(array5{4,0,2,1,3})
                            .reshape(array3{d, chiL*mpoDimL, chiR*mpoDimR});

            mps.set_G(G_temp);
            Eigen::Tensor<Scalar,1> L_temp = mps.get_L().broadcast(array1{mpoDimR});
//            Eigen::Tensor<Scalar,1> L_temp = mps.get_L().contract(I_tensor.broadcast(array1{mpoDimR}), idx()).reshape(array1{chiR * mpoDimR});
//            Eigen::Tensor<Scalar,1> L_temp = I_tensor.broadcast(array1{mpoDimR}).contract(mps.get_L(), idx()).reshape(array1{chiR * mpoDimR});
            double Lnorm = Textra::Tensor1_to_Vector(L_temp).cwiseAbs2().sum();
            if (Lnorm < 1e-1) {
                std::cout << "Ltemp\n" << L_temp<< std::endl;
                std::cout << "chiL = " << chiL << std::endl;
                std::cout << "mpoDimL = " << mpoDimR << std::endl;
                throw std::runtime_error("Schmidt values became too small when applying MPO");
            }
            mps.set_L(L_temp);
//            std::cout << "applied chi  : " << mps.get_chiR() << std::endl;
//            std::cout << "Ltemp\n" << L_temp<< std::endl;
            mpo++;
        }
    }

    // Take care of the edges. Apply the left and right MPO-edges on A's and B's
    // so the left and right most lambdas become ones
    {
        long mpoDimL = mpos.front().dimension(0);
        auto Ldim = Ledge.dimension(0);
        Eigen::Tensor<Scalar, 3> G_temp =
                Ledge
                        .shuffle(Textra::array3{0,2,1})
                        .reshape(Textra::array2{Ldim*mpoDimL,Ldim})
                        .contract(state.get_MPS_L().front().get_A(),Textra::idx({0},{1}))
                        .shuffle(Textra::array3{1,0,2});
        state.get_MPS_L().front().set_G(G_temp);
        state.get_MPS_L().front().set_L(Eigen::Tensor<Scalar,1>(Ldim).constant(1.0));
    }
    {
        long mpoDimR = mpos.back().dimension(1);
        auto Rdim = Redge.dimension(0);
        Eigen::Tensor<Scalar, 3> G_temp =
                Redge
                        .shuffle(Textra::array3{0,2,1})
                        .reshape(Textra::array2{Rdim*mpoDimR,Rdim})
                        .contract(state.get_MPS_R().back().get_B(),Textra::idx({0},{2}))
                        .shuffle(Textra::array3{1,2,0});
        state.get_MPS_R().back().set_G(G_temp);
        state.get_MPS_R().back().set_L(Eigen::Tensor<Scalar,1>(Rdim).constant(1.0));
    }


}

void MPS_Tools::Finite::Ops::normalize_chain(class_finite_chain_state & state){

    // First the left side
    spdlog::trace("Starting normalization");
    class_SVD<Scalar> svd;
//    svd.setThreshold(1e-8);
    {
        auto mpsL_0            = state.get_MPS_L().begin();
        auto mpsL_1            = state.get_MPS_L().begin();
        auto mpsL_2            = state.get_MPS_L().begin();
        std::advance(mpsL_1,1);
        std::advance(mpsL_2,2);
        // Make sure the first lambda is a 1
        Eigen::Tensor<Scalar,1> L0(1);
        L0.setConstant(1.0);
        mpsL_0->set_L(L0);

        while(mpsL_1 != state.get_MPS_L().end()) {
            const auto & L0 = mpsL_0->get_L();
            const auto & A0 = mpsL_0->get_A();
            const auto & A1 = mpsL_1->get_A();
            const auto & L2 = mpsL_2 != state.get_MPS_L().end() ? mpsL_2->get_L() : state.get_MPS_C();
            Eigen::Tensor<Scalar,4> theta = A0.contract(A1, idx({2},{1})).contract(Textra::asDiagonal(L2), idx({3},{0}));
            auto [U,S,V] = svd.schmidt(theta);
            Eigen::Tensor<Scalar,3> L_U = asDiagonalInversed(L0).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
            Eigen::Tensor<Scalar,3> V_L = V.contract(asDiagonalInversed(L2), idx({2},{0}));

            // Now you can update your matrices
            mpsL_0->set_G(L_U);
            mpsL_1->set_L(S);
            mpsL_1->set_G(V_L);

            mpsL_0++;
            mpsL_1++;
            if(mpsL_2 != state.get_MPS_L().end()){mpsL_2++;}
        }
    }

    {
        auto mpsR_0            = state.get_MPS_R().rbegin();
        auto mpsR_1            = state.get_MPS_R().rbegin();
        auto mpsR_2            = state.get_MPS_R().rbegin();
        std::advance(mpsR_1,1);
        std::advance(mpsR_2,2);
        // Make sure the first lambda is a 1
        Eigen::Tensor<Scalar,1> L0(1);
        L0.setConstant(1.0);
        mpsR_0->set_L(L0);

        while(mpsR_1 != state.get_MPS_R().rend()) {
            const auto & L0 = mpsR_0->get_L();
            const auto & B0 = mpsR_0->get_B();
            const auto & B1 = mpsR_1->get_B();
            const auto & L2 = mpsR_2 != state.get_MPS_R().rend() ? mpsR_2->get_L() : state.get_MPS_C();
            Eigen::Tensor<Scalar,4> theta =B0.contract(B1, idx({1},{2})).contract(Textra::asDiagonal(L2), idx({3},{1})).shuffle(array4{2,3,0,1});
            auto [U,S,V] = svd.schmidt(theta);
            Eigen::Tensor<Scalar,3> L_U = asDiagonalInversed(L2).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
            Eigen::Tensor<Scalar,3> V_L = V.contract(asDiagonalInversed(L0), idx({2},{0}));

            // Now you can update your matrices
            mpsR_0->set_G(V_L);
            mpsR_1->set_L(S);
            mpsR_1->set_G(L_U);

            mpsR_0++;
            mpsR_1++;
            if(mpsR_2 != state.get_MPS_R().rend()){mpsR_2++;}
        }

    }

    // Now take care of the middle
    {
        const auto & A = state.get_MPS_L().back().get_A();
        const auto & C = Textra::asDiagonal(state.get_MPS_C());
        const auto & B = state.get_MPS_R().front().get_B();
        Eigen::Tensor<Scalar,4> theta = A.contract(C,idx({2},{0})).contract(B, idx({2},{1}));
        auto [U,S,V] = svd.schmidt(theta);
        state.get_MPS_L().back().set_G (Textra::asDiagonalInversed(state.get_MPS_L().back().get_L()).contract(U, idx({1},{1})).shuffle(array3{1,0,2}));
        state.get_MPS_R().front().set_G(V.contract(Textra::asDiagonalInversed(state.get_MPS_R().front().get_L()), idx({2},{0})));
        state.get_MPS_C() = S;
    }

    // Print bond dimensions:
    //    for(auto & mps : state.get_MPS_L()){
    //        std::cout << "chi : " << mps.get_chiL() << std::endl;
    //    }
    //    std::cout << "chiC: " << state.get_MPS_C().dimension(0) << std::endl;
    //    for(auto & mps : state.get_MPS_R()){
    //        std::cout << "chi : " << mps.get_chiR() << std::endl;
    //    }
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

class_finite_chain_state MPS_Tools::Finite::Ops::get_parity_projected_state(const class_finite_chain_state &state, const Eigen::MatrixXcd  paulimatrix, const int sign) {
    if (std::abs(sign) != 1) throw std::runtime_error("Expected 'sign' +1 or -1. Got: " + std::to_string(sign));
    spdlog::trace("Generating parity projected state");
    class_finite_chain_state state_projected = state;
    const auto [mpo,L,R]    = class_mpo::parity_projector_mpos(paulimatrix,state_projected.get_length(), sign);
    apply_mpos(state_projected,mpo, L,R);
    normalize_chain(state_projected);
    rebuild_environments(state_projected);
    MPS_Tools::Finite::Debug::check_integrity_of_mps(state_projected);
    return state_projected;
}



void MPS_Tools::Finite::Ops::rebuild_superblock(class_finite_chain_state &state, class_superblock &superblock) {
    MPS_Tools::Finite::Chain::copy_state_to_superblock(state,superblock);
}

void MPS_Tools::Finite::Ops::rebuild_environments(class_finite_chain_state &state){

    // Generate new environments
    assert(not state.get_MPS_L().empty() and "ERROR: The MPS L list is empty:");
    assert(not state.get_MPS_R().empty() and "ERROR: The MPS R list is empty:");
    assert(not state.get_MPO_L().empty() and "ERROR: The MPO L list is empty:");
    assert(not state.get_MPO_R().empty() and "ERROR: The MPO R list is empty:");


    state.get_ENV_L().clear();
    state.get_ENV2_L().clear();
    state.get_ENV_R().clear();
    state.get_ENV2_R().clear();

    {
        auto ENV_L = class_environment("L",state.get_MPS_L().front().get_chiL(), state.get_MPO_L().front()->MPO.dimension(0));
        auto ENV2_L = class_environment_var("L",state.get_MPS_L().front().get_chiL(), state.get_MPO_L().front()->MPO.dimension(0));
        auto mpsL_it   = state.get_MPS_L().begin();
        auto mpoL_it   = state.get_MPO_L().begin();

        while(mpsL_it != state.get_MPS_L().end() and mpoL_it != state.get_MPO_L().end()) {
            state.get_ENV_L().push_back(ENV_L);
            state.get_ENV2_L().push_back(ENV2_L);
            if (mpsL_it == state.get_MPS_L().end()) break;
            ENV_L.enlarge(mpsL_it->get_A(), mpoL_it->get()->MPO);
            ENV2_L.enlarge(mpsL_it->get_A(), mpoL_it->get()->MPO);
            mpsL_it++;
            mpoL_it++;
        }
    }



    {
        auto ENV_R  = class_environment("R",state.get_MPS_R().back().get_chiR(), state.get_MPO_L().back()->MPO.dimension(1));
        auto ENV2_R = class_environment_var("R",state.get_MPS_R().back().get_chiR(), state.get_MPO_L().back()->MPO.dimension(1));
        auto mpsR_it   = state.get_MPS_R().rbegin();
        auto mpoR_it   = state.get_MPO_R().rbegin();
        while(mpsR_it != state.get_MPS_R().rend() and mpoR_it != state.get_MPO_R().rend()){
            state.get_ENV_R().push_front(ENV_R);
            state.get_ENV2_R().push_front(ENV2_R);
            if (mpsR_it == state.get_MPS_R().rend()) break;
            ENV_R.enlarge(mpsR_it->get_B(), mpoR_it->get()->MPO);
            ENV2_R.enlarge(mpsR_it->get_B(), mpoR_it->get()->MPO);
            mpsR_it++;
            mpoR_it++;
        }
    }


    state.set_positions();
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
