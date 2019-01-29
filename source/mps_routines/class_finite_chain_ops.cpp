//
// Created by david on 2019-01-21.
//

#include "class_finite_chain.h"
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mpo.h>
#include <general/nmspc_tensor_extra.h>
#include <general/class_svd_wrapper.h>
#include <general/nmspc_quantum_mechanics.h>

#include <iomanip>
using namespace Textra;

void class_finite_chain::apply_mpo(const Eigen::Tensor<Scalar,4> mpo,const Eigen::Tensor<Scalar,3> Ledge, const Eigen::Tensor<Scalar,3> Redge){
    // First the left side
    auto mps            = MPS_L.begin();
    auto mpsR_it            = MPS_R.rbegin();
    long mpoDim = mpo.dimension(0);
    auto I_tensor1 = Eigen::Tensor<Scalar,1>(mpoDim);
    I_tensor1.setConstant(1);


    // Apply MPO's on Gamma matrices and
    // increase the size on all Lambdas except at the edges

    for (auto & mps : MPS_L){
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

    for (auto & mps : MPS_R){
        if (mps.get_position() == max_length - 1) continue;
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
        auto[d, chiL, chiR] = MPS_L.front().get_dims();
        Eigen::Tensor<Scalar, 3> G_temp =
                MPS_L.front().get_G()
                        .contract(mpo, idx({0}, {2}))
                        .shuffle(array5{4, 2, 0, 3, 1})
                        .contract(Ledge, idx({1, 2}, {2, 0}))
                        .reshape(array3{d, mpoDim*chiR, Ledge.dimension(1)})
                        .shuffle(array3{0, 2, 1});
        MPS_L.front().set_G(G_temp);

    }
    {
        auto[d, chiL, chiR] = MPS_R.back().get_dims();
        Eigen::Tensor<Scalar, 3> G_temp =
                MPS_R.back().get_G()
                        .contract(mpo, idx({0}, {2}))
                        .shuffle(array5{4, 2, 0, 3, 1})
                        .contract(Redge, idx({3, 4}, {2, 0}))
                        .reshape(array3{d, mpoDim*chiL, Redge.dimension(1)});
        MPS_R.back().set_G(G_temp);
    }


    //Take care of the middle
    auto chiC = MPS_C.dimension(0);
    Eigen::Tensor<Scalar,1> C_temp = MPS_C.contract(I_tensor1, idx()).reshape(array1{chiC * mpoDim});
    MPS_C = C_temp;


}


void class_finite_chain::normalize_chain(){

    // First the left side
//    std::cout << "Starting Normalization" << std::endl;
    class_SVD<Scalar> svd;
    Eigen::Tensor<Scalar,2> C_temp;
    auto mpsL_it            = MPS_L.begin();
    Eigen::Tensor<Scalar,3> A_temp = mpsL_it->get_A();

    while(mpsL_it != MPS_L.end()) {
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
        if (mpsL_it == MPS_L.end()){
            C_temp = Textra::asDiagonal(S)
                    .contract(V               , idx({1},{0}))
                    .contract(Textra::asDiagonal(MPS_C), idx({1},{0}));
            break;
        }
        A_temp =
                Textra::asDiagonal(S)
                .contract(V               , idx({1},{0}))
                .contract(mpsL_it->get_A(), idx({1},{1}))
                .shuffle(array3{1,0,2});

        mpsL_it->set_L(S);
    }


    auto mpsR_it            = MPS_R.rbegin();
    Eigen::Tensor<Scalar,3> B_temp = mpsR_it->get_B().shuffle(array3{1,0,2});
    while(mpsR_it != MPS_R.rend()) {
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
        if (mpsR_it == MPS_R.rend()){
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
            MPS_L.back().get_A()
            .contract(C_temp , idx({2},{0}))
            .contract(MPS_R.front().get_B(), idx({2},{1}));
    auto [U,S,V] = svd.schmidt(ACB);
    MPS_L.back().set_G (Textra::asDiagonalInversed(MPS_L.back().get_L()).contract(U, idx({1},{1})).shuffle(array3{1,0,2}));
    MPS_R.front().set_G(V.contract(Textra::asDiagonalInversed(MPS_R.front().get_L()), idx({2},{0})));
    MPS_C = S;


}


void class_finite_chain::set_parity_projected_mps(const Eigen::MatrixXcd  paulimatrix) {
//    std::cout << "Dimensions before projecting \n";
//    for (auto & mps : MPS_L) std::cout << "MPS_L dims: " << mps.get_G().dimensions()<< std::endl;
//    std::cout << "MPS_C dims: " << MPS_C.dimensions() << std::endl;
//    for (auto & mps : MPS_R) std::cout << "MPS_R dims: " << mps.get_G().dimensions()<< std::endl;



    const auto [mpo,L,R]    = class_mpo::parity_selector_mpo(paulimatrix);
//    const auto [mpo,L,R]    = class_mpo::pauli_mpo(paulimatrix);

    apply_mpo(mpo,L,R);

    normalize_chain();

    refresh_environments();
    refresh_superblock();


}


void class_finite_chain::refresh_superblock() {
    superblock->Lblock->block = ENV_L.back().block;
    superblock->Rblock->block = ENV_R.front().block;
    superblock->Lblock2->block = ENV2_L.back().block;
    superblock->Rblock2->block = ENV2_R.front().block;
    superblock->MPS->MPS_A->set_mps(MPS_L.back().get_G(), MPS_L.back().get_L());
    superblock->MPS->MPS_B->set_mps(MPS_R.front().get_G(), MPS_R.front().get_L());
    superblock->MPS->LC = MPS_C;
}


void class_finite_chain::set_parity_projected_mps2(const Eigen::MatrixXcd  paulimatrix) {
    std::cout << "Dimensions before projecting \n";
    for (auto & mps : MPS_L) std::cout << "MPS_L dims: " << mps.get_G().dimensions()<< std::endl;
    std::cout << "MPS_C dims: " << MPS_C.dimensions() << std::endl;
    for (auto & mps : MPS_R) std::cout << "MPS_R dims: " << mps.get_G().dimensions()<< std::endl;
    class_SVD<Scalar> svd;
    auto mpsL_it            = MPS_L.begin();
    auto mpsR_it            = MPS_R.rbegin();
    const auto [mpo,L,R]    = class_mpo::pauli_mpo(paulimatrix);
    long mpoDimL = mpo.dimension(0);
    long mpoDimR = mpo.dimension(1);


    // First the left side

    while(mpsL_it != MPS_L.end()) {
        auto [d,chiL,chiR] = mpsL_it->get_dims();
        std::cout << "Current Lambda: \n" << mpsL_it->get_L() << std::endl;
        std::cout << "Current dims L: [ " << d << " " << chiL << " " << chiR << " ]"<< std::endl;
        Eigen::Tensor<Scalar,3> temp = mpsL_it->get_A()
                .contract(mpo, idx({0},{2}))
                .shuffle(array5{4,2,0,3,1})
                .reshape(array3{d, chiL*mpoDimL, chiR*mpoDimR});
        std::cout << "New     dims L: " << temp.dimensions() << std::endl;

        auto [U,S,V] = svd.decompose(temp,temp.dimension(0)*temp.dimension(1),temp.dimension(2));
        std::cout << "  U dims: " << U.dimensions() << std::endl;
        std::cout << "  S dims: " << S.dimensions() << std::endl;
        std::cout << "  V dims: " << V.dimensions() << std::endl;
        long chiC = S.dimension(0);
        temp = Textra::asDiagonalInversed(mpsL_it->get_L())
                .contract(U.reshape(Textra::array3{d,chiL,chiC}), idx({1},{1}))
                .shuffle(array3{1,0,2});
        mpsL_it->set_G(temp);
        mpsL_it++;
        if (mpsL_it == MPS_L.end()) break;
        temp =  Textra::asDiagonal(S)
                .contract(V               , idx({1},{0}))
                .contract(mpsL_it->get_A(), idx({1},{1}))
                .shuffle(array3{1,0,2});
        mpsL_it->set_L(S);
        mpsL_it->set_G(temp);
    }

    // Check that the matrices in MPS_L are left-unitary
    for (auto & mps : MPS_L){
        Eigen::Tensor<Scalar,2> ID = mps.get_A().contract(mps.get_A().conjugate(), idx({0,1},{0,1}));
        std::cout << "ID "<< mps.get_position() << "\n" << ID << std::endl;
    }



    // Increase the size of MPS_C (the center bond) by doing an outer product
    // with an identity matrix that has the same dimension as the MPO, contracting over no indices,
    // mpsC       =  0 ------ 1
    // I_matrix   =  0 ------ 1
    // C_matrix   = 01 ====== 23
    long mpsCdim = MPS_C.dimension(0);
    auto I_tensor = Textra::Matrix_to_Tensor2(Eigen::MatrixXcd::Identity(mpoDimL,mpoDimR));
    Eigen::Tensor<Scalar,2> C_matrix = Textra::asDiagonal(MPS_C).contract(I_tensor, idx()).shuffle(array4{0,2,1,3}).reshape(array2{mpsCdim*mpoDimL, mpsCdim*mpoDimR});
    std::cout << std::setprecision(16) << C_matrix.real() << std::endl;



    // Now the right side
    while(mpsR_it != MPS_R.rend()) {
        auto [d,chiL,chiR] = mpsR_it->get_dims();
        std::cout << "Current Lambda: \n" << mpsR_it->get_L() << std::endl;
        std::cout << "Current dims R: [ " << d << " " << chiL << " " << chiR << " ]"<< std::endl;
        Eigen::Tensor<Scalar,3> temp = mpsR_it->get_B()
                .contract(mpo, idx({0},{2}))
                .shuffle(array5{4,2,0,3,1})
                .reshape(array3{d, chiL*mpoDimL, chiR*mpoDimR})
                .shuffle(array3{1,0,2});
        std::cout << "New     dims R: " << temp.dimensions() << std::endl;
        auto [U,S,V] = svd.decompose(temp,temp.dimension(0),temp.dimension(1)*temp.dimension(2));
        std::cout << "Current S: \n" << S << std::endl;
        std::cout << "  U dims: " << U.dimensions() << std::endl;
        std::cout << "  S dims: " << S.dimensions() << std::endl;
        std::cout << "  V dims: " << V.dimensions() << std::endl;
        long chiC = S.dimension(0);

        temp = V.reshape(Textra::array3{d,chiC,chiR})
                .contract(Textra::asDiagonalInversed(mpsR_it->get_L()), idx({2},{0}));
        mpsR_it->set_G(temp);
        mpsR_it++;
        if (mpsR_it == MPS_R.rend()) break;
        temp = mpsR_it->get_B()
                .contract(U, idx({2},{0}))
                .contract(Textra::asDiagonal(S),idx({2},{0}));

        mpsR_it->set_L(S);
        mpsR_it->set_G(temp);
    }


    // Check that the matrices in MPS_R are left-unitary
    for (auto & mps : MPS_R){
        Eigen::Tensor<Scalar,2> ID = mps.ref_B().contract(mps.ref_B().conjugate(), idx({0,2},{0,2}));
        std::cout << "ID "<< mps.get_position() << "\n" << ID << std::endl;
    }


    // Now we just have the middle part left
    Eigen::Tensor<Scalar,4> theta =
            MPS_L.back().get_A()
                    .contract(C_matrix, idx({2},{0}))
                    .contract(MPS_R.front().get_B(), idx({2},{1}));
    auto [U,S,V] = svd.schmidt(theta);
    MPS_L.back().set_G(Textra::asDiagonalInversed(MPS_L.back().get_L()).contract(U, idx({1},{1})).shuffle(array3{1,0,2}));
    MPS_R.front().set_G(V.contract(Textra::asDiagonalInversed(MPS_R.front().get_L()), idx({2},{0})));
    MPS_C = S;


    std::cout << "Dimensions after projecting \n";
    for (auto & mps : MPS_L) std::cout << "MPS_L dims: " << mps.get_G().dimensions()<< std::endl;
    std::cout << "MPS_C dims: " << MPS_C.dimensions() << std::endl;
    for (auto & mps : MPS_R) std::cout << "MPS_R dims: " << mps.get_G().dimensions()<< std::endl;


}




void class_finite_chain::set_parity_projected_mps3(const Eigen::MatrixXcd  paulimatrix){
    std::cout << "Dimensions before projecting \n ";
    for (auto & mps : MPS_L) std::cout << "MPS_L dims: " << mps.get_G().dimensions()<< std::endl;
    std::cout << "MPS_C dims: " << MPS_C.dimensions() << std::endl;
    for (auto & mps : MPS_R) std::cout << "MPS_R dims: " << mps.get_G().dimensions()<< std::endl;


    long A_length           = get_MPS_L().size();
    long B_length           = get_MPS_R().size();

    // Get handles for the current mps
    const auto &mpsL        = get_MPS_L();
    const auto &mpsC        = get_MPS_C();
    const auto &mpsR        = get_MPS_R();
    const auto [mpo,L,R]    = class_mpo::pauli_mpo(paulimatrix);

    // Make containers for the new updated mps
    std::vector<Eigen::Tensor<Scalar,3>> A_matrices(A_length);
    std::vector<Eigen::Tensor<Scalar,3>> B_matrices(B_length);
    Eigen::Tensor<Scalar,2> C_matrix;

    // Start updating  A matrices
    int pos = 0;
    long mpoDimL = mpo.dimension(0);
    long mpoDimR = mpo.dimension(1);
    for (auto &mps : mpsL){
//        assert(pos == mps.get_position() and "ERROR: Position and pos do not match");
        long d    = mps.get_spin_dim();
        long chiL = mps.get_chiL();
        long chiR = mps.get_chiR();
        A_matrices[pos++] =
                mps.ref_A()
                        .contract(mpo, idx({0},{2}))
                        .shuffle(array5{4,2,0,3,1})
                        .reshape(array3{d, chiL*mpoDimL, chiR*mpoDimR});
    }
    // Increase the size of mpsC (the center bond) by doing an outer product
    // with an identity matrix that has the same dimension as the MPO, contracting over no indices,
    // mpsC       =  0 ------ 1
    // I_matrix   =  0 ------ 1
    // C_matrix   = 01 ====== 23
    long mpsCdim = mpsC.dimension(0);
    long mpodim = mpo.dimension(0);
    auto I_tensor = Textra::Matrix_to_Tensor2(Eigen::MatrixXcd::Identity(mpodim,mpodim));
    C_matrix = Textra::asDiagonal(mpsC).contract(I_tensor, idx()).shuffle(array4{0,2,1,3}).reshape(array2{mpsCdim*mpodim, mpsCdim*mpodim});

    // Start updating  B matrices
//    pos = 0;
    int i = 0;
    for (auto &mps : mpsR){
//        assert(pos == mps.get_position() and "ERROR: Position and pos do not match");
        long d    = mps.get_spin_dim();
        long chiL = mps.get_chiL();
        long chiR = mps.get_chiR();
        B_matrices[i++] =
                mps.ref_B()
                        .contract(mpo, idx({0},{2}))
                        .shuffle(array5{4,2,0,3,1})
                        .reshape(array3{d, chiL*mpoDimL, chiR*mpoDimR});
    }

    // Now we have updated versions of all matrices. We need to split them
    // up into Gamma Lambda form, and because the updated matrices do not yield
    // a normalized wave function, we need to do SVD from the edges of the chain
    std::list<class_vidal_mps> mpsL_new;
    std::list<class_vidal_mps> mpsR_new;
    Eigen::Tensor<Scalar,1>    mpsC_new;

    class_SVD<Scalar> svd;
//    svd.setThreshold(1e-8);
    Eigen::Tensor<Scalar,3> A_next;
    Eigen::Tensor<Scalar,3> B_next;
    Eigen::Tensor<Scalar,1> LambdaA_next(1);
    Eigen::Tensor<Scalar,1> LambdaB_next(1);
    LambdaA_next.setConstant(1);
    LambdaB_next.setConstant(1);
    for (size_t pos = 0; pos < A_matrices.size()-1; pos++){
        long d    = A_matrices[pos].dimension(0);
        long chiL = A_matrices[pos].dimension(1);
        long chiR = A_matrices[pos].dimension(2);
        long rows = d*chiL;
        long cols = chiR;
        auto [U,S,V] = svd.decompose(A_matrices[pos],rows,cols);
        long chiC = S.size();
        A_next = Textra::asDiagonalInversed(LambdaA_next).contract(U.reshape(Textra::array3{d,chiL,chiC}), idx({1},{1})).shuffle(array3{1,0,2});
        mpsL_new.emplace_back(A_next,LambdaA_next,pos);

        A_next = Textra::asDiagonalInversed(S)
                .contract(V, idx({1},{0}))
                .contract(A_matrices[pos+1], idx({1},{1}))
                .shuffle(array3{1,0,2});

        A_matrices[pos+1] = A_next;
        LambdaA_next = S;
    }

    for (size_t pos = B_matrices.size()-1; pos > 0; pos--){
        long d    = B_matrices[pos].dimension(0);
        long chiL = B_matrices[pos].dimension(1);
        long chiR = B_matrices[pos].dimension(2);
        long rows = chiL;
        long cols = d*chiR;
        auto [U,S,V] = svd.decompose(B_matrices[pos],rows,cols);
        long chiC = S.size();
        B_next = V.reshape(Textra::array3{d,chiC,chiR}).contract(Textra::asDiagonalInversed(LambdaB_next), idx({2},{0}));
        mpsL_new.emplace_front(B_next,LambdaB_next,pos);
        B_next = B_matrices[pos-1]
                .contract(U, idx({2},{0}))
                .contract(Textra::asDiagonalInversed(S),idx({2},{0}));
        B_matrices[pos-1] = B_next;
        LambdaB_next = S;
    }

    // We have now worked our way from the edges of the chain in to the middle,
    // so we have something of the form
    // L G L G A_{i} C B_{i+1} G L G L
    //
    // The middle chunk A C B needs a 2-site SVD, and then
    // G_i     = LambdaA_next^-1  U
    // C       = S
    // G_i+1   = V LambdaB_next^-1
    pos = A_matrices.size();
    Eigen::Tensor<Scalar,4> ACB = A_matrices.back()
            .contract(C_matrix, idx({2},{0}))
            .contract(B_matrices.front(), idx({2},{1}));
    auto [U,S,V] = svd.schmidt(ACB);
    mpsL_new.emplace_back (Textra::asDiagonalInversed(LambdaA_next).contract(U, idx({1},{1})).shuffle(array3{1,0,2}),LambdaA_next, pos);
    mpsR_new.emplace_front(V.contract(Textra::asDiagonalInversed(LambdaB_next), idx({2},{0})), LambdaB_next, pos+1);
    mpsC_new = S;

    // Now update the original MPS in env_storage
    MPS_L = mpsL_new;
    MPS_C = mpsC_new;
    MPS_R = mpsR_new;


    std::cout << "Dimensions after projecting \n ";
    for (auto & mps : MPS_L) std::cout << "MPS_L dims: " << mps.get_G().dimensions()<< std::endl;
                             std::cout << "MPS_C dims: " << MPS_C.dimensions() << std::endl;
    for (auto & mps : MPS_R) std::cout << "MPS_R dims: " << mps.get_G().dimensions()<< std::endl;


    // Environments need updating now
    refresh_environments();




//    assert(L.dimensions() == R.dimensions());
//    Eigen::Tensor<Scalar,0> parity_tmp = L.contract(R, idx({0,1,2},{0,1,2}));
//    double parity = std::real(parity_tmp(0));
//    std::cout << setprecision(16) << "Parity: " << parity_tmp << std::endl;
}



void class_finite_chain::refresh_environments(){

    // Generate new environments

    auto mpsL_it  = MPS_L.begin();
    auto mpoL_it  = MPO_L.begin();
    auto envL_it  = ENV_L.begin();
    auto env2L_it  = ENV2_L.begin();
    int i = 0;
    while(mpsL_it != MPS_L.end() and mpoL_it != MPO_L.end()) {
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
        if (mpsL_it == MPS_L.end()) break;
        envL_it->block = temp;
        env2L_it->block = temp2;
        i++;
    }

    auto mpsR_it  = MPS_R.rbegin();
    auto mpoR_it  = MPO_R.rbegin();
    auto envR_it  = ENV_R.rbegin();
    auto env2R_it  = ENV2_R.rbegin();
    while(mpsR_it != MPS_R.rend() and mpoR_it != MPO_R.rend()){
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
        if (mpsR_it == MPS_R.rend()) break;
        envR_it->block = temp;
        env2R_it->block = temp2;
    }

}
