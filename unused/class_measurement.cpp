//
// Created by david on 2017-11-12.
//



#include <iomanip>
#include <complex>
#include <mps_routines/class_measurement.h>
#include <mps_routines/nmspc_mps_tools.h>
#include <general/nmspc_tensor_extra.h>
#include <general/class_svd_wrapper.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_mps_2site.h>
#include <mps_routines/class_finite_chain_state.h>
#include <mps_routines/class_mpo.h>
#include <general/nmspc_quantum_mechanics.h>
using namespace std;
using namespace Textra;




class_measurement::class_measurement(SimulationType sim_)
        : sim_type(sim_)
{
//    mps_util = std::make_unique<class_mps_util>();
    set_profiling_labels();

}




void class_measurement::compute_all_observables_from_superblock(class_superblock & superblock)
{
    if (is_measured){return;}
//    auto theta = MPS_Tools::Common::Views::get_theta(superblock);
    switch(sim_type){
        case SimulationType::iDMRG:
            compute_energy_mpo(superblock);
            compute_energy_variance_mpo(superblock);
            if (superblock.MPS->chiA() != superblock.MPS->chiB()) { return;}
            if (superblock.MPS->chiA() != superblock.MPS->chiC()) { return;}
            MPS_Tools::Common::Views::compute_mps_components(superblock);
            compute_energy_ham(superblock);
            compute_energy_variance_ham(superblock);
            if (superblock.MPS->chiC() < 2) { break; }
            compute_energy_and_variance_mom(superblock,a);
            break;
        case SimulationType::xDMRG:
        case SimulationType::fDMRG:
            compute_energy_mpo(superblock);
            compute_energy_variance_mpo(superblock);
            break;
        case SimulationType::iTEBD:
            if (superblock.MPS->chiA() != superblock.MPS->chiB()) { return;}
            if (superblock.MPS->chiA() != superblock.MPS->chiC()) { return;}
            MPS_Tools::Common::Views::compute_mps_components(superblock);
            compute_energy_ham(superblock);
            compute_energy_variance_ham(superblock);
            if (superblock.MPS->chiC() < 2) { break; }
            compute_energy_and_variance_mom(superblock,a);
            break;
    }

    compute_entanglement_entropy(superblock);
    is_measured = true;

}


void class_measurement::compute_all_observables_from_state(class_finite_chain_state & state){
    compute_finite_chain_energy(state);
    compute_finite_chain_energy_variance(state);

    if(sim_type == SimulationType::xDMRG and settings::xdmrg::store_wavefn){
        compute_finite_chain_mps_wavefn(state); // Runs out of memory for L > 20 or so, because the wave vector is doubled in size for each additional site.
    }
    compute_finite_chain_norm(state);
    compute_finite_chain_norm2(state);
    parity_x = compute_parity(state,qm::spinOneHalf::sx);
    parity_y = compute_parity(state,qm::spinOneHalf::sy);
    parity_z = compute_parity(state,qm::spinOneHalf::sz);
}




void class_measurement::compute_energy_mpo(const class_superblock & superblock){
    t_ene_mpo.tic();
    auto theta = MPS_Tools::Common::Views::get_theta(superblock);
    double L   = superblock.Lblock->size + superblock.Rblock->size + 2;
    Eigen::Tensor<Scalar, 0>  E =
            superblock.Lblock->block
                    .contract(theta,                                      idx({0},{1}))
                    .contract(superblock.HA->MPO,                        idx({1,2},{0,2}))
                    .contract(superblock.HB->MPO,                        idx({3,1},{0,2}))
                    .contract(theta.conjugate(),                          idx({0,2,4},{1,0,2}))
                    .contract(superblock.Rblock->block,                  idx({0,2,1},{0,1,2}));
    assert(abs(imag(E(0))) < 1e-10 and "Energy has an imaginary part!!!");
    if(sim_type == SimulationType::iDMRG ){
        //iDMRG uses reduced mpo in the environments!
        energy_mpo_all_sites = std::real(E(0))/2.0*L;
        energy_mpo           = std::real(E(0))/2.0;

    }else{
        energy_mpo_all_sites = std::real(E(0));
        energy_mpo           = std::real(energy_mpo_all_sites)/L;
    }

    t_ene_mpo.toc();

}

template<typename T>
double class_measurement::compute_energy_mpo(const class_superblock & superblock, const T * theta_ptr, Eigen::DSizes<long,4> dsizes) {
    Eigen::Tensor<std::complex<double>, 4> theta = Eigen::TensorMap<const Eigen::Tensor<const T, 4>>(theta_ptr, dsizes).template cast<std::complex<double>>();
//    t_ene_mpo.tic();
    Eigen::Tensor<Scalar, 0>  E =
                     superblock.Lblock->block
                    .contract(theta,                                      idx({0},{1}))
                    .contract(superblock.HA->MPO,                        idx({1,2},{0,2}))
                    .contract(superblock.HB->MPO,                        idx({3,1},{0,2}))
                    .contract(theta.conjugate(),                          idx({0,2,4},{1,0,2}))
                    .contract(superblock.Rblock->block,                  idx({0,2,1},{0,1,2}));
//    t_ene_mpo.toc();
    return std::real(E(0));
}

template double class_measurement::compute_energy_mpo(const class_superblock & superblock, const double *,Eigen::DSizes<long,4>);
template double class_measurement::compute_energy_mpo(const class_superblock & superblock, const std::complex<double> *, Eigen::DSizes<long,4>);




void class_measurement::compute_energy_ham(const class_superblock & superblock){
    t_ene_ham.tic();
    auto SX = qm::gen_manybody_spin(qm::spinOneHalf::sx,2);
    auto SY = qm::gen_manybody_spin(qm::spinOneHalf::sy,2);
    auto SZ = qm::gen_manybody_spin(qm::spinOneHalf::sz,2);
    auto h_evn = superblock.HA->single_site_hamiltonian(0,2,SX,SY, SZ);
    auto h_odd = superblock.HB->single_site_hamiltonian(1,2,SX,SY, SZ);
    MPS_Tools::Common::Views::compute_mps_components(superblock);
    using namespace MPS_Tools::Common::Views;

    E_evn = theta_evn_normalized
                .contract(Matrix_to_Tensor(h_evn,2,2,2,2),     idx({0, 2}, {0, 1}))
                .contract(theta_evn_normalized.conjugate(), idx({2, 3}, {0, 2}))
                .contract(l_evn,                            idx({0, 2}, {0, 1}))
                .contract(r_evn,                            idx({0, 1}, {0, 1}));
    E_odd  = theta_odd_normalized
                    .contract(Matrix_to_Tensor(h_odd,2,2,2,2),    idx({0, 2}, {0, 1}))
                    .contract(theta_odd_normalized.conjugate(),idx({2, 3}, {0,2}))
                    .contract(l_odd,                           idx({0,2}, {0,1}))
                    .contract(r_odd,                           idx({0,1},{0,1}));
    assert(abs(imag(E_evn(0)+ E_odd(0))) < 1e-10 and "Energy has an imaginary part!!!" );

    t_ene_ham.toc();
    energy_ham = 0.5*std::real(E_evn(0) + E_odd(0));
}


void class_measurement::compute_entanglement_entropy(const class_superblock & superblock){
    t_entropy.tic();
    Eigen::Tensor<Scalar,0> SA  = -superblock.MPS->LC.square()
            .contract(superblock.MPS->LC.square().log().eval(), idx({0},{0}));
    entanglement_entropy =  std::real(SA(0));
    t_entropy.toc();
}

//
double class_measurement::compute_parity(const class_finite_chain_state & state,const Eigen::Matrix2cd  paulimatrix){

    Eigen::TensorRef<Eigen::Tensor<Scalar,3>> temp;

    auto mpsL        = state.get_MPS_L().begin();
    auto endL        = state.get_MPS_L().end();
    auto [mpo,L,R]   = class_mpo::pauli_mpo(paulimatrix);

    int iter = 0;
    while(mpsL != endL){
        const Eigen::Tensor<Scalar,1> &LA = mpsL->get_L(); // std::get<0>(*mpsL);
        const Eigen::Tensor<Scalar,3> &GA = mpsL->get_G(); // std::get<1>(*mpsL);
        assert(LA.dimension(0) == L.dimension(0));
        assert(LA.dimension(0) == GA.dimension(1));

        temp = L.contract(asDiagonal(LA), idx({0},{0}))
                .contract(asDiagonal(LA), idx({0},{0}))
                .contract(mpo           ,idx({0},{0}))
                .contract(GA,                  idx({0,3},{1,0}))
                .contract(GA.conjugate(),      idx({0,2},{1,0}))
                .shuffle(array3{1,2,0});
        L = temp;
        mpsL++;
        iter++;
    }


    //Contract the center point
    auto &MPS_C = state.get_MPS_C();
    temp = L.contract(asDiagonal(MPS_C) , idx({0},{0}))
            .contract(asDiagonal(MPS_C) , idx({0},{0}))
            .shuffle(array3{1,2,0});
    L = temp;

    //Contract the right half of the chain
    auto mpsR  = state.get_MPS_R().begin();
    auto endR  = state.get_MPS_R().end();
    while(mpsR != endR){
        const Eigen::Tensor<Scalar,3> &GB = mpsR->get_G(); // std::get<0>(*mpsR);
        const Eigen::Tensor<Scalar,1> &LB = mpsR->get_L(); // std::get<1>(*mpsR);
        assert(GB.dimension(1) == L.dimension(0));
        assert(LB.dimension(0) == GB.dimension(2));
        temp = L.contract(GB,            idx({0},{1}))
                .contract(GB.conjugate(),idx({0},{1}))
                .contract(mpo           ,idx({0,1,3},{0,2,3}))
                .contract(asDiagonal(LB),idx({0},{0}))
                .contract(asDiagonal(LB),idx({0},{0}))
                .shuffle(array3{1,2,0});
        L = temp;
        mpsR++;
    }

    assert(L.dimensions() == R.dimensions());
    Eigen::Tensor<Scalar,0> parity_tmp = L.contract(R, idx({0,1,2},{0,1,2}));
    double parity = std::real(parity_tmp(0));
    std::cout << setprecision(16) << "Parity: " << parity_tmp << std::endl;
    return parity;
}


//void class_measurement::compute_parity_projected_mps(const Eigen::Matrix2cd  paulimatrix){
//
//    long A_length           = env_storage->get_MPS_L().size();
//    long B_length           = env_storage->get_MPS_R().size();
//
//    // Get handles for the current mps
//    const auto &mpsL        = env_storage->get_MPS_L();
//    const auto &mpsC        = env_storage->get_MPS_C();
//    const auto &mpsR        = env_storage->get_MPS_R();
//    const auto [mpo,L,R]   = class_mpo::pauli_mpo(paulimatrix);
//
//    // Make containers for the new updated mps
//    std::vector<Eigen::Tensor<Scalar,3>> A_matrices(A_length);
//    std::vector<Eigen::Tensor<Scalar,3>> B_matrices(B_length);
//    Eigen::Tensor<Scalar,1> C_matrix;
//
//    // Start updating  A matrices
//    int pos = 0;
//    long mpoDimL = mpo.dimension(0);
//    long mpoDimR = mpo.dimension(1);
//    for (auto &mps : mpsL){
//        assert(pos == mps.get_position() and "ERROR: Position and pos do not match");
//        long d    = mps.get_spin_dim();
//        long chiL = mps.get_chiL();
//        long chiR = mps.get_chiR();
//        A_matrices[pos++] =
//                mps.ref_A()
//                .contract(mpo, idx({0},{2}))
//                .shuffle(array5{4,2,0,3,1})
//                .reshape(array3{d, chiL*mpoDimL, chiR*mpoDimR});
//    }
//    // Increase the size of mpsC (the center bond) by doing an outer product
//    // with an identity matrix that has the same dimension as the MPO, contracting over no indices,
//    // mpsC       =  0 ------ 1
//    // I_matrix   =  0 ------ 1
//    // C_matrix   = 01 ====== 23
//    long mpsCdim = mpsC.dimension(0);
//    long mpodim = mpo.dimension(0);
//    auto I_tensor = Textra::Matrix_to_Tensor2(Eigen::MatrixXcd::Identity(mpodim,mpodim));
//    C_matrix = mpsC.contract(I_tensor, idx()).shuffle(array4{0,2,1,3}).reshape(array2{mpsCdim*mpodim, mpsCdim*mpodim});
//
//    // Start updating  B matrices
//    pos = 0;
//    for (auto &mps : mpsR){
//        assert(pos == mps.get_position() and "ERROR: Position and pos do not match");
//        long d    = mps.get_spin_dim();
//        long chiL = mps.get_chiL();
//        long chiR = mps.get_chiR();
//        B_matrices[pos++] =
//                mps.ref_B()
//                        .contract(mpo, idx({0},{2}))
//                        .shuffle(array5{4,2,0,3,1})
//                        .reshape(array3{d, chiL*mpoDimL, chiR*mpoDimR});
//    }
//
//    // Now we have updated versions of all matrices. We need to split them
//    // up into Gamma Lambda form, and because the updated matrices do not yield
//    // a normalized wave function, we need to do SVD from the edges of the chain
//    std::list<class_vidal_mps> mpsL_new;
//    std::list<class_vidal_mps> mpsR_new;
//    Eigen::Tensor<Scalar,1>    mpsC_new;
//
//    class_SVD<Scalar> svd;
////    svd.setThreshold(1e-8);
//    Eigen::Tensor<Scalar,3> A_next;
//    Eigen::Tensor<Scalar,3> B_next;
//    Eigen::Tensor<Scalar,1> LambdaA_next(1);
//    Eigen::Tensor<Scalar,1> LambdaB_next(1);
//    LambdaA_next.setConstant(1);
//    LambdaB_next.setConstant(1);
//    for (size_t pos = 0; pos < A_matrices.size()-1; pos++){
//        long d    = A_matrices[pos].dimension(0);
//        long chiL = A_matrices[pos].dimension(1);
//        long chiR = A_matrices[pos].dimension(2);
//        long rows = d*chiL;
//        long cols = chiR;
//        auto [U,S,V] = svd.decompose(A_matrices[pos],rows,cols);
//        A_next = Textra::asDiagonalInversed(LambdaA_next).contract(U.reshape(Textra::array3{chiL,d,chiR}), idx({1},{1})).shuffle(array3{1,0,2});
//        mpsL_new.emplace_back(A_next,LambdaA_next);
//
//        A_next = Textra::asDiagonalInversed(S)
//                .contract(V, idx({1},{0}))
//                .contract(A_matrices[pos+1], idx({1},{1}))
//                .shuffle(array3{1,0,2});
//
//        A_matrices[pos+1] = A_next;
//        LambdaA_next = S;
//    }
//
//    for (size_t pos = B_matrices.size(); pos > 0; pos--){
//        long d    = B_matrices[pos].dimension(0);
//        long chiL = B_matrices[pos].dimension(1);
//        long chiR = B_matrices[pos].dimension(2);
//        long rows = chiL;
//        long cols = d*chiR;
//        auto [U,S,V] = svd.decompose(B_matrices[pos],rows,cols);
//        B_next = V.reshape(Textra::array3{chiL,d,chiR}).contract(Textra::asDiagonalInversed(LambdaB_next), idx({2},{0}));
//        mpsL_new.emplace_front(B_next,LambdaB_next);
//        B_next = B_matrices[pos-1]
//                .contract(U, idx({2},{0}))
//                .contract(Textra::asDiagonalInversed(S),idx({2},{0}));
//        B_matrices[pos-1] = B_next;
//        LambdaB_next = S;
//    }
//
//    // We have now worked our way from the edges of the chain in to the middle,
//    // so we have something of the form
//    // L G L G A_{i} C B_{i+1} G L G L
//    //
//    // The middle chunk A C B needs a 2-site SVD, and then
//    // G_i     = LambdaA_next^-1  U
//    // C       = S
//    // G_i+1   = V LambdaB_next^-1
//
//    Eigen::Tensor<Scalar,4> ACB = A_matrices.back()
//            .contract(C_matrix, idx({2},{0}))
//            .contract(B_matrices.front(), idx({2},{1}));
//    auto [U,S,V] = svd.schmidt(ACB);
//    mpsL_new.emplace_back (Textra::asDiagonalInversed(LambdaA_next).contract(U, idx({1},{1})).shuffle(array3{1,0,2}));
//    mpsR_new.emplace_front(V.contract(Textra::asDiagonalInversed(LambdaB_next), idx({2},{0})));
//    mpsC_new = S;
//
//    // Now update the original MPS in env_storage
//    env_storage->replace_mps(mpsL_new, mpsC_new, mpsR_new);
//
//
//
//
//
////    assert(L.dimensions() == R.dimensions());
////    Eigen::Tensor<Scalar,0> parity_tmp = L.contract(R, idx({0,1,2},{0,1,2}));
////    double parity = std::real(parity_tmp(0));
////    std::cout << setprecision(16) << "Parity: " << parity_tmp << std::endl;
//}




void class_measurement::compute_finite_chain_energy(const class_finite_chain_state & state){
    Eigen::Tensor<Scalar,3> L = state.get_ENV_L().front().block;
    auto mpsL  = state.get_MPS_L().begin();
    auto mpoL  = state.get_MPO_L().begin();
    auto endL  = state.get_MPS_L().end();
    while(mpsL != endL){
        const Eigen::Tensor<Scalar,3> A = mpsL->get_A();
        Eigen::Tensor<Scalar,3> temp =
                L
                .contract(A                   , idx({0},{1}))
                .contract(mpoL->get()->MPO    , idx({1,2},{0,2}))
                .contract(A.conjugate()       , idx({0,3},{1,0}))
                .shuffle(array3{0,2,1});

        L = temp;
        mpsL++;
        mpoL++;
    }

    {
        //Contract the center point
        auto &MPS_C = state.get_MPS_C();
        Eigen::Tensor<Scalar,3> temp =
                L
                .contract(asDiagonal(MPS_C), idx({0}, {0}))
                .contract(asDiagonal(MPS_C), idx({0}, {0}))
                .shuffle(array3{1, 2, 0});
        L = temp;
    }
    //Contract the right half of the chain
    Eigen::Tensor<Scalar,3> R = state.get_ENV_R().back().block;
    auto mpsR  = state.get_MPS_R().begin();
    auto mpoR  = state.get_MPO_R().begin();
    auto endR  = state.get_MPS_R().end();
    while(mpsR != endR){
        const Eigen::Tensor<Scalar,3> B = mpsR->get_B();
        Eigen::Tensor<Scalar,3> temp =
                L
                .contract(B                   , idx({0},{1}))
                .contract(mpoR->get()->MPO    , idx({1,2},{0,2}))
                .contract(B.conjugate()       , idx({0,3},{1,0}))
                .shuffle(array3{0,2,1});
        L = temp;
        mpsR++;
        mpoR++;
    }

    assert(L.dimensions() == R.dimensions());
    Eigen::Tensor<Scalar,0> E_all_sites = L.contract(R, idx({0,1,2},{0,1,2}));
    energy_chain = std::real(E_all_sites(0));
    std::cout << setprecision(16) << "E all sites: " << energy_chain << " | per site: " << energy_chain / state.get_length() << std::endl;
}


void class_measurement::compute_finite_chain_norm(const class_finite_chain_state & state){

    auto mpsL  = state.get_MPS_L().begin();
    auto endL  = state.get_MPS_L().end();
    const Eigen::Tensor<Scalar,1>  & LA_ = mpsL->get_L(); // std::get<0>(*mpsL);
    const Eigen::Tensor<Scalar,3>  & GA_ = mpsL->get_G(); // std::get<1>(*mpsL);
    Eigen::Tensor<Scalar,2> chain =
             asDiagonal(LA_)
            .contract(asDiagonal(LA_), idx({0},{1}))
            .contract(GA_                 , idx({1},{1}))
            .contract(GA_.conjugate()     , idx({0,1},{1,0}));
    Eigen::TensorRef<Eigen::Tensor<Scalar,2>> temp;
    mpsL++;
    while(mpsL != endL){
        const Eigen::Tensor<Scalar,1>  & LA = mpsL->get_L() ; // std::get<0>(*mpsL);
        const Eigen::Tensor<Scalar,3>  & GA = mpsL->get_G() ; // std::get<1>(*mpsL);
        assert(LA.dimension(0) == chain.dimension(0));
        assert(LA.dimension(0) == GA.dimension(1));
        temp = chain
                .contract(asDiagonal(LA), idx({0},{0}))
                .contract(asDiagonal(LA), idx({0},{0}))
                .contract(GA            , idx({0},{1}))
                .contract(GA.conjugate(), idx({0,1},{1,0}));
        chain = temp;
        mpsL++;
    }

//    Contract the center point
    auto &MPS_C = state.get_MPS_C();
    temp = chain
            .contract(asDiagonal(MPS_C), idx({0},{0}))
            .contract(asDiagonal(MPS_C), idx({0},{0}));
    chain = temp;

    //Contract the right half of the chain
    auto mpsR  = state.get_MPS_R().begin();
    auto endR  = state.get_MPS_R().end();

    while(mpsR != endR){
        const Eigen::Tensor<Scalar,3> &GB  = mpsR->get_G(); //std::get<0>(*mpsR);
        const Eigen::Tensor<Scalar,1> &LB  = mpsR->get_L(); //std::get<1>(*mpsR);
        assert(GB.dimension(1) == chain.dimension(0));
        assert(LB.dimension(0) == GB.dimension(2));

        temp = chain
                .contract(GB              , idx({0},{1}))
                .contract(GB.conjugate()  , idx({0,1},{1,0}))
                .contract(asDiagonal(LB)  , idx({0},{0}))
                .contract(asDiagonal(LB)  , idx({0},{0}));
        chain = temp;
        mpsR++;
    }
    norm_chain = std::abs(Textra::Tensor2_to_Matrix(chain).trace());
    std::cout << setprecision(16) << "Norm: " << norm_chain << std::endl;
    assert(std::abs(norm_chain - 1.0) < 1e-10 );
}

void class_measurement::compute_finite_chain_norm2(const class_finite_chain_state & state){

    auto mpsL  = state.get_MPS_L().begin();
    auto endL  = state.get_MPS_L().end();
    const Eigen::Tensor<Scalar,3>  & A = mpsL->ref_A(); // std::get<1>(*mpsL);
    Eigen::Tensor<Scalar,2> chain = A.contract(A.conjugate(), idx({0,1},{0,1}));
    Eigen::TensorRef<Eigen::Tensor<Scalar,2>> temp;
    mpsL++;
    while(mpsL != endL){
        const Eigen::Tensor<Scalar,3>  & A = mpsL->ref_A() ; // std::get<1>(*mpsL);
        assert(A.dimension(1) == chain.dimension(0));
        temp = chain
                .contract(A,             idx({0},{1}))
                .contract(A.conjugate(), idx({0,1},{1,0}));
        chain = temp;
        mpsL++;
    }

//    Contract the center point
    auto &MPS_C = state.get_MPS_C();
    temp = chain
            .contract(asDiagonal(MPS_C), idx({0},{0}))
            .contract(asDiagonal(MPS_C), idx({0},{0}));
    chain = temp;
    //Contract the right half of the chain
    auto mpsR  = state.get_MPS_R().begin();
    auto endR  = state.get_MPS_R().end();

    while(mpsR != endR){
        const Eigen::Tensor<Scalar,3> &B  = mpsR->ref_B(); //std::get<0>(*mpsR);
        assert(B.dimension(1) == chain.dimension(0));
        Eigen::Tensor<Scalar,3>  tempB = chain.contract(B, idx({0},{1}));
        temp = chain
                .contract(B               , idx({0},{1}))
                .contract(B.conjugate()   , idx({0,1},{1,0}));
        chain = temp;
        mpsR++;
    }
    norm_chain = std::abs(Textra::Tensor2_to_Matrix(chain).trace());
    std::cout << setprecision(16) << "Norm: " << norm_chain << std::endl;
    assert(std::abs(norm_chain - 1.0) < 1e-10 );
}

void class_measurement::compute_finite_chain_mps_wavefn(const class_finite_chain_state & state){

    auto mpsL  = state.get_MPS_L().begin();
    auto endL  = state.get_MPS_L().end();
    Eigen::Tensor<Scalar,2> chain(1,1);
    chain.setConstant(1.0);
    Eigen::TensorRef<Eigen::Tensor<Scalar,2>> temp;
    // The "chain" is a matrix whose 0 index keeps growing.
    // For each site that passes, it grows by GA.dimension(0) = phys dim
    // Say the chain is a 16x7 matrix (having contracted 4 particles, and the latest
    // chi was 7). Then contracting the next site, with dimensions 2x7x9 will get you a
    // 16x2x9 tensor. Now the reshaping convert it into a 32 x 9 matrix. Because
    // Eigen is column major, the doubling 16->32 will stack the third index twice.

    while(mpsL != endL){
        const Eigen::Tensor<Scalar,1>  & LA = mpsL->get_L(); //std::get<0>(*mpsL);
        const Eigen::Tensor<Scalar,3>  & GA = mpsL->get_G(); //std::get<1>(*mpsL);
        assert(LA.dimension(0) == GA.dimension(1));

        temp = chain
                .contract(asDiagonal(LA), idx({1},{0}))
                .contract(GA            , idx({1},{1}))
                .reshape(array2{GA.dimension(0) * chain.dimension(0), GA.dimension(2)});
        chain = temp;
        mpsL++;
    }

//    Contract the center point
    auto &MPS_C = state.get_MPS_C();
    temp = chain.contract(asDiagonal(MPS_C), idx({1},{0}));
    chain = temp;

    //Contract the right half of the chain
    auto mpsR  = state.get_MPS_R().begin();
    auto endR  = state.get_MPS_R().end();

    while(mpsR != endR){
        const Eigen::Tensor<Scalar,3> &GB  = mpsR->get_G(); // std::get<0>(*mpsR);
        const Eigen::Tensor<Scalar,1> &LB  = mpsR->get_L(); // std::get<1>(*mpsR);
        assert(LB.dimension(0) == GB.dimension(2));
        temp = chain
                .contract(GB              , idx({1},{1}))
                .contract(asDiagonal(LB)  , idx({2},{0}))
                 .reshape(array2{GB.dimension(0) * chain.dimension(0), GB.dimension(2)});
        chain = temp;
        mpsR++;
    }
    mps_chain = chain.reshape(array1{chain.dimension(0)});
    norm_chain = Textra::Tensor2_to_Matrix(chain).norm();
    assert(std::abs(norm_chain - 1.0) < 1e-10 );
}


double class_measurement::get_energy_mpo(){return energy_mpo;}
double class_measurement::get_energy_ham(){return energy_ham;}
double class_measurement::get_energy_mom(){return energy_mom;}

void   class_measurement::set_variance(double new_variance){
    variance_mpo = new_variance;
    variance_ham = new_variance;
    variance_mom = new_variance;
}

double class_measurement::get_variance_mpo(){return std::abs(variance_mpo);}
double class_measurement::get_variance_ham(){return std::abs(variance_ham);}
double class_measurement::get_variance_mom(){return std::abs(variance_mom);}

double class_measurement::get_entanglement_entropy(const class_superblock & superblock){
    compute_entanglement_entropy(superblock);
    return entanglement_entropy;
}

double class_measurement::get_truncation_error(const class_superblock & superblock){
    return superblock.MPS->truncation_error;
}

std::vector<double> class_measurement::get_parity(){
    return {parity_x,parity_y,parity_z};
}
double class_measurement::get_parity_x(){return parity_x;}
double class_measurement::get_parity_y(){return parity_y;}
double class_measurement::get_parity_z(){return parity_z;}

long class_measurement::get_chi(const class_superblock & superblock){
    return superblock.MPS->chiC();
}

long class_measurement::get_chi(const class_finite_chain_state & state){
    return state.get_MPS_C().size();
}

long class_measurement::get_chain_length(const class_superblock & superblock){
    return superblock.environment_size + 2;
}
long class_measurement::get_chain_length(const class_finite_chain_state & state){
    return state.get_length();
}


// Profiling

void class_measurement::set_profiling_labels() {
    using namespace settings::profiling;
    t_ene_mpo.set_properties(on, precision,"↳ Energy (MPO)           ");
    t_ene_ham.set_properties(on, precision,"↳ Energy (HAM)           ");
    t_ene_gen.set_properties(on, precision,"↳ Energy (MOM)           ");
    t_var_mpo.set_properties(on, precision,"↳ Variance (MPO)         ");
    t_var_ham.set_properties(on, precision,"↳ Variance (HAM)         ");
    t_var_gen.set_properties(on, precision,"↳ Variance (MOM)         ");
    t_entropy.set_properties(on, precision,"↳ Ent. Entropy           ");
    t_temp1.set_properties(on, precision,  "↳ Temp1                  ");
    t_temp2.set_properties(on, precision,  "↳ Temp2                  ");
    t_temp3.set_properties(on, precision,  "↳ Temp3                  ");
    t_temp4.set_properties(on, precision,  "↳ Temp4                  ");

}

void class_measurement::print_profiling(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        std::cout << "\nComputing observables breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_ene_mpo.print_time_w_percent(t_parent);
        t_ene_ham.print_time_w_percent(t_parent);
        t_ene_gen.print_time_w_percent(t_parent);
        t_var_mpo.print_time_w_percent(t_parent);
        t_var_ham.print_time_w_percent(t_parent);
        t_var_gen.print_time_w_percent(t_parent);
        t_entropy.print_time_w_percent(t_parent);
        t_temp1.print_time_w_percent(t_parent);
        t_temp2.print_time_w_percent(t_parent);
        t_temp3.print_time_w_percent(t_parent);
        t_temp4.print_time_w_percent(t_parent);
    }
}

