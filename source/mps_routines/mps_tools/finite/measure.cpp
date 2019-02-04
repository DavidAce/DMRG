//
// Created by david on 2019-02-01.
//

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
using Scalar = std::complex<double>;



namespace MPS_Tools::Finite::Measure::Results {
    size_t length;
    size_t bond_dimension                       = 0;
    double norm                                 = 0;
    double truncation_error                     = 0;
    double energy_mpo                           = 0;
    double energy_per_site_mpo                  = 0;
    double energy_variance_mpo                  = 0;
    double midchain_entanglement_entropy        = 0;
    std::vector<double> entanglement_entropy    = std::vector<double>();
    double parity_sx                            = 0;
    double parity_sy                            = 0;
    double parity_sz                            = 0;
    Eigen::Tensor<Scalar,1> mps_wavefn          = Eigen::Tensor<Scalar,1>();
    bool state_measured                         = false;
}






void MPS_Tools::Finite::Measure::do_all_measurements(const class_finite_chain_state & state){
    using namespace MPS_Tools::Finite::Measure;
    if (Results::state_measured){return;}

    Results::length                         = length(state);
    Results::bond_dimension                 = bond_dimension(state);
    Results::norm                           = norm(state);
    Results::energy_mpo                     = energy_mpo(state);  //This number is needed for variance calculation!
    Results::energy_per_site_mpo            = Results::energy_mpo/Results::length;
    Results::energy_variance_mpo            = energy_variance_mpo(state, Results::energy_mpo);
    Results::midchain_entanglement_entropy  = midchain_entanglement_entropy(state);
    Results::entanglement_entropy           = entanglement_entropy(state);
    Results::parity_sx                      = parity(state,qm::spinOneHalf::sx);
    Results::parity_sy                      = parity(state,qm::spinOneHalf::sy);
    Results::parity_sz                      = parity(state,qm::spinOneHalf::sz);

    if (Results::length <= 14){
        Results::mps_wavefn = mps_wavefn(state);
    }

    Results::state_measured = true;
}


size_t MPS_Tools::Finite::Measure::length(const class_finite_chain_state & state){
    return state.get_length();
}

double MPS_Tools::Finite::Measure::norm(const class_finite_chain_state & state){
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
    double norm_chain = std::real(Textra::Tensor2_to_Matrix(chain).trace());
    std::cout << setprecision(16) << "Norm: " << norm_chain << std::endl;
    assert(std::abs(norm_chain - 1.0) < 1e-10 );
    return norm_chain;
}


size_t MPS_Tools::Finite::Measure::bond_dimension(const class_finite_chain_state & state){
    return state.get_MPS_C().dimension(0);
}


double MPS_Tools::Finite::Measure::energy_mpo(const class_finite_chain_state & state){
    auto theta = MPS_Tools::Common::Views::get_theta(state);
    Eigen::Tensor<Scalar, 0>  E =
            state.get_ENV_L().back().block
                    .contract(theta,                                     idx({0},{1}))
                    .contract(state.get_MPO_L().back()->MPO,             idx({1,2},{0,2}))
                    .contract(state.get_MPO_R().front()->MPO,            idx({3,1},{0,2}))
                    .contract(theta.conjugate(),                         idx({0,2,4},{1,0,2}))
                    .contract(state.get_ENV_R().front().block,           idx({0,2,1},{0,1,2}));
    assert(abs(imag(E(0))) < 1e-10 and "Energy has an imaginary part!!!");

    return std::real(E(0)) ;
}



double MPS_Tools::Finite::Measure::energy_variance_mpo(const class_finite_chain_state & state){
//    t_var_mpo.tic();
    double energy = MPS_Tools::Finite::Measure::energy_mpo(state);
    auto theta = MPS_Tools::Common::Views::get_theta(state);
    Eigen::Tensor<Scalar, 0> H2 =
            state.get_ENV2_L().back().block
                    .contract(theta                              , idx({0}  ,{1}))
                    .contract(state.get_MPO_L().back()->MPO      , idx({1,3},{0,2}))
                    .contract(state.get_MPO_R().front()->MPO     , idx({4,2},{0,2}))
                    .contract(state.get_MPO_L().back()->MPO      , idx({1,3},{0,2}))
                    .contract(state.get_MPO_R().front()->MPO     , idx({4,3},{0,2}))
                    .contract(theta.conjugate()                  , idx({0,3,5},{1,0,2}))
                    .contract(state.get_ENV2_R().front().block   , idx({0,3,1,2},{0,1,2,3}));

    return std::abs(H2(0) - energy*energy);

//
//    if (sim_type == SimulationType::iDMRG) {
//        t_var_mpo.toc();
//        variance_mpo_all_sites = std::real(H2(0))/2.0;
//        variance_mpo           = std::real(H2(0))/2.0/L;
//    }else{
//        variance_mpo_all_sites = std::abs(H2(0) - energy_mpo_all_sites*energy_mpo_all_sites);
//        variance_mpo           = variance_mpo_all_sites/L;
//    }
//    t_var_mpo.toc();
}


double MPS_Tools::Finite::Measure::energy_variance_mpo(const class_finite_chain_state & state, double energy_mpo){
//    t_var_mpo.tic();
    auto theta = MPS_Tools::Common::Views::get_theta(state);
    Eigen::Tensor<Scalar, 0> H2 =
            state.get_ENV2_L().back().block
                    .contract(theta                              , idx({0}  ,{1}))
                    .contract(state.get_MPO_L().back()->MPO      , idx({1,3},{0,2}))
                    .contract(state.get_MPO_R().front()->MPO     , idx({4,2},{0,2}))
                    .contract(state.get_MPO_L().back()->MPO      , idx({1,3},{0,2}))
                    .contract(state.get_MPO_R().front()->MPO     , idx({4,3},{0,2}))
                    .contract(theta.conjugate()                  , idx({0,3,5},{1,0,2}))
                    .contract(state.get_ENV2_R().front().block   , idx({0,3,1,2},{0,1,2,3}));

    return std::abs(H2(0) - energy_mpo*energy_mpo);

//
//    if (sim_type == SimulationType::iDMRG) {
//        t_var_mpo.toc();
//        variance_mpo_all_sites = std::real(H2(0))/2.0;
//        variance_mpo           = std::real(H2(0))/2.0/L;
//    }else{
//        variance_mpo_all_sites = std::abs(H2(0) - energy_mpo_all_sites*energy_mpo_all_sites);
//        variance_mpo           = variance_mpo_all_sites/L;
//    }
//    t_var_mpo.toc();
}







double MPS_Tools::Finite::Measure::midchain_entanglement_entropy(const class_finite_chain_state & state){
//    t_entropy.tic();
    auto & LC = state.get_MPS_C();

    Eigen::Tensor<Scalar,0> SA  = -LC.square()
            .contract(LC.square().log().eval(), idx({0},{0}));
    return std::real(SA(0));
//    t_entropy.toc();
}


std::vector<double> MPS_Tools::Finite::Measure::entanglement_entropy(const class_finite_chain_state & state){
//    t_entropy.tic();
    std::vector<double> SA;
    for (auto & mps : state.get_MPS_L()) {
        auto &L = mps.get_L();
        Eigen::Tensor<Scalar, 0> SA_L = -L.square().contract(L.square().log().eval(), idx({0}, {0}));
        SA.emplace_back(std::real(SA_L(0)));
    }
    SA.emplace_back(midchain_entanglement_entropy(state));
    for (auto & mps : state.get_MPS_R()) {
        auto &L = mps.get_L();
        Eigen::Tensor<Scalar, 0> SA_R = -L.square().contract(L.square().log().eval(), idx({0}, {0}));
        SA.emplace_back(std::real(SA_R(0)));
    }
    return SA;
//    t_entropy.toc();
}



//
double MPS_Tools::Finite::Measure::parity(const class_finite_chain_state & state,const Eigen::Matrix2cd  paulimatrix){

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


Eigen::Tensor<Scalar,1> MPS_Tools::Finite::Measure::mps_wavefn(const class_finite_chain_state & state){

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
    Eigen::Tensor<Scalar,1> mps_chain = chain.reshape(array1{chain.dimension(0)});
    double norm_chain = Textra::Tensor2_to_Matrix(chain).norm();
    assert(std::abs(norm_chain - 1.0) < 1e-10 );
    return mps_chain;
}

//
//double class_measurement::get_energy_mpo(){return energy_mpo;}
//double class_measurement::get_energy_ham(){return energy_ham;}
//double class_measurement::get_energy_mom(){return energy_mom;}
//
//void   class_measurement::set_variance(double new_variance){
//    variance_mpo = new_variance;
//    variance_ham = new_variance;
//    variance_mom = new_variance;
//}
//
//double class_measurement::get_variance_mpo(){return std::abs(variance_mpo);}
//double class_measurement::get_variance_ham(){return std::abs(variance_ham);}
//double class_measurement::get_variance_mom(){return std::abs(variance_mom);}
//
//double class_measurement::get_entanglement_entropy(const class_superblock & superblock){
//    compute_entanglement_entropy(superblock);
//    return entanglement_entropy;
//}
//
//double class_measurement::get_truncation_error(const class_superblock & superblock){
//    return superblock.MPS->truncation_error;
//}


//// Profiling
//
//void class_measurement::set_profiling_labels() {
//    using namespace settings::profiling;
//    t_ene_mpo.set_properties(on, precision,"↳ Energy (MPO)           ");
//    t_ene_ham.set_properties(on, precision,"↳ Energy (Ham)           ");
//    t_ene_gen.set_properties(on, precision,"↳ Energy (Gen)           ");
//    t_var_mpo.set_properties(on, precision,"↳ Variance (MPO)         ");
//    t_var_ham.set_properties(on, precision,"↳ Variance (Ham)         ");
//    t_var_gen.set_properties(on, precision,"↳ Variance (Gen)         ");
//    t_entropy.set_properties(on, precision,"↳ Ent. Entropy           ");
//    t_temp1.set_properties(on, precision,  "↳ Temp1                  ");
//    t_temp2.set_properties(on, precision,  "↳ Temp2                  ");
//    t_temp3.set_properties(on, precision,  "↳ Temp3                  ");
//    t_temp4.set_properties(on, precision,  "↳ Temp4                  ");
//
//}
//
//void class_measurement::print_profiling(class_tic_toc &t_parent){
//    if (settings::profiling::on) {
//        std::cout << "\nComputing observables breakdown:" << std::endl;
//        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
//        t_ene_mpo.print_time_w_percent(t_parent);
//        t_ene_ham.print_time_w_percent(t_parent);
//        t_ene_gen.print_time_w_percent(t_parent);
//        t_var_mpo.print_time_w_percent(t_parent);
//        t_var_ham.print_time_w_percent(t_parent);
//        t_var_gen.print_time_w_percent(t_parent);
//        t_entropy.print_time_w_percent(t_parent);
//        t_temp1.print_time_w_percent(t_parent);
//        t_temp2.print_time_w_percent(t_parent);
//        t_temp3.print_time_w_percent(t_parent);
//        t_temp4.print_time_w_percent(t_parent);
//    }
//}
//
