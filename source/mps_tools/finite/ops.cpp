//
// Created by david on 2019-01-30.
//


#include <mps_tools/nmspc_mps_tools.h>
#include <mps_state/class_finite_chain_state.h>
#include <mps_state/class_superblock.h>
#include <mps_state/class_mpo.h>
#include <model/class_hamiltonian_base.h>
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
#include <iomanip>
#include <spdlog/spdlog.h>
#include <general/nmspc_random_numbers.h>


using Scalar         = std::complex<double>;
using namespace Textra;

std::list<Eigen::Tensor<Scalar,4>> mpstools::finite::ops::make_mpo_list (const std::list<std::shared_ptr<class_hamiltonian_base>> &mpos_L, const std::list<std::shared_ptr<class_hamiltonian_base>> &mpos_R){
    std::list<Eigen::Tensor<Scalar,4>> mpos;
    for(auto &mpo_L : mpos_L){
        mpos.push_back(mpo_L->MPO());
    }
    for(auto &mpo_R : mpos_R){
        mpos.push_back(mpo_R->MPO());
    }
    return mpos;
}


void mpstools::finite::ops::apply_mpo(class_finite_chain_state &state,const Eigen::Tensor<Scalar,4> mpo,const Eigen::Tensor<Scalar,3> Ledge, const Eigen::Tensor<Scalar,3> Redge){
    std::list<Eigen::Tensor<Scalar,4>> mpos(state.get_length(), mpo);
    apply_mpos(state,mpos,Ledge,Redge);
}



void mpstools::finite::ops::apply_mpos(class_finite_chain_state &state,const std::list<Eigen::Tensor<Scalar,4>> &mpos,const Eigen::Tensor<Scalar,3> Ledge, const Eigen::Tensor<Scalar,3> Redge){
    if (mpos.size() != state.get_length()) throw std::runtime_error("Number of mpo's doesn't match the number of sites on the system");
    // Apply MPO's on Gamma matrices and
    // increase the size on all Lambdas by chi*mpoDim
    mpstools::log->trace("Applying MPO's");
//    state.unset_measurements();
//    std::cout << "Norm              (before mpos): " << mpstools::finite::measure::norm(state)  << std::endl;
//    std::cout << "Spin component sx (before mpos): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sx)  << std::endl;
//    std::cout << "Spin component sy (before mpos): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sy)  << std::endl;
//    std::cout << "Spin component sz (before mpos): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sz)  << std::endl;
    auto mpo = mpos.begin();

    {
        for (auto & mps : state.MPS_L){
            long mpoDimL = mpo->dimension(0);
            long mpoDimR = mpo->dimension(1);
            auto [d,chiL,chiR] = mps.get_dims();
            Eigen::Tensor<Scalar,3> G_temp =
                    mps.get_G()
                            .contract(*mpo, idx({0},{2}))
                            .shuffle(array5{4,0,2,1,3})
                            .reshape(array3{d, chiL*mpoDimL, chiR*mpoDimR});

            Eigen::Tensor<Scalar,1> L_temp = mps.get_L().broadcast(array1{mpoDimL});
            mps.set_G(G_temp);
            mps.set_L(L_temp);
            mpo++;
        }
    }

    //Take care of the middle
    {
        long mpoDim = mpo->dimension(0);
        Eigen::Tensor<Scalar, 1> C_temp = state.MPS_C.broadcast(array1{mpoDim});
        state.MPS_C = C_temp;
    }
    {
        for (auto & mps : state.MPS_R){
            long mpoDimL = mpo->dimension(0);
            long mpoDimR = mpo->dimension(1);
            auto [d,chiL,chiR] = mps.get_dims();
            Eigen::Tensor<Scalar,3> G_temp =
                    mps.get_G()
                            .contract(*mpo, idx({0},{2}))
                            .shuffle(array5{4,0,2,1,3})
                            .reshape(array3{d, chiL*mpoDimL, chiR*mpoDimR});

            Eigen::Tensor<Scalar,1> L_temp = mps.get_L().broadcast(array1{mpoDimR});
            mps.set_G(G_temp);
            mps.set_L(L_temp);
            mpo++;
        }
    }

    // Take care of the edges. Apply the left and right MPO-edges on A's and B's
    // so the left and right most lambdas become ones
    {
        long mpoDimL = mpos.front().dimension(0);
        auto Ldim    = Ledge.dimension(0);
        Eigen::Tensor<Scalar, 3> G_temp =
                Ledge
                        .shuffle(Textra::array3{0,2,1})
                        .reshape(Textra::array2{Ldim*mpoDimL,Ldim})
                        .contract(state.MPS_L.front().get_A(),Textra::idx({0},{1}))
                        .shuffle(Textra::array3{1,0,2});
        state.MPS_L.front().set_G(G_temp);
        state.MPS_L.front().set_L(Eigen::Tensor<Scalar,1>(Ldim).constant(1.0));
    }
    {
        long mpoDimR = mpos.back().dimension(1);
        auto Rdim    = Redge.dimension(0);
        Eigen::Tensor<Scalar, 3> G_temp =
                Redge
                        .shuffle(Textra::array3{0,2,1})
                        .reshape(Textra::array2{Rdim*mpoDimR,Rdim})
                        .contract(state.MPS_R.back().get_B(),Textra::idx({0},{2}))
                        .shuffle(Textra::array3{1,2,0});
        state.MPS_R.back().set_G(G_temp);
        state.MPS_R.back().set_L(Eigen::Tensor<Scalar,1>(Rdim).constant(1.0));
    }
    state.unset_measurements();
//    std::cout << "Norm              (after mpos): " << mpstools::finite::measure::norm(state)  << std::endl;
//    std::cout << "Spin component sx (after mpos): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sx)  << std::endl;
//    std::cout << "Spin component sy (after mpos): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sy)  << std::endl;
//    std::cout << "Spin component sz (after mpos): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sz)  << std::endl;
}


class_finite_chain_state mpstools::finite::ops::get_parity_projected_state(const class_finite_chain_state &state, const Eigen::MatrixXcd  paulimatrix, const int sign) {
    if (std::abs(sign) != 1) throw std::runtime_error("Expected 'sign' +1 or -1. Got: " + std::to_string(sign));
    mpstools::log->trace("Generating parity projected state");
    class_finite_chain_state state_projected = state;
    state_projected.unset_measurements();

    const auto [mpo,L,R]    = class_mpo::parity_projector_mpos(paulimatrix,state_projected.get_length(), sign);
    apply_mpos(state_projected,mpo, L,R);
    mpstools::finite::mps::normalize(state_projected);
    mpstools::finite::mps::rebuild_environments(state_projected);
    mpstools::finite::debug::check_integrity_of_mps(state_projected);
    return state_projected;
}

class_finite_chain_state mpstools::finite::ops::get_closest_parity_state(const class_finite_chain_state &state, const Eigen::MatrixXcd  paulimatrix) {
    mpstools::log->trace("Finding closest projection");
    double measured_spin_component = mpstools::finite::measure::spin_component(state, paulimatrix);
    if (measured_spin_component > 0){
        return get_parity_projected_state(state, paulimatrix, 1);
    }else{
        return get_parity_projected_state(state, paulimatrix,-1);
    }
}

class_finite_chain_state mpstools::finite::ops::get_closest_parity_state(const class_finite_chain_state &state, const std::string paulistring) {
    mpstools::log->trace("Finding closest projection");
    if      (paulistring == "sx"){return get_closest_parity_state(state,qm::spinOneHalf::sx);}
    else if (paulistring == "sy"){return get_closest_parity_state(state,qm::spinOneHalf::sy);}
    else if (paulistring == "sz"){return get_closest_parity_state(state,qm::spinOneHalf::sz);}
    else{
        mpstools::log->warn(R"(Wrong pauli string. Expected one of  "sx","sy" or "sz". Got: )" + paulistring);
        auto spin_components = mpstools::finite::measure::spin_components(state);
        auto max_idx = std::distance(spin_components.begin(), std::max_element(spin_components.begin(),spin_components.end()));
        if(max_idx == 0)      {return get_closest_parity_state(state,"sx"); }
        else if(max_idx == 1) {return get_closest_parity_state(state,"sy"); }
        else if(max_idx == 2) {return get_closest_parity_state(state,"sz"); }
        else {throw std::runtime_error("Wrong pauli string and could not find closest parity state");}
    }
}


template<typename Scalar, auto rank>
void throw_if_has_imaginary_part(const Eigen::Tensor<Scalar,rank> &tensor, double threshold = 1e-14) {
    Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>> vector (tensor.data(),tensor.size());
    if constexpr (std::is_same<Scalar, std::complex<double>>::value){
        auto imagSum = vector.imag().cwiseAbs().sum();
        if (imagSum > threshold){
            std::cout << vector << std::endl;
            throw std::runtime_error("Has imaginary part. Sum: " + std::to_string(imagSum));
        }
    }
}


double mpstools::finite::ops::overlap(const class_finite_chain_state &state1, const class_finite_chain_state &state2){

    assert(state1.get_length() == state2.get_length() and "ERROR: States have different lengths! Can't do overlap.");
    assert(state1.get_position() == state2.get_position() and "ERROR: States need to be at the same position! Can't do overlap.");

    Eigen::Tensor<Scalar,2> overlap =
            state1.MPS_L.front().get_A()
            .contract(state2.MPS_L.front().get_A().conjugate(), Textra::idx({0,1},{0,1}));
    auto mps1_it_L = state1.MPS_L.begin();
    auto mps2_it_L = state2.MPS_L.begin();
    mps1_it_L++;
    mps2_it_L++;

    while (mps1_it_L != state1.MPS_L.end() ){
        Eigen::Tensor<Scalar,2> temp = overlap
                .contract(mps1_it_L->get_A()            , Textra::idx({0},{1}))
                .contract(mps2_it_L->get_A().conjugate(), Textra::idx({0,1},{1,0}));
        overlap = temp;
        mps1_it_L++;
        mps2_it_L++;
    }

    Eigen::Tensor<Scalar,2> temp = overlap
            .contract(Textra::asDiagonal(state1.MPS_C), Textra::idx({0},{0}))
            .contract(Textra::asDiagonal(state2.MPS_C), Textra::idx({0},{0}));
    overlap = temp;



    auto mps1_it_R = state1.MPS_R.begin();
    auto mps2_it_R = state2.MPS_R.begin();
    while (mps1_it_R != state1.MPS_R.end() ){
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

double mpstools::finite::ops::expectation_value(const class_finite_chain_state &state1, const class_finite_chain_state &state2,const std::list<Eigen::Tensor<std::complex<double>,4>>  &mpos, Eigen::Tensor<std::complex<double>,3> Ledge, Eigen::Tensor<std::complex<double>,3> Redge){

    assert(state1.get_length() == state2.get_length() and "ERROR: States have different lengths! Can't do overlap.");
    assert(state1.get_position() == state2.get_position() and "ERROR: States need to be at the same position! Can't do overlap.");
    auto mps1_it_L = state1.MPS_L.begin();
    auto mps2_it_L = state2.MPS_L.begin();
    auto mpo_it    = mpos.begin();
    Eigen::Tensor<Scalar,3> L = Ledge;
    while(mps1_it_L != state1.MPS_L.end()){
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
        auto &MPS_C1 = state1.MPS_C;
        auto &MPS_C2 = state2.MPS_C;
        Eigen::Tensor<Scalar,3> temp =
                L
                .contract(asDiagonal(MPS_C1), idx({0}, {0}))
                .contract(asDiagonal(MPS_C2), idx({0}, {0}))
                .shuffle(array3{1, 2, 0});
        L = temp;
    }
    //Contract the right half of the state
    Eigen::Tensor<Scalar,3> R = Redge;
    auto mps1_it_R = state1.MPS_R.begin();
    auto mps2_it_R = state2.MPS_R.begin();
    while(mps1_it_R != state1.MPS_R.end()){
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

double mpstools::finite::ops::exp_sq_value(const class_finite_chain_state &state1, const class_finite_chain_state &state2,const std::list<Eigen::Tensor<std::complex<double>,4>>  &mpos, const Eigen::Tensor<std::complex<double>,4> Ledge, const Eigen::Tensor<std::complex<double>,4> Redge){

    assert(state1.get_length() == state2.get_length() and "ERROR: States have different lengths! Can't do overlap.");
    assert(state1.get_position() == state2.get_position() and "ERROR: States need to be at the same position! Can't do overlap.");
    auto mps1_it_L = state1.MPS_L.begin();
    auto mps2_it_L = state2.MPS_L.begin();
    auto mpo_it    = mpos.begin();
    Eigen::Tensor<Scalar,4> L = Ledge;
    while(mps1_it_L != state1.MPS_L.end()){
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
        auto &MPS_C1 = state1.MPS_C;
        auto &MPS_C2 = state2.MPS_C;
        Eigen::Tensor<Scalar,4> temp =
                L
                        .contract(asDiagonal(MPS_C1), idx({0}, {0}))
                        .contract(asDiagonal(MPS_C2), idx({0}, {0}))
                        .shuffle(array4{2,3,0,1});
        L = temp;
    }
    //Contract the right half of the state
    Eigen::Tensor<Scalar,4> R = Redge;
    auto mps1_it_R = state1.MPS_R.begin();
    auto mps2_it_R = state2.MPS_R.begin();
    while(mps1_it_R != state1.MPS_R.end()){
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
