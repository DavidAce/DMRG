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
#include <general/nmspc_random_numbers.h>


using Scalar         = std::complex<double>;
using namespace Textra;

std::list<Eigen::Tensor<Scalar,4>> MPS_Tools::Finite::Ops::make_mpo_list (const std::list<std::shared_ptr<class_hamiltonian_base>> &mpos_L, const std::list<std::shared_ptr<class_hamiltonian_base>> &mpos_R){
    std::list<Eigen::Tensor<Scalar,4>> mpos;
    for(auto &mpo_L : mpos_L){
        mpos.push_back(mpo_L->MPO());
    }
    for(auto &mpo_R : mpos_R){
        mpos.push_back(mpo_R->MPO());
    }
    return mpos;
}


void MPS_Tools::Finite::Ops::apply_mpo(class_finite_chain_state &state,const Eigen::Tensor<Scalar,4> mpo,const Eigen::Tensor<Scalar,3> Ledge, const Eigen::Tensor<Scalar,3> Redge){
    std::list<Eigen::Tensor<Scalar,4>> mpos(state.get_length(), mpo);
    apply_mpos(state,mpos,Ledge,Redge);
}



void MPS_Tools::Finite::Ops::apply_mpos(class_finite_chain_state &state,const std::list<Eigen::Tensor<Scalar,4>> &mpos,const Eigen::Tensor<Scalar,3> Ledge, const Eigen::Tensor<Scalar,3> Redge){
    if (mpos.size() != state.get_length()) throw std::runtime_error("Number of mpo's doesn't match the number of sites on the system");
    // Apply MPO's on Gamma matrices and
    // increase the size on all Lambdas by chi*mpoDim
    MPS_Tools::log->trace("Applying MPO's");
    auto mpo = mpos.begin();

    {
        for (auto & mps : state.get_MPS_L()){
            long mpoDimL = mpo->dimension(0);
            long mpoDimR = mpo->dimension(1);
            auto [d,chiL,chiR] = mps.get_dims();
//            std::cout << "L     " << mps.get_L().dimensions() <<"\n" << mps.get_L() << std::endl;
//            std::cout << "G     " << mps.get_G().dimensions() <<"\n" << mps.get_G() << std::endl;
            Eigen::Tensor<Scalar,3> G_temp =
                    mps.get_G()
                            .contract(*mpo, idx({0},{2}))
                            .shuffle(array5{4,0,2,1,3})
                            .reshape(array3{d, chiL*mpoDimL, chiR*mpoDimR});

            Eigen::Tensor<Scalar,1> L_temp = mps.get_L().broadcast(array1{mpoDimL});
            Scalar Lnorm = 1.0;//Textra::Tensor1_to_Vector(L_temp).norm();

            mps.set_G(G_temp);
            mps.set_L(L_temp/Lnorm);
            mpo++;
        }
    }

    //Take care of the middle
    {
        long mpoDim = mpo->dimension(0);
        Eigen::Tensor<Scalar, 1> C_temp = state.get_MPS_C().broadcast(array1{mpoDim});
        Scalar Cnorm = 1.0;//Textra::Tensor1_to_Vector(C_temp).norm();
        state.get_MPS_C() = C_temp / Cnorm;
    }
    {
        for (auto & mps : state.get_MPS_R()){
            long mpoDimL = mpo->dimension(0);
            long mpoDimR = mpo->dimension(1);
            auto [d,chiL,chiR] = mps.get_dims();
//            std::cout << "G     " << mps.get_G().dimensions() <<"\n" << mps.get_G() << std::endl;
//            std::cout << "L     " << mps.get_L().dimensions() <<"\n" << mps.get_L() << std::endl;
            Eigen::Tensor<Scalar,3> G_temp =
                    mps.get_G()
                            .contract(*mpo, idx({0},{2}))
                            .shuffle(array5{4,0,2,1,3})
                            .reshape(array3{d, chiL*mpoDimL, chiR*mpoDimR});

            Eigen::Tensor<Scalar,1> L_temp = mps.get_L().broadcast(array1{mpoDimR});
            Scalar Lnorm = 1.0;//Textra::Tensor1_to_Vector(L_temp).norm();

            mps.set_G(G_temp);
            mps.set_L(L_temp / Lnorm);
            mpo++;
        }
    }

    // Take care of the edges. Apply the left and right MPO-edges on A's and B's
    // so the left and right most lambdas become ones
    {
        long mpoDimL = mpos.front().dimension(0);
        auto Ldim    = Ledge.dimension(0);
        Scalar norm  = std::sqrt(Ldim);
        Eigen::Tensor<Scalar, 3> G_temp =
                Ledge
                        .shuffle(Textra::array3{0,2,1})
                        .reshape(Textra::array2{Ldim*mpoDimL,Ldim})
                        .contract(state.get_MPS_L().front().get_A(),Textra::idx({0},{1}))
                        .shuffle(Textra::array3{1,0,2});
        state.get_MPS_L().front().set_G(G_temp);
        state.get_MPS_L().front().set_L(Eigen::Tensor<Scalar,1>(Ldim).constant(1.0)/norm);
    }
    {
        long mpoDimR = mpos.back().dimension(1);
        auto Rdim    = Redge.dimension(0);
        Scalar norm  = std::sqrt(Rdim);
        Eigen::Tensor<Scalar, 3> G_temp =
                Redge
                        .shuffle(Textra::array3{0,2,1})
                        .reshape(Textra::array2{Rdim*mpoDimR,Rdim})
                        .contract(state.get_MPS_R().back().get_B(),Textra::idx({0},{2}))
                        .shuffle(Textra::array3{1,2,0});
        state.get_MPS_R().back().set_G(G_temp);
        state.get_MPS_R().back().set_L(Eigen::Tensor<Scalar,1>(Rdim).constant(1.0)/norm);
    }


}


void MPS_Tools::Finite::Ops::normalize_chain(class_finite_chain_state & state){

    MPS_Tools::log->trace("Normalizing chain");

    // Sweep back and forth once on the chain

    class_SVD<Scalar> svd;
    size_t pos_LA = 0; // Lambda left of GA
    size_t pos_LC = 1; //Center Lambda
    size_t pos_LB = 2; // Lambda right of GB
    size_t pos_A  = 0; // Lambda * G
    size_t pos_B  = 1; // G * Lambda
//    bool finished = false;
    size_t num_traversals = 0;
    int direction = 1;
    Eigen::Tensor<Scalar,3> U;
    Eigen::Tensor<Scalar,1> S;
    Eigen::Tensor<Scalar,3> V;
    double norm;
    while(num_traversals < 2){
        Eigen::Tensor<Scalar,4> theta = state.get_theta(pos_A);
        try {std::tie(U,S,V,norm) = svd.schmidt_with_norm(theta);}
        catch(std::exception &ex){
            std::cerr << "A:\n" << state.get_A(pos_A) << std::endl;
            std::cerr << "L:\n" << state.get_L(pos_LC) << std::endl;
            std::cerr << "B:\n" << state.get_B(pos_B) << std::endl;
            throw std::runtime_error("Normalization failed at step " + std::to_string(pos_A) + ": " + std::string(ex.what()) );
        }

        Eigen::Tensor<Scalar,3> LA_U = Textra::asDiagonalInversed(state.get_L(pos_LA)).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
        Eigen::Tensor<Scalar,3> V_LB = V.contract(asDiagonalInversed(state.get_L(pos_LB)), idx({2},{0}));
        norm = pos_B == state.get_length()-1 or pos_A == 0 ? 1.0 : norm;

        if (direction == 1){
            state.get_G(pos_A)  = LA_U;
            state.get_L(pos_LC) = S;
            state.get_G(pos_B)  = Scalar(norm)*V_LB;

        }
        if (direction == -1){
            state.get_G(pos_A)  = LA_U * Scalar(norm);
            state.get_L(pos_LC) = S;
            state.get_G(pos_B)  = V_LB;
        }

        if (direction ==  1 and pos_B == state.get_length()-1)  {num_traversals++; direction *= -1;}
        if (direction == -1 and pos_A == 0)                     {num_traversals++; direction *= -1;}

        pos_LA += direction;
        pos_LC += direction;
        pos_LB += direction;
        pos_A  += direction;
        pos_B  += direction;


    }




}



void MPS_Tools::Finite::Ops::reset_to_random_product_state(class_finite_chain_state &state, const std::string parity){
    std::array<Eigen::Tensor <Scalar,1>,2> paulieigvec;
    Eigen::Matrix2cd paulimatrix;
    bool no_parity = false;
    if(parity == "sx"){
        paulimatrix = qm::spinOneHalf::sx;
        paulieigvec[0] = Textra::Matrix_to_Tensor1(qm::spinOneHalf::sx_eigvecs[0]);
        paulieigvec[1] = Textra::Matrix_to_Tensor1(qm::spinOneHalf::sx_eigvecs[1]);
    }else if (parity == "sy"){
        paulimatrix = qm::spinOneHalf::sy;
        paulieigvec[0] = Textra::Matrix_to_Tensor1(qm::spinOneHalf::sy_eigvecs[0]);
        paulieigvec[1] = Textra::Matrix_to_Tensor1(qm::spinOneHalf::sy_eigvecs[1]);
    }else if (parity == "sz"){
        paulimatrix = qm::spinOneHalf::sz;
        paulieigvec[0] = Textra::Matrix_to_Tensor1(qm::spinOneHalf::sz_eigvecs[0]);
        paulieigvec[1] = Textra::Matrix_to_Tensor1(qm::spinOneHalf::sz_eigvecs[1]);
    }else if (parity == "random" or parity == "none"){
        no_parity = true;
    }else {
        throw std::runtime_error("Invalid spin parity name: " + parity);
    }

    Eigen::Tensor<Scalar,3> G (2,1,1);
    Eigen::Tensor<Scalar,1> L (1);
    L.setConstant(1);

    for (auto &mpsL : state.get_MPS_L() ){
        if(no_parity){
            G.setRandom();
        }else{
            auto rnd = rn::uniform_integer_1();
            G = paulieigvec[rnd].reshape(Textra::array3{2,1,1});
        }
        mpsL.set_mps(G,L);
    }
    auto C = Eigen::Tensor<Scalar,1>(1);
    C.setConstant(1);
    state.get_MPS_C() = C;
    for (auto &mpsR : state.get_MPS_R() ){
        if(no_parity){
            G.setRandom();
        }else{
            auto rnd = rn::uniform_integer_1();
            G = paulieigvec[rnd].reshape(Textra::array3{2,1,1});
        }
        mpsR.set_mps(G,L);
    }

    normalize_chain(state);
    rebuild_environments(state);
    MPS_Tools::Finite::Debug::check_integrity_of_mps(state);
}




class_finite_chain_state MPS_Tools::Finite::Ops::get_parity_projected_state(const class_finite_chain_state &state, const Eigen::MatrixXcd  paulimatrix, const int sign) {
    if (std::abs(sign) != 1) throw std::runtime_error("Expected 'sign' +1 or -1. Got: " + std::to_string(sign));
    MPS_Tools::log->trace("Generating parity projected state");
    class_finite_chain_state state_projected = state;
    const auto [mpo,L,R]    = class_mpo::parity_projector_mpos(paulimatrix,state_projected.get_length(), sign);
    apply_mpos(state_projected,mpo, L,R);
    normalize_chain(state_projected);
    rebuild_environments(state_projected);
    MPS_Tools::Finite::Debug::check_integrity_of_mps(state_projected);
    return state_projected;
}

class_finite_chain_state MPS_Tools::Finite::Ops::get_closest_parity_state(const class_finite_chain_state &state, const Eigen::MatrixXcd  paulimatrix) {
    MPS_Tools::log->trace("Finding closest projection");
    double measured_spin_component = MPS_Tools::Finite::Measure::spin_component(state, paulimatrix);
    if (measured_spin_component > 0){
        return get_parity_projected_state(state, paulimatrix, 1);
    }else{
        return get_parity_projected_state(state, paulimatrix,-1);
    }
}

class_finite_chain_state MPS_Tools::Finite::Ops::get_closest_parity_state(const class_finite_chain_state &state, const std::string paulistring) {
    MPS_Tools::log->trace("Finding closest projection");
    Eigen::MatrixXcd paulimatrix;
    if (paulistring == "sx"){paulimatrix = qm::spinOneHalf::sx;}
    if (paulistring == "sy"){paulimatrix = qm::spinOneHalf::sy;}
    if (paulistring == "sz"){paulimatrix = qm::spinOneHalf::sz;}

    double measured_spin_component = MPS_Tools::Finite::Measure::spin_component(state, paulimatrix);
    if (measured_spin_component > 0){
        return get_parity_projected_state(state, paulimatrix, 1);
    }else{
        return get_parity_projected_state(state, paulimatrix,-1);
    }
}

void MPS_Tools::Finite::Ops::rebuild_superblock(class_finite_chain_state &state, class_superblock &superblock) {
    MPS_Tools::log->trace("Rebuilding superblock");
    MPS_Tools::Finite::Chain::copy_state_to_superblock(state,superblock);
}

void MPS_Tools::Finite::Ops::rebuild_environments(class_finite_chain_state &state){
    MPS_Tools::log->trace("Rebuilding environments");

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
        auto ENV_L = class_environment("L",state.get_MPS_L().front().get_chiL(), state.get_MPO_L().front()->MPO().dimension(0));
        auto ENV2_L = class_environment_var("L",state.get_MPS_L().front().get_chiL(), state.get_MPO_L().front()->MPO().dimension(0));
        auto mpsL_it   = state.get_MPS_L().begin();
        auto mpoL_it   = state.get_MPO_L().begin();

        while(mpsL_it != state.get_MPS_L().end() and mpoL_it != state.get_MPO_L().end()) {
            state.get_ENV_L().push_back(ENV_L);
            state.get_ENV2_L().push_back(ENV2_L);
            if (mpsL_it == state.get_MPS_L().end()) break;
            ENV_L.enlarge(mpsL_it->get_A(), mpoL_it->get()->MPO());
            ENV2_L.enlarge(mpsL_it->get_A(), mpoL_it->get()->MPO());
            mpsL_it++;
            mpoL_it++;
        }
    }



    {
        auto ENV_R  = class_environment("R",state.get_MPS_R().back().get_chiR(), state.get_MPO_L().back()->MPO().dimension(1));
        auto ENV2_R = class_environment_var("R",state.get_MPS_R().back().get_chiR(), state.get_MPO_L().back()->MPO().dimension(1));
        auto mpsR_it   = state.get_MPS_R().rbegin();
        auto mpoR_it   = state.get_MPO_R().rbegin();
        while(mpsR_it != state.get_MPS_R().rend() and mpoR_it != state.get_MPO_R().rend()){
            state.get_ENV_R().push_front(ENV_R);
            state.get_ENV2_R().push_front(ENV2_R);
            if (mpsR_it == state.get_MPS_R().rend()) break;
            ENV_R.enlarge(mpsR_it->get_B(), mpoR_it->get()->MPO());
            ENV2_R.enlarge(mpsR_it->get_B(), mpoR_it->get()->MPO());
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
