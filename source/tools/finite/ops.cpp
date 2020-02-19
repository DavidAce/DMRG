//
// Created by david on 2019-01-30.
//

#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <math/nmspc_random.h>
#include <model/class_model_base.h>
#include <state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>
#include <tools/finite/debug.h>

using Scalar         = std::complex<double>;
using namespace Textra;

std::list<Eigen::Tensor<Scalar,4>> tools::finite::ops::make_mpo_list (const std::list<std::unique_ptr<class_model_base>> & mpos_L, const std::list<std::unique_ptr<class_model_base>> & mpos_R){
    std::list<Eigen::Tensor<Scalar,4>> mpos;
    for(auto &mpo_L : mpos_L){
        mpos.push_back(mpo_L->MPO());
    }
    for(auto &mpo_R : mpos_R){
        mpos.push_back(mpo_R->MPO());
    }
    return mpos;
}


void tools::finite::ops::apply_mpo(class_state_finite & state, const Eigen::Tensor<Scalar,4> &mpo, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge){
    std::list<Eigen::Tensor<Scalar,4>> mpos(state.get_length(), mpo);
    apply_mpos(state,mpos,Ledge,Redge);
}



void tools::finite::ops::apply_mpos(class_state_finite & state, const std::list<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge){
    if (mpos.size() != state.get_length()) throw std::runtime_error("Number of mpo's doesn't match the number of sites on the system");
    // Apply MPO's on Gamma matrices and
    // increase the size on all Lambdas by chi*mpoDim
    tools::log->trace("Applying MPO's");
    state.clear_measurements();
    if(state.hasNaN()) throw std::runtime_error("State has NAN's before applying MPO's");
    tools::log->trace("Bond dimensions before applying MPO's = {}", tools::finite::measure::bond_dimensions(state));
    tools::log->trace("Norm before applying MPO's = {:.16f}", tools::finite::measure::norm(state));

    auto mpo = mpos.begin();
    for(size_t pos = 0; pos < state.get_length() ; pos++){
        state.get_MPS(pos).apply_mpo(*mpo);
        mpo++;
    }


    // Take care of the edges. Apply the left and right MPO-edges on A's and B's
    // so the left and right most lambdas become ones
    {
        long mpoDimL = mpos.front().dimension(0);
        auto Ldim    = Ledge.dimension(0);
//        auto dims = state.MPS_L.front().get_M().dimensions();
//        std::cout << "Mdims = " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
        Eigen::Tensor<Scalar, 3> M_temp =
                Ledge
                        .shuffle(Textra::array3{0,2,1})
                        .reshape(Textra::array2{Ldim*mpoDimL,Ldim})
                        .contract(state.MPS_L.front().get_M_bare(),Textra::idx({0},{1}))
                        .shuffle(Textra::array3{1,0,2});
        state.MPS_L.front().set_M(M_temp);
        state.MPS_L.front().set_L(Eigen::Tensor<Scalar,1>(Ldim).constant(1.0));
    }
    {
        long mpoDimR = mpos.back().dimension(1);
        auto Rdim    = Redge.dimension(0);
        Eigen::Tensor<Scalar, 3> M_temp =
                Redge
                        .shuffle(Textra::array3{0,2,1})
                        .reshape(Textra::array2{Rdim*mpoDimR,Rdim})
                        .contract(state.MPS_R.back().get_M_bare(),Textra::idx({0},{2}))
                        .shuffle(Textra::array3{1,2,0});
        state.MPS_R.back().set_M(M_temp);
        state.MPS_R.back().set_L(Eigen::Tensor<Scalar,1>(Rdim).constant(1.0));
    }
    state.clear_measurements();
    if(state.hasNaN()) throw std::runtime_error("State has NAN's after applying MPO's");
    tools::log->trace("Bond dimensions after applying MPO's = {}", tools::finite::measure::bond_dimensions(state));
    tools::log->trace("Norm after  applying MPO's = {:.16f}", tools::finite::measure::norm(state));//    std::cout << "Norm              (after mpos): " << tools::finite::measure::norm(state)  << std::endl;
//    std::cout << "Spin component sx (after mpos): " << tools::finite::measure::spin_component(state, qm::spinOneHalf::sx)  << std::endl;
//    std::cout << "Spin component sy (after mpos): " << tools::finite::measure::spin_component(state, qm::spinOneHalf::sy)  << std::endl;
//    std::cout << "Spin component sz (after mpos): " << tools::finite::measure::spin_component(state, qm::spinOneHalf::sz)  << std::endl;
}


class_state_finite tools::finite::ops::get_projection_to_parity_sector(const class_state_finite & state, const Eigen::MatrixXcd  & paulimatrix, int sign) {
    if (std::abs(sign) != 1) throw std::runtime_error("Expected 'sign' +1 or -1. Got: " + std::to_string(sign));
    tools::log->debug("Generating parity projected state with sign {}", sign);
    auto spin_components = tools::finite::measure::spin_components(state);
    double requested_spin_component = tools::finite::measure::spin_component(state, paulimatrix);
    tools::log->debug("Current global spin components : X = {:.16f}  Y = {:.16f}  Z = {:.16f}",spin_components[0],spin_components[1],spin_components[2] );
    tools::log->debug("Current reqstd spin component  :     {:.16f}", requested_spin_component );

    tools::common::profile::t_prj->tic();
    class_state_finite state_projected = state;
    state_projected.clear_measurements();
    state_projected.clear_cache();

    const auto [mpo,L,R]    = qm::mpo::parity_projector_mpos(paulimatrix,state_projected.get_length(), sign);
    apply_mpos(state_projected,mpo, L,R);
    tools::common::profile::t_prj->toc();
    tools::finite::mps::normalize(state_projected);
    tools::finite::mps::rebuild_environments(state_projected);
    tools::finite::debug::check_integrity_of_mps(state_projected);
    state_projected.tag_all_sites_have_been_updated(true); // All sites change in this operation
    spin_components          = tools::finite::measure::spin_components(state_projected);
    requested_spin_component = tools::finite::measure::spin_component(state_projected, paulimatrix);
    tools::log->debug("Resulting global spin components : X = {:.16f}  Y = {:.16f}  Z = {:.16f}",spin_components[0],spin_components[1],spin_components[2] );
    tools::log->debug("Resulting reqstd spin component  :     {:.16f}", requested_spin_component );
    return state_projected;
}

class_state_finite tools::finite::ops::get_projection_to_closest_parity_sector(const class_state_finite &state, const Eigen::MatrixXcd & paulimatrix) {
    tools::log->trace("Finding closest projection");
    double requested_spin_component = tools::finite::measure::spin_component(state, paulimatrix);
    if (requested_spin_component > 0){
        return get_projection_to_parity_sector(state, paulimatrix, 1);
    }else{
        return get_projection_to_parity_sector(state, paulimatrix, -1);
    }
}

class_state_finite tools::finite::ops::get_projection_to_closest_parity_sector(const class_state_finite &state, std::string parity_sector) {
    tools::log->trace("Finding closest projection in parity sector {}", parity_sector );
    if      (parity_sector == "x")  {return get_projection_to_closest_parity_sector(state, qm::spinOneHalf::sx);}
    else if (parity_sector == "y")  {return get_projection_to_closest_parity_sector(state, qm::spinOneHalf::sy);}
    else if (parity_sector == "z")  {return get_projection_to_closest_parity_sector(state, qm::spinOneHalf::sz);}
    else if (parity_sector == "+x") {return get_projection_to_parity_sector(state, qm::spinOneHalf::sx, 1);}
    else if (parity_sector == "-x") {return get_projection_to_parity_sector(state, qm::spinOneHalf::sx,-1);}
    else if (parity_sector == "+y") {return get_projection_to_parity_sector(state, qm::spinOneHalf::sy, 1);}
    else if (parity_sector == "-y") {return get_projection_to_parity_sector(state, qm::spinOneHalf::sy,-1);}
    else if (parity_sector == "+z") {return get_projection_to_parity_sector(state, qm::spinOneHalf::sz, 1);}
    else if (parity_sector == "-z") {return get_projection_to_parity_sector(state, qm::spinOneHalf::sz,-1);}
    else if (parity_sector == "randomAxis"){
        std::vector<std::string> possibilities = {"x","y","z"};
        std::string chosen_axis = possibilities[rn::uniform_integer_box(0, 2)];
        get_projection_to_closest_parity_sector(state, chosen_axis);
    }
    else if (parity_sector == "random"){
        auto coeffs = Eigen::Vector3d::Random().normalized();
        Eigen::Matrix2cd random_c2 =
                    coeffs(0) * qm::spinOneHalf::sx
                +   coeffs(1) * qm::spinOneHalf::sy
                +   coeffs(2) * qm::spinOneHalf::sz;
        return get_projection_to_closest_parity_sector(state, random_c2);
    }
    else if (parity_sector == "none"){return state;}
    else{
        tools::log->warn(R"(Wrong pauli string. Expected one of (+-) "x","y","z", "randomAxis", "random" or "none". Got: )" + parity_sector);
        tools::log->warn("Taking whichever is closest to current state!");
        auto spin_components = tools::finite::measure::spin_components(state);
        auto max_idx = std::distance(spin_components.begin(), std::max_element(spin_components.begin(),spin_components.end()));
        if(max_idx == 0)      {return get_projection_to_closest_parity_sector(state, "x"); }
        else if(max_idx == 1) {return get_projection_to_closest_parity_sector(state, "y"); }
        else if(max_idx == 2) {return get_projection_to_closest_parity_sector(state, "z"); }
        else {throw std::runtime_error("Wrong parity_sector string and could not find closest parity state");}
    }
    throw std::runtime_error(fmt::format(R"(Wrong pauli string. Expected one of (+-) "x","y","z", "randomAxis", "random" or "none". Got: )" + parity_sector));
}


double tools::finite::ops::overlap(const class_state_finite & state1, const class_state_finite & state2){

    assert(state1.get_length() == state2.get_length() and "ERROR: States have different lengths! Can't do overlap.");
    assert(state1.get_position() == state2.get_position() and "ERROR: States need to be at the same position! Can't do overlap.");

    Eigen::Tensor<Scalar,2> overlap =
            state1.MPS_L.front().get_M()
            .contract(state2.MPS_L.front().get_M().conjugate(), Textra::idx({0,1},{0,1}));
    for(size_t pos = 0; pos < state1.get_length(); pos++){
        Eigen::Tensor<Scalar,2> temp = overlap
                .contract(state1.get_MPS(pos).get_M()            , Textra::idx({0},{1}))
                .contract(state2.get_MPS(pos).get_M().conjugate(), Textra::idx({0,1},{1,0}));
        overlap = temp;
    }

    double norm_chain = std::real(Textra::TensorMatrixMap(overlap).trace());
    return norm_chain;
}

double tools::finite::ops::expectation_value(const class_state_finite & state1, const class_state_finite & state2, const std::list<Eigen::Tensor<std::complex<double>,4>>  & mpos, const Eigen::Tensor<std::complex<double>,3> & Ledge, const Eigen::Tensor<std::complex<double>,3> & Redge){

    assert(state1.get_length() == state2.get_length() and "ERROR: States have different lengths! Can't do overlap.");
    assert(state1.get_position() == state2.get_position() and "ERROR: States need to be at the same position! Can't do overlap.");
    auto mpo_it    = mpos.begin();
    Eigen::Tensor<Scalar,3> L = Ledge;
    for(size_t pos = 0; pos < state1.get_length(); pos++) {
        Eigen::Tensor<Scalar,3> temp =
                L
                .contract(state1.get_MPS(pos).get_M()              , idx({0},{1}))
                .contract(*mpo_it++                                , idx({1,2},{0,2}))
                .contract(state2.get_MPS(pos).get_M().conjugate()  , idx({0,3},{1,0}))
                .shuffle(array3{0,2,1});

        L = temp;
    }
    assert(L.dimensions() == Redge.dimensions());
    Eigen::Tensor<Scalar,0> E_all_sites = L.contract(Redge, idx({0,1,2},{0,1,2}));
    double energy_chain = std::real(E_all_sites(0));
    return energy_chain;
}

double tools::finite::ops::exp_sq_value     (const class_state_finite & state1, const class_state_finite & state2, const std::list<Eigen::Tensor<std::complex<double>,4>>  & mpos, const Eigen::Tensor<std::complex<double>,4> & Ledge, const Eigen::Tensor<std::complex<double>,4> & Redge){

    assert(state1.get_length() == state2.get_length() and "ERROR: States have different lengths! Can't do overlap.");
    assert(state1.get_position() == state2.get_position() and "ERROR: States need to be at the same position! Can't do overlap.");
    auto mpo_it    = mpos.begin();
    Eigen::Tensor<Scalar,4> L = Ledge;
    for(size_t pos = 0; pos < state1.get_length(); pos++) {
        Eigen::Tensor<Scalar,4> temp =
                L
                .contract(state1.get_MPS(pos).get_M()              , idx({0},{1}))
                .contract(*mpo_it                                  , idx({1,3},{0,2}))
                .contract(*mpo_it++                                , idx({1,4},{0,2}))
                .contract(state2.get_MPS(pos).get_M().conjugate()  , idx({0,4},{1,0}))
                .shuffle(array4{0,3,1,2});

        L = temp;
    }
    assert(L.dimensions() == Redge.dimensions());
    Eigen::Tensor<Scalar,0> H2_all_sites = L.contract(Redge, idx({0,1,2,3},{0,1,2,3}));
    return std::real(H2_all_sites(0));
}
