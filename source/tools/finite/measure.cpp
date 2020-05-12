//
// Created by david on 2019-02-01.
//

//
// Created by david on 2017-11-12.
//


#include "measure.h"

#include <iomanip>
#include <simulation/nmspc_settings.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>
#include <simulation/class_simulation_status.h>
#include <state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>

using namespace std;
using namespace Textra;
using Scalar = class_state_finite::Scalar;



size_t tools::finite::measure::length(const class_state_finite & state){
    return state.get_length();
}

double tools::finite::measure::norm(const class_state_finite & state){
    if (state.measurements.norm){return state.measurements.norm.value();}
    Eigen::Tensor<Scalar,2> chain;
    Eigen::Tensor<Scalar,2> temp;
    bool first = true;
    for(size_t pos = 0; pos < state.get_length(); pos++){
        const Eigen::Tensor<Scalar,3>  &M = state.get_MPS(pos).get_M(); // std::get<1>(*mpsL);
        if(first) {chain = M.contract(M.conjugate(), idx({0,1},{0,1})); first=false;continue;}
        temp =
            chain
                .contract(M,             idx({0},{1}))
                .contract(M.conjugate(), idx({0,1},{1,0}));
        chain = temp;
    }
    double norm_chain = std::abs(Textra::TensorMatrixMap(chain).trace());
    if(std::abs(norm_chain - 1.0) > settings::precision::max_norm_error)
        tools::log->debug("Norm far from unity: {:.16f}", norm_chain);
    state.measurements.norm = norm_chain;
    return state.measurements.norm.value();
}

size_t tools::finite::measure::bond_dimension_current(const class_state_finite & state){
    if (state.measurements.bond_dimension_current){return state.measurements.bond_dimension_current.value();}
    if (state.MPS_L.back().get_chiR() != state.current_bond().dimension(0)) throw std::runtime_error("Center bond dimension mismatch!");
    state.measurements.bond_dimension_current = state.current_bond().dimension(0);
    return state.measurements.bond_dimension_current.value();
}


size_t tools::finite::measure::bond_dimension_midchain(const class_state_finite & state){
    if (state.measurements.bond_dimension_midchain){return state.measurements.bond_dimension_midchain.value();}
    state.measurements.bond_dimension_midchain = state.midchain_bond().dimension(0);
    return state.measurements.bond_dimension_midchain.value();
}


std::vector<size_t> tools::finite::measure::bond_dimensions(const class_state_finite & state){
    if (state.measurements.bond_dimensions){return state.measurements.bond_dimensions.value();}
    state.measurements.bond_dimensions = std::vector<size_t>{};
    for (size_t pos = 0; pos < state.get_length(); pos++){
        state.measurements.bond_dimensions.value().emplace_back(state.get_MPS(pos).get_L().dimension(0));
        if(state.get_MPS(pos).isCenter()){
            state.measurements.bond_dimensions.value().emplace_back(state.get_MPS(pos).get_LC().dimension(0));
        }
    }
    return state.measurements.bond_dimensions.value();
}



double tools::finite::measure::energy(const class_state_finite &state){
    if (state.measurements.energy)         return state.measurements.energy.value();
    if (state.active_sites.empty()){
        tools::common::profile::t_ene->tic();
        auto theta = state.get_theta();
        tools::common::profile::t_ene->toc();
        state.measurements.energy = twosite::energy(state,theta);
    }else
    {
        tools::common::profile::t_ene->tic();
        auto theta = state.get_multisite_mps();
        tools::common::profile::t_ene->toc();
        state.measurements.energy = multisite::energy(state,theta);
    }
    return state.measurements.energy.value();
}


double tools::finite::measure::energy_per_site(const class_state_finite &state){
    if (state.measurements.energy_per_site)return state.measurements.energy_per_site.value();
    if (state.active_sites.size() > 2)     return multisite::energy_per_site(state);
    state.measurements.energy_per_site = energy(state)/static_cast<double>(state.get_length());
    return state.measurements.energy_per_site.value();

}



double tools::finite::measure::energy_variance(const class_state_finite &state){
    if (state.measurements.energy_variance) return state.measurements.energy_variance.value();
    if (state.active_sites.empty()){
        tools::common::profile::t_var->tic();
        auto theta = state.get_theta();
        tools::common::profile::t_var->toc();
        state.measurements.energy_variance = twosite::energy_variance(state, theta);
    }else{
        tools::common::profile::t_var->tic();
        auto theta = state.get_multisite_mps();
        tools::common::profile::t_var->toc();
        state.measurements.energy_variance = multisite::energy_variance(state, theta);
    }
    return state.measurements.energy_variance.value();
}


double tools::finite::measure::energy_variance_per_site(const class_state_finite &state){
    if (state.measurements.energy_variance_per_site) return state.measurements.energy_variance_per_site.value();
    state.measurements.energy_variance_per_site = tools::finite::measure::energy_variance(state)/static_cast<double>(state.get_length());
    return state.measurements.energy_variance_per_site.value();
}


double tools::finite::measure::energy_normalized(const class_state_finite &state, const class_simulation_status &sim_status) {
    return  (tools::finite::measure::energy_per_site(state) - sim_status.energy_min ) / (sim_status.energy_max - sim_status.energy_min);

}


double tools::finite::measure::entanglement_entropy_current(const class_state_finite & state){
    if (state.measurements.entanglement_entropy_current){return state.measurements.entanglement_entropy_current.value();}
    tools::common::profile::t_ent->tic();
    auto & LC = state.current_bond();
    Eigen::Tensor<Scalar,0> SE  = -LC.square()
            .contract(LC.square().log().eval(), idx({0},{0}));
    state.measurements.entanglement_entropy_current = std::real(SE(0));
    tools::common::profile::t_ent->toc();
    return state.measurements.entanglement_entropy_current.value();
}

double tools::finite::measure::entanglement_entropy_midchain(const class_state_finite & state){
    if (state.measurements.entanglement_entropy_midchain){return state.measurements.entanglement_entropy_midchain.value();}
    tools::common::profile::t_ent->tic();
    auto & LC = state.midchain_bond();
    Eigen::Tensor<Scalar,0> SE  = -LC.square()
            .contract(LC.square().log().eval(), idx({0},{0}));
    state.measurements.entanglement_entropy_midchain =  std::real(SE(0));
    tools::common::profile::t_ent->toc();
    return state.measurements.entanglement_entropy_midchain.value();
}

std::vector<double> tools::finite::measure::entanglement_entropies(const class_state_finite & state){
    if (state.measurements.entanglement_entropies){return state.measurements.entanglement_entropies.value();}
    tools::common::profile::t_ent->tic();
    std::vector<double> entanglement_entropies;
    for (size_t pos = 0; pos < state.get_length(); pos++){
        auto &L = state.get_MPS(pos).get_L();
        Eigen::Tensor<Scalar, 0> SE = -L.square().contract(L.square().log().eval(), idx({0}, {0}));
        entanglement_entropies.emplace_back(std::real(SE(0)));
        if(state.get_MPS(pos).isCenter()){
            auto &LC = state.get_MPS(pos).get_LC();
            SE = -LC.square().contract(LC.square().log().eval(), idx({0}, {0}));
            entanglement_entropies.emplace_back(std::real(SE(0)));
            state.measurements.entanglement_entropy_current =  std::real(SE(0));
        }
    }
    state.measurements.entanglement_entropies = entanglement_entropies;
    tools::common::profile::t_ent->toc();
    return state.measurements.entanglement_entropies.value();
}


std::array<double,3> tools::finite::measure::spin_components(const class_state_finite &state){
    if (state.measurements.spin_components){return state.measurements.spin_components.value();}
    state.measurements.spin_component_sx                      = measure::spin_component(state, qm::spinOneHalf::sx);
    state.measurements.spin_component_sy                      = measure::spin_component(state, qm::spinOneHalf::sy);
    state.measurements.spin_component_sz                      = measure::spin_component(state, qm::spinOneHalf::sz);
    state.measurements.spin_components =  {state.measurements.spin_component_sx.value(),
                                           state.measurements.spin_component_sy.value(),
                                           state.measurements.spin_component_sz.value()};
    return state.measurements.spin_components.value();
}


double tools::finite::measure::spin_component(const class_state_finite &state,
                                                  const Eigen::Matrix2cd & paulimatrix){

    auto [mpo,L,R]   = qm::mpo::pauli_mpo(paulimatrix);
    Eigen::TensorRef<Eigen::Tensor<Scalar,3>> temp;
    for (size_t pos = 0; pos < state.get_length(); pos++){
        temp = L.contract(state.get_MPS(pos).get_M()             , idx({0},{1}))
                .contract(state.get_MPS(pos).get_M().conjugate() , idx({0},{1}))
                .contract(mpo                                    , idx({0,1,3},{0,2,3}));
        L = temp;
    }

    assert(L.dimensions() == R.dimensions());
    Eigen::Tensor<Scalar,0> parity_tmp = L.contract(R, idx({0,1,2},{0,1,2}));
    double parity = std::real(parity_tmp(0));
    return parity;

}


double tools::finite::measure::spin_component(const class_state_finite &state,
                                              const std::string & parity_sector){
    if (parity_sector.find('x') != parity_sector.npos)
        return measure::spin_component(state, qm::spinOneHalf::sx);
    if (parity_sector.find('y') != parity_sector.npos)
        return measure::spin_component(state, qm::spinOneHalf::sy);
    if (parity_sector.find('z') != parity_sector.npos)
        return measure::spin_component(state, qm::spinOneHalf::sz);

    throw std::runtime_error("Unexpected parity sector: " + parity_sector);

}



Eigen::Tensor<Scalar,1> tools::finite::measure::mps_wavefn(const class_state_finite & state){

    Eigen::Tensor<Scalar,2> chain(1,1);
    chain.setConstant(1.0);
    Eigen::TensorRef<Eigen::Tensor<Scalar,2>> temp;
    // The "state" is a matrix whose 0 index keeps growing.
    // For each site that passes, it grows by GA.dimension(0) = phys dim
    // Say the state is a 16x7 matrix (having contracted 4 particles, and the latest
    // chi was 7). Then contracting the next site, with dimensions 2x7x9 will get you a
    // 16x2x9 tensor. Now the reshaping convert it into a 32 x 9 matrix. Because
    // Eigen is column major, the doubling 16->32 will stack the third index twice.

    for(auto & mpsL : state.MPS_L){
        long dim0 = mpsL.get_spin_dim();
        long dimR = mpsL.get_chiR();
        long dimL = chain.dimension(0);
        temp = chain
                .contract(mpsL.get_M(), idx({1},{1}))
                .reshape(array2{dimL * dim0, dimR});
        chain = temp;
    }
    for(auto & mpsR : state.MPS_R){
        long dim0 = mpsR.get_spin_dim();
        long dimR = mpsR.get_chiR();
        long dimL = chain.dimension(0);
        temp = chain
                .contract(mpsR.get_M(), idx({1},{1}))
                .reshape(array2{dimL * dim0, dimR});
        chain = temp;
    }

    Eigen::Tensor<Scalar,1> mps_chain = chain.reshape(array1{chain.dimension(0)});
    double norm_chain = Textra::TensorVectorMap(chain).norm();
    if(std::abs(norm_chain - 1.0) > settings::precision::max_norm_error){
        tools::log->warn("Norm far from unity: {}", norm_chain);
        throw std::runtime_error("Norm too far from unity: " + std::to_string(norm_chain));
    }
    return mps_chain;
}
