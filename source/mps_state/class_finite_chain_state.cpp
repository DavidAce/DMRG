//
// Created by david on 2019-01-29.
//


#include "class_finite_chain_state.h"
#include <mps_tools/nmspc_mps_tools.h>
#include <general/nmspc_quantum_mechanics.h>

void class_finite_chain_state::clear(){
    *this = class_finite_chain_state();
}


void class_finite_chain_state::do_all_measurements(){
    using namespace mpstools::finite;
    measurements.length                         = measure::length(*this);
    measurements.bond_dimension_current         = measure::bond_dimension_current(*this);
    measurements.bond_dimension_midchain        = measure::bond_dimension_midchain(*this);
    measurements.bond_dimensions                = measure::bond_dimensions(*this);
    measurements.norm                           = measure::norm(*this);
    measurements.energy                     = measure::energy(*this);  //This number is needed for variance calculation!
    measurements.energy_per_site                = measure::energy_per_site(*this);
    measurements.energy_variance_mpo            = measure::energy_variance(*this);
    measurements.energy_variance_per_site       = measure::energy_variance_per_site(*this);
    measurements.entanglement_entropy_current   = measure::entanglement_entropy_current (*this);
    measurements.entanglement_entropy_midchain  = measure::entanglement_entropy_midchain(*this);
    measurements.entanglement_entropies         = measure::entanglement_entropies(*this);
    measurements.spin_components                = measure::spin_components(*this);
}


void class_finite_chain_state::set_positions(){

    size_t pos = 0;
    for (auto &MPS: MPS_L){MPS.set_position(pos++);}
    for (auto &MPS: MPS_R){MPS.set_position(pos++);}
    pos = 0;
    for (auto &ENV: ENV_L){ENV.set_position(pos++);}
    for (auto &ENV: ENV_R){ENV.set_position(pos++);}
    pos = 0;
    for (auto &ENV2: ENV2_L){ENV2.set_position(pos++);}
    for (auto &ENV2: ENV2_R){ENV2.set_position(pos++);}
    pos = 0;
    for (auto &MPO : MPO_L){MPO->set_position(pos++);}
    for (auto &MPO : MPO_R){MPO->set_position(pos++);}

}

size_t class_finite_chain_state::get_length()    const {return MPS_L.size() + MPS_R.size();}
size_t class_finite_chain_state::get_position()  const {return MPS_L.size() - 1u;}

int class_finite_chain_state::get_sweeps()    const {return num_sweeps;}
int class_finite_chain_state::reset_sweeps()  {num_sweeps = 0; return num_sweeps;}

int  class_finite_chain_state::get_direction() const {return direction;}
void class_finite_chain_state::flip_direction() {direction *= -1;}


bool class_finite_chain_state::position_is_the_middle() const {
    return (size_t) get_position() + 1 == (size_t)(get_length() / 2.0) and direction == 1 ;
}
bool class_finite_chain_state::position_is_the_middle_any_direction() const {
    return (size_t) get_position() + 1 == (size_t)(get_length() / 2.0);
}

bool class_finite_chain_state::position_is_the_left_edge() const {
    return get_position() == 0;
}

bool class_finite_chain_state::position_is_the_right_edge() const {
    return get_position() == get_length() - 2;
}

bool class_finite_chain_state::position_is_any_edge() const {
    return position_is_the_left_edge() or position_is_the_right_edge();
}

bool class_finite_chain_state::position_is_at(size_t pos)const{
    return get_position() == pos;
}


//With indices
const class_hamiltonian_base & class_finite_chain_state::get_MPO(size_t pos) const {
    if (pos >= get_length()){throw std::range_error("get_MPO(pos) pos out of range: " + std::to_string(pos));}
    if(pos <= MPO_L.size()-1){
        return *std::next(MPO_L.begin(), pos)->get();
    }else{
        pos -= MPS_L.size();
        return *std::next(MPO_R.begin(), pos)->get();
    }
}


const Eigen::Tensor<class_finite_chain_state::Scalar,3> & class_finite_chain_state::get_G(size_t pos)const{
    if (pos >= get_length()){throw std::range_error("get_G(pos) pos out of range: " + std::to_string(pos));}
    if(pos <= MPS_L.size()-1){
        return std::next(MPS_L.begin(), pos)->get_G();
    }else{
        pos -= MPS_L.size();
        return std::next(MPS_R.begin(), pos)->get_G();
    }
}

const Eigen::Tensor<class_finite_chain_state::Scalar,1> & class_finite_chain_state::get_L(size_t pos) const {
    if (pos >= get_length()+1 ){throw std::range_error("get_L(pos):  pos out of range: " + std::to_string(pos));}
    if(pos <= MPS_L.size()-1){
        return std::next(MPS_L.begin(), pos)->get_L();
    }else if (pos == MPS_L.size()){
        return MPS_C;
    }
    else {
        pos -= (MPS_L.size() + 1);
        return std::next(MPS_R.begin(), pos)->get_L();
    }
}


Eigen::Tensor<class_finite_chain_state::Scalar,3> & class_finite_chain_state::get_G(size_t pos){
    if (pos >= get_length()){throw std::range_error("get_G(pos) pos out of range: " + std::to_string(pos));}
    if(pos <= MPS_L.size()-1){
        return std::next(MPS_L.begin(), pos)->get_G();
    }else{
        pos -= MPS_L.size();
        return std::next(MPS_R.begin(), pos)->get_G();
    }
}

Eigen::Tensor<class_finite_chain_state::Scalar,1> & class_finite_chain_state::get_L(size_t pos) {
    if (pos >= get_length()+1 ){throw std::range_error("get_L(pos):  pos out of range: " + std::to_string(pos));}
    if(pos <= MPS_L.size()-1){
        return std::next(MPS_L.begin(), pos)->get_L();
    }else if (pos == MPS_L.size()){
        return MPS_C;
    }
    else {
        pos -= (MPS_L.size() + 1);
        return std::next(MPS_R.begin(), pos)->get_L();
    }
}




Eigen::Tensor<class_finite_chain_state::Scalar,3> class_finite_chain_state::get_A(size_t pos) const {
    if(pos > get_length()-2 )  {throw std::range_error("get_A(pos): pos larger than length-2");}
    return Textra::asDiagonal(get_L(pos)).contract(get_G(pos), Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
}

Eigen::Tensor<class_finite_chain_state::Scalar,3> class_finite_chain_state::get_B(size_t pos) const {
    if(pos > get_length()-1 )  {throw std::range_error("get_B(pos): pos larger than length-1");}
    return get_G(pos).contract(Textra::asDiagonal(get_L(pos+1)), Textra::idx({2},{0}));
}



//std::tuple<long,long,long> class_finite_chain_state::get_dims(size_t pos) const{
//    return {get_G(pos).dimension(0),get_G(pos).dimension(1),get_G(pos).dimension(2)};
//}



Eigen::Tensor<class_finite_chain_state::Scalar,3> class_finite_chain_state::get_A() const{
    return Textra::asDiagonal(MPS_L.back().get_L()).contract(MPS_L.back().get_G(), Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
}

Eigen::Tensor<class_finite_chain_state::Scalar,3> class_finite_chain_state::get_B() const{
    return MPS_R.front().get_G().contract(Textra::asDiagonal(MPS_R.front().get_L()), Textra::idx({2},{0}));
}


Eigen::Tensor<class_finite_chain_state::Scalar,4> class_finite_chain_state::get_theta() const{
    return
    get_A()
    .contract(Textra::asDiagonal(MPS_C), Textra::idx({2},{0}))
    .contract(get_B(), Textra::idx({2},{1}));
}

Eigen::Tensor<class_finite_chain_state::Scalar,4> class_finite_chain_state::get_theta(size_t pos) const {
    if(pos > get_length()-2 )  {throw std::range_error("get_theta(pos): pos larger than length-2");}
    return get_A(pos)
            .contract(Textra::asDiagonal(get_L(pos+1)), Textra::idx({2},{0}))
            .contract(get_B(pos+1), Textra::idx({2},{1}));
}


std::vector<size_t> class_finite_chain_state::activate_sites(long threshold){
    return active_sites = mpstools::finite::multisite::generate_site_list(*this,threshold);
}

std::vector<size_t> class_finite_chain_state::active_dimensions() const{
    std::vector<size_t> dimensions;
    dimensions.emplace_back(get_G(active_sites.front()).dimension(1));
    for (auto & site : active_sites){
        dimensions.emplace_back(get_G(site).dimension(0));
    }
    dimensions.emplace_back(get_G(active_sites.back()).dimension(2));
    return dimensions;
}


void class_finite_chain_state::unset_measurements()const {
    measurements = Measurements();
}

void class_finite_chain_state::do_all_measurements()const {
    measurements.length                           = mpstools::finite::measure::length                        (*this);
    measurements.bond_dimension_current           = mpstools::finite::measure::bond_dimension_current        (*this);
    measurements.bond_dimension_midchain          = mpstools::finite::measure::bond_dimension_midchain       (*this);
    measurements.bond_dimensions                  = mpstools::finite::measure::bond_dimensions               (*this);
    measurements.norm                             = mpstools::finite::measure::norm                          (*this);
    measurements.energy                           = mpstools::finite::measure::energy                        (*this);
    measurements.energy_per_site                  = mpstools::finite::measure::energy_per_site               (*this);
    measurements.energy_variance_mpo              = mpstools::finite::measure::energy_variance               (*this);
    measurements.energy_variance_per_site         = mpstools::finite::measure::energy_variance_per_site      (*this);
    measurements.spin_components                  = mpstools::finite::measure::spin_components               (*this); // This will automatically measure sx,sy and sz as well
    measurements.entanglement_entropy_current     = mpstools::finite::measure::entanglement_entropy_current  (*this);
    measurements.entanglement_entropy_midchain    = mpstools::finite::measure::entanglement_entropy_midchain (*this);
    measurements.entanglement_entropies           = mpstools::finite::measure::entanglement_entropies        (*this);
}

