//
// Created by david on 2019-01-29.
//


#include "class_finite_chain_state.h"
#include <mps_tools/nmspc_mps_tools.h>
#include <general/nmspc_quantum_mechanics.h>


void class_finite_chain_state::do_all_measurements(){
    using namespace mpstools::finite;
    measurements.length                         = measure::length(*this);
    measurements.bond_dimension_current         = measure::bond_dimension_current(*this);
    measurements.bond_dimension_midchain        = measure::bond_dimension_midchain(*this);
    measurements.bond_dimensions                = measure::bond_dimensions(*this);
    measurements.norm                           = measure::norm(*this);
    measurements.energy                         = measure::energy(*this);  //This number is needed for variance calculation!
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

bool class_finite_chain_state::isReal() const{
    bool mps_real = true;
    bool mpo_real = true;
    for(auto & mps : MPS_L ){mps_real = mps_real and mps.isReal();}
    for(auto & mps : MPS_R ){mps_real = mps_real and mps.isReal();}
    for(auto & mpo : MPO_L ){mpo_real = mpo_real and mpo->isReal();}
    for(auto & mpo : MPO_R ){mpo_real = mpo_real and mpo->isReal();}
    return mps_real and mpo_real;
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
    if(pos > get_length()-1 )  {throw std::range_error("get_A(pos): pos larger than length-1");}
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

Eigen::Tensor<class_finite_chain_state::Scalar,3> class_finite_chain_state::get_ENVL(size_t pos) const {
    if (pos >= ENV_L.size()){throw std::range_error("get_ENVL(pos) pos out of range: " + std::to_string(pos));}
    return std::next(ENV_L.begin(), pos)->block;
}

Eigen::Tensor<class_finite_chain_state::Scalar,3> class_finite_chain_state::get_ENVR(size_t pos) const {
    pos -= (ENV_L.size() + 1);
    if (pos >= ENV_R.size()){throw std::range_error("get_ENVR(pos) pos out of range: " + std::to_string(pos));}
    return std::next(ENV_R.begin(), pos)->block;
}

Eigen::Tensor<class_finite_chain_state::Scalar,4> class_finite_chain_state::get_ENV2L(size_t pos) const {
    if (pos >= ENV2_L.size()){throw std::range_error("get_ENV2L(pos) pos out of range: " + std::to_string(pos));}
    return std::next(ENV2_L.begin(), pos)->block;
}

Eigen::Tensor<class_finite_chain_state::Scalar,4> class_finite_chain_state::get_ENV2R(size_t pos) const {
    pos -= (ENV2_L.size() + 1);
    if (pos >= ENV2_R.size()){throw std::range_error("get_ENV2R(pos) pos out of range: " + std::to_string(pos));}
    return std::next(ENV2_R.begin(), pos)->block;
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

Eigen::DSizes<long,3> class_finite_chain_state::active_dimensions() const{
    Eigen::DSizes<long,3> dimensions;
    dimensions[1] = get_G(active_sites.front()).dimension(1);
    dimensions[2] = get_G(active_sites.back()).dimension(2);
    dimensions[0] = 1;
    for (auto & site : active_sites){
        dimensions[0] *= get_G(site).dimension(0);
    }
    return dimensions;
}




Eigen::Tensor<class_finite_chain_state::Scalar,3>   class_finite_chain_state::get_multitheta()    const{
    if(active_sites.empty()){throw std::runtime_error("No active sites on which to build multitheta");}
    Eigen::Tensor<Scalar,3> multitheta = get_A(active_sites[0]);
    Eigen::Tensor<Scalar,3> temp;
    for (size_t i = 1; i < active_sites.size(); i++ ){
        auto A = get_A(active_sites[i]);
        long dim0 = multitheta.dimension(0) * A.dimension(0);
        long dim1 = multitheta.dimension(1);
        long dim2 = A.dimension(2);
        temp = multitheta
               .contract(A, Textra::idx({2},{1}))
               .shuffle(Textra::array4{0,2,1,3})
               .reshape(Textra::array3{dim0,dim1,dim2});
        multitheta = temp;
    }
    auto & L = get_L(active_sites.back()+1);
    return multitheta.contract(Textra::asDiagonal(L), Textra::idx({2},{0}));
}

Eigen::Tensor<class_finite_chain_state::Scalar,4>   class_finite_chain_state::get_multimpo()    const{
    if(active_sites.empty()){throw std::runtime_error("No active sites on which to build multitheta");}
    Eigen::Tensor<Scalar,4> multimpo = get_MPO(active_sites[0]).MPO();
    Eigen::Tensor<Scalar,4> temp;
    for (size_t i = 1; i < active_sites.size(); i++ ){
        auto mpo  = get_MPO(active_sites[i]).MPO();
        long dim0 = multimpo.dimension(0);
        long dim1 = mpo.dimension(2);
        long dim2 = multimpo.dimension(2) * mpo.dimension(2);
        long dim3 = multimpo.dimension(3) * mpo.dimension(3);
        temp = multimpo
                .contract(mpo, Textra::idx({1},{0}))
                .shuffle(Textra::array6{0,3,1,4,2,5})
                .reshape(Textra::array4{dim0,dim1,dim2,dim3});
        multimpo = temp;
    }
    return multimpo;
}


std::pair<Eigen::Tensor<class_finite_chain_state::Scalar,3>, Eigen::Tensor<class_finite_chain_state::Scalar,3>>
class_finite_chain_state::get_multienv ()     const{
    return std::make_pair(get_ENVL(active_sites.back()), get_ENVR(active_sites.front()));
}
std::pair<Eigen::Tensor<class_finite_chain_state::Scalar,4>, Eigen::Tensor<class_finite_chain_state::Scalar,4>>
class_finite_chain_state::get_multienv2()     const{
    return std::make_pair(get_ENV2L(active_sites.back()), get_ENV2R(active_sites.front()));
}


Eigen::Tensor<class_finite_chain_state::Scalar,6>   class_finite_chain_state::get_multi_hamiltonian() const{
    auto [envL,envR] = get_multienv();
    return envL
           .contract(get_multimpo(), Textra::idx({2},{0}))
           .contract(envR          , Textra::idx({2},{2}))
           .shuffle(Textra::array6{0,2,4,1,3,5});
}

Eigen::Tensor<class_finite_chain_state::Scalar,6>   class_finite_chain_state::get_multi_hamiltonian2() const{
    auto [env2L,env2R] = get_multienv2();
    auto mpo = get_multimpo();
    return env2L
            .contract(mpo       , Textra::idx({2},{0}))
            .contract(mpo       , Textra::idx({5,2},{2,0}))
            .contract(env2R     , Textra::idx({2,4},{2,3}))
            .shuffle(Textra::array6{0,2,4,1,3,5});
}
Eigen::Matrix<class_finite_chain_state::Scalar,Eigen::Dynamic,Eigen::Dynamic> class_finite_chain_state::get_multi_hamiltonian_matrix() const{
    auto dims = active_dimensions();
    long shape = dims[0] * dims[1] * dims[2];
    return Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(get_multi_hamiltonian().data(),shape,shape);
}
Eigen::Matrix<class_finite_chain_state::Scalar,Eigen::Dynamic,Eigen::Dynamic> class_finite_chain_state::get_multi_hamiltonian2_matrix() const{
    auto dims = active_dimensions();
    long shape = dims[0] * dims[1] * dims[2];
    return Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(get_multi_hamiltonian2().data(),shape,shape);
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

