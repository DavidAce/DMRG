//
// Created by david on 2019-01-29.
//


#include "class_finite_state.h"
#include <tools/nmspc_tools.h>
#include <general/nmspc_quantum_mechanics.h>


void class_finite_state::do_all_measurements(){
    using namespace tools::finite;
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


void class_finite_state::set_positions(){

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

size_t class_finite_state::get_length()    const {return MPS_L.size() + MPS_R.size();}
size_t class_finite_state::get_position()  const {return MPS_L.size() - 1u;}

int  class_finite_state::get_sweeps()    const       {return num_sweeps;}
int  class_finite_state::reset_sweeps()              {num_sweeps = 0; return num_sweeps;}
void class_finite_state::set_sweeps(int num_sweeps_) {num_sweeps = num_sweeps_;}
void class_finite_state::increment_sweeps()          {num_sweeps++;}

int  class_finite_state::get_moves()    const       {return num_moves;}
int  class_finite_state::reset_moves()              { num_moves = 0; return num_moves;}
void class_finite_state::set_moves(int num_moves_)  { num_moves = num_moves_;}
void class_finite_state::increment_moves()          {num_moves++;}


long class_finite_state::get_chi_max()  const {return chi_max;}
void class_finite_state::set_chi_max(long chi_max_){ chi_max = chi_max_;}
int  class_finite_state::get_direction() const {return direction;}
void class_finite_state::flip_direction() {direction *= -1;}


Eigen::DSizes<long,3>
        class_finite_state::dimensions_2site()  const{
    Eigen::DSizes<long,3> dimensions;
    dimensions[1] = MPS_L.back().get_chiL();
    dimensions[2] = MPS_R.front().get_chiR();
    dimensions[0] = MPS_L.back().get_spin_dim() *  MPS_R.front().get_spin_dim();
    return dimensions;

}
size_t  class_finite_state::size_2site() const{
    auto dims = dimensions_2site();
    return dims[0]*dims[1]*dims[2];
}

bool class_finite_state::position_is_the_middle() const {
    return (size_t) get_position() + 1 == (size_t)(get_length() / 2.0) and direction == 1 ;
}
bool class_finite_state::position_is_the_middle_any_direction() const {
    return (size_t) get_position() + 1 == (size_t)(get_length() / 2.0);
}

bool class_finite_state::position_is_the_left_edge() const {
    return get_position() == 0;
}

bool class_finite_state::position_is_the_right_edge() const {
    return get_position() == get_length() - 2;
}

bool class_finite_state::position_is_any_edge() const {
    return position_is_the_left_edge() or position_is_the_right_edge();
}

bool class_finite_state::position_is_at(size_t pos)const{
    return get_position() == pos;
}

bool class_finite_state::isReal() const{
    bool mps_real = true;
    bool mpo_real = true;
    for(auto & mps : MPS_L ){mps_real = mps_real and mps.isReal();}
    for(auto & mps : MPS_R ){mps_real = mps_real and mps.isReal();}
    for(auto & mpo : MPO_L ){mpo_real = mpo_real and mpo->isReal();}
    for(auto & mpo : MPO_R ){mpo_real = mpo_real and mpo->isReal();}
    return mps_real and mpo_real;
}



Eigen::Tensor<class_finite_state::Scalar,3> class_finite_state::get_A() const{
    return Textra::asDiagonal(MPS_L.back().get_L()).contract(MPS_L.back().get_G(), Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
}

Eigen::Tensor<class_finite_state::Scalar,3> class_finite_state::get_B() const{
    return MPS_R.front().get_G().contract(Textra::asDiagonal(MPS_R.front().get_L()), Textra::idx({2},{0}));
}

Eigen::Tensor<class_finite_state::Scalar,4> class_finite_state::get_theta() const{
    return get_A()
           .contract(Textra::asDiagonal(MPS_C), Textra::idx({2},{0}))
           .contract(get_B(), Textra::idx({2},{1}));
}


const class_vidal_site & class_finite_state::get_MPS(size_t pos) const {
    if (pos >= get_length())                 throw std::range_error(fmt::format("get_MPS(pos) pos out of range: {}", pos));
    if(pos <= MPS_L.back().get_position()){
        auto mps_it = std::next(MPS_L.begin(),pos);
        if (mps_it->get_position() != pos)   throw std::range_error(fmt::format("get_MPS(pos): Mismatch in mps L position and pos: {} != {}", mps_it->get_position(), pos));
        return *mps_it;
    }else{
        if(pos < MPS_R.front().get_position()) throw std::range_error(fmt::format("get_MPS(pos): Mismatch in pos and MPSR front position: {} < {}", pos,  MPS_R.front().get_position()));
        auto mps_it = std::next(MPS_R.begin(), pos - MPS_R.front().get_position());
        if (mps_it->get_position() != pos)   throw std::range_error(fmt::format("get_MPS(pos): Mismatch in mps R position and pos: {} != {}", mps_it->get_position(), pos));
        return *mps_it;
    }
}

class_vidal_site & class_finite_state::get_MPS(size_t pos){
    return const_cast<class_vidal_site &>(static_cast<const class_finite_state &>(*this).get_MPS(pos));
}


const class_model_base & class_finite_state::get_MPO(size_t pos) const{
    if (pos >= get_length())throw std::range_error(fmt::format("get_MPO(pos) pos out of range: {}", pos));
    if(pos <= MPO_L.back()->get_position()){
        auto mpo_it = std::next(MPO_L.begin(),pos)->get();
        if (mpo_it->get_position() != pos)throw std::range_error(fmt::format("get_MPO(pos): Mismatch in mpo position and pos: {} != {}", mpo_it->get_position(), pos));
        return *mpo_it;
    }else{
        if(pos < MPO_R.front()->get_position()) throw std::range_error(fmt::format("get_MPS(pos): Mismatch in pos and MPOR front position: {} < {}", pos,  MPO_R.front()->get_position()));
        auto mpo_it = std::next(MPO_R.begin(), pos - MPO_R.front()->get_position())->get();
        if (mpo_it->get_position() != pos)throw std::range_error(fmt::format("get_MPO(pos): Mismatch in mpo position and pos: {} != {}", mpo_it->get_position(), pos));
        return *mpo_it;
    }
}

class_model_base & class_finite_state::get_MPO(size_t pos){
    return const_cast<class_model_base &>(static_cast<const class_finite_state &>(*this).get_MPO(pos));
}




const Eigen::Tensor<class_finite_state::Scalar,3> & class_finite_state::get_G(size_t pos) const{
    return std::as_const(get_MPS(pos).get_G());
}

Eigen::Tensor<class_finite_state::Scalar,3> & class_finite_state::get_G(size_t pos){
    return const_cast<Eigen::Tensor<class_finite_state::Scalar,3> &>(static_cast<const class_finite_state &>(*this).get_G(pos));
}

const Eigen::Tensor<class_finite_state::Scalar,1> & class_finite_state::get_L(size_t pos) const {
    if      (pos == MPS_L.back().get_position() + 1){return MPS_C;}
    else if (pos <= MPS_L.back().get_position())    {return get_MPS(pos).get_L();}
    else if (pos >= MPS_R.front().get_position())   {return get_MPS(pos-1).get_L();}
    else {throw std::runtime_error("Unhandled position");}
}

Eigen::Tensor<class_finite_state::Scalar,1> & class_finite_state::get_L(size_t pos) {
    return const_cast<Eigen::Tensor<class_finite_state::Scalar,1> &>(static_cast<const class_finite_state &>(*this).get_L(pos));
}




Eigen::Tensor<class_finite_state::Scalar,3> class_finite_state::get_A(size_t pos) const {
    return Textra::asDiagonal(get_L(pos)).contract(get_G(pos), Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
}

Eigen::Tensor<class_finite_state::Scalar,3> class_finite_state::get_B(size_t pos) const {
    return get_G(pos).contract(Textra::asDiagonal(get_L(pos+1)), Textra::idx({2},{0}));
}



const class_environment & class_finite_state::get_ENVL(size_t pos) const {
    if (pos > ENV_L.back().get_position() )  throw std::range_error(fmt::format("get_ENVL(pos):  pos is not in left side: {}",pos));
    if (pos >= ENV_L.size())                 throw std::range_error(fmt::format("get_ENVL(pos) pos out of range: {}",pos));
    auto env_it = std::next(ENV_L.begin(), pos);
    if (env_it->get_position() != pos) throw std::range_error(fmt::format("get_ENVL(pos): Mismatch in env position and pos: {} != {}", env_it->get_position(), pos));
    return *env_it;
}

const class_environment & class_finite_state::get_ENVR(size_t pos) const {
    if (pos < ENV_R.front().get_position() ){throw std::range_error(fmt::format("get_ENVR(pos):  pos is not in right side: {}" , pos));}
    if (pos >= get_length() )               {throw std::range_error(fmt::format("get_ENVR(pos):  pos out of range: {}" , pos));}

    if(pos < ENV_R.front().get_position()) throw std::range_error(fmt::format("get_ENVR(pos): Mismatch in pos and ENVR front position: {} < {}", pos,  ENV_R.front().get_position()));
    auto env_it = std::next(ENV_R.begin(), pos - ENV_R.front().get_position());
    if (env_it->get_position() != pos)      throw std::range_error(fmt::format("get_ENVR(pos): Mismatch in env position and pos: {} != {}", env_it->get_position(), pos));
    return *env_it;
}

const class_environment_var & class_finite_state::get_ENV2L(size_t pos) const {
    if (pos > ENV2_L.back().get_position() )  throw std::range_error(fmt::format("get_ENV2L(pos):  pos is not in left side: {}",pos));
    if (pos >= ENV2_L.size())                 throw std::range_error(fmt::format("get_ENV2L(pos) pos out of range: {}",pos));
    auto env2_it = std::next(ENV2_L.begin(), pos);
    if (env2_it->get_position() != pos)       throw std::range_error(fmt::format("get_ENV2L(pos): Mismatch in env position and pos: {} != {}", env2_it->get_position(), pos));
    return *env2_it;
}

const class_environment_var & class_finite_state::get_ENV2R(size_t pos) const {
    if (pos < ENV2_R.front().get_position() )throw std::range_error(fmt::format("get_ENV2R(pos):  pos is not in right side: {}" , pos));
    if (pos > ENV2_R.back().get_position() ) throw std::range_error(fmt::format("get_ENV2R(pos):  pos is not in right side: {}" , pos));
    if (pos >= get_length() )                throw std::range_error(fmt::format("get_ENV2R(pos):  pos out of range: {}" , pos));

    if(pos < ENV2_R.front().get_position()) throw std::range_error(fmt::format("get_ENV2R(pos): Mismatch in pos and ENV2R front position: {} < {}", pos,  ENV2_R.front().get_position()));
    auto env2_it = std::next(ENV2_R.begin(), pos - ENV2_R.front().get_position());
    if (env2_it->get_position() != pos)      throw std::range_error(fmt::format("get_ENV2R(pos): Mismatch in env2 position and pos: {} != {}", env2_it->get_position(), pos));
    return *env2_it;
}





Eigen::Tensor<class_finite_state::Scalar,4> class_finite_state::get_theta(size_t pos) const {
    return get_A(pos)
            .contract(Textra::asDiagonal(get_L(pos+1)), Textra::idx({2},{0}))
            .contract(get_B(pos+1), Textra::idx({2},{1}));
}


std::list<size_t> class_finite_state::activate_sites(long threshold){
    clear_cache();
    return active_sites = tools::finite::multisite::generate_site_list(*this,threshold);
}

Eigen::DSizes<long,3> class_finite_state::active_dimensions() const{
    Eigen::DSizes<long,3> dimensions;
    dimensions[1] = get_G(active_sites.front()).dimension(1);
    dimensions[2] = get_G(active_sites.back()).dimension(2);
    dimensions[0] = 1;
    for (auto & site : active_sites){
        dimensions[0] *= get_G(site).dimension(0);
    }
    return dimensions;
}

size_t class_finite_state::active_size() const {
    if (active_dimensions().empty()) return 0;
    auto dims = active_dimensions();
    return dims[0]*dims[1]*dims[2];
}


Eigen::Tensor<class_finite_state::Scalar,3>   class_finite_state::get_multitheta()    const{
    if(cache.multitheta) return cache.multitheta.value();
    tools::log->trace("Generating multitheta");
    if(active_sites.empty()){throw std::runtime_error("No active sites on which to build multitheta");}
    Eigen::Tensor<Scalar,3> multitheta;
    Eigen::Tensor<Scalar,3> temp;
    bool first = true;
    for (auto &site : active_sites){
        if (first){multitheta = get_A(site); first = false; continue;}
        auto A    = get_A(site);
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
    temp = multitheta.contract(Textra::asDiagonal(L), Textra::idx({2},{0}));
    cache.multitheta = temp;
    return cache.multitheta.value();
}

Eigen::Tensor<class_finite_state::Scalar,4>   class_finite_state::get_multimpo()    const{
    if(cache.multimpo) return cache.multimpo.value();
    tools::log->trace("Generating multimpo");
    if(active_sites.empty()){throw std::runtime_error("No active sites on which to build multimpo");}
    Eigen::Tensor<Scalar,4> multimpo;
    Eigen::Tensor<Scalar,4> temp;
    bool first = true;
    for (auto &site : active_sites){
        if (first){multimpo = get_MPO(site).MPO(); first = false; continue;}
        auto &mpo = get_MPO(site).MPO();
        long dim0 = multimpo.dimension(0);
        long dim1 = mpo.dimension(1);
        long dim2 = multimpo.dimension(2) * mpo.dimension(2);
        long dim3 = multimpo.dimension(3) * mpo.dimension(3);
        temp = multimpo
                .contract(mpo, Textra::idx({1},{0}))
                .shuffle(Textra::array6{0,3,1,4,2,5})
                .reshape(Textra::array4{dim0,dim1,dim2,dim3});
        multimpo = temp;
    }
    cache.multimpo = multimpo;
    return cache.multimpo.value();
}


std::pair<std::reference_wrapper<const class_environment> , std::reference_wrapper<const class_environment>>
class_finite_state::get_multienv ()     const{
    return std::make_pair(get_ENVL(active_sites.front()), get_ENVR(active_sites.back()));
}

std::pair<std::reference_wrapper<const class_environment_var> , std::reference_wrapper<const class_environment_var>>
class_finite_state::get_multienv2()     const{
    return std::make_pair(get_ENV2L(active_sites.front()), get_ENV2R(active_sites.back()));
}


class_finite_state::TType<6> class_finite_state::get_multi_hamiltonian() const{
    tools::log->trace("Generating multi hamiltonian");

//    auto [envL,envR] = get_multienv();
    auto & envL = get_ENVL(active_sites.front());
    auto & envR = get_ENVR(active_sites.back());
    if (envL.get_position() != active_sites.front()) throw std::runtime_error(fmt::format("Mismatch in ENVL and active site positions: {} != {}", envL.get_position() , active_sites.front()));
    if (envR.get_position() != active_sites.back())  throw std::runtime_error(fmt::format("Mismatch in ENVR and active site positions: {} != {}", envR.get_position() , active_sites.back()));
    return envL.block
           .contract(get_multimpo(), Textra::idx({2},{0}))
           .contract(envR.block    , Textra::idx({2},{2}))
            .shuffle(Textra::array6{2,0,4,3,1,5});
}

class_finite_state::TType<6>   class_finite_state::get_multi_hamiltonian2() const{
    tools::log->trace("Generating multi hamiltonian squared");
//    auto [env2L,env2R] = get_multienv2();
    auto & env2L = get_ENV2L(active_sites.front());
    auto & env2R = get_ENV2R(active_sites.back());
    if (env2L.get_position() != active_sites.front()) throw std::runtime_error(fmt::format("Mismatch in ENVL and active site positions: {} != {}", env2L.get_position() , active_sites.front()));
    if (env2R.get_position() != active_sites.back())  throw std::runtime_error(fmt::format("Mismatch in ENVR and active site positions: {} != {}", env2R.get_position() , active_sites.back()));
    auto mpo = get_multimpo();
    tools::log->trace("Starting contraction of multi hamiltonian squared");
    return env2L.block
            .contract(mpo             , Textra::idx({2},{0}))
            .contract(mpo             , Textra::idx({5,2},{2,0}))
            .contract(env2R.block     , Textra::idx({2,4},{2,3}))
            .shuffle(Textra::array6{2,0,4,3,1,5});
}

class_finite_state::MType class_finite_state::get_multi_hamiltonian_matrix() const{
    auto dims = active_dimensions();
    long shape = dims[0] * dims[1] * dims[2];
    return Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(get_multi_hamiltonian().data(),shape,shape).transpose();
}
class_finite_state::MType class_finite_state::get_multi_hamiltonian2_matrix() const{
    auto dims = active_dimensions();
    long shape = dims[0] * dims[1] * dims[2];
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> ham_squared = Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(get_multi_hamiltonian2().data(),shape,shape);
    tools::log->trace("Finished contraction of multi hamiltonian squared");
    return ham_squared.transpose();
}

class_finite_state::MType  class_finite_state::get_multi_hamiltonian2_subspace_matrix(const MType & eigvecs ) const{
    tools::log->trace("Generating multi hamiltonian squared");

    auto dims = active_dimensions();
    Eigen::DSizes<long,4> eigvecs_dims {dims[0],dims[1],dims[2],eigvecs.cols()};
    auto eigvecs_tensor = Eigen::TensorMap<const Eigen::Tensor<const Scalar,4>>(eigvecs.data(), eigvecs_dims );
    auto & env2L = get_ENV2L(active_sites.front());
    auto & env2R = get_ENV2R(active_sites.back());
    if (env2L.get_position() != active_sites.front()) throw std::runtime_error(fmt::format("Mismatch in ENVL and active site positions: {} != {}", env2L.get_position() , active_sites.front()));
    if (env2R.get_position() != active_sites.back())  throw std::runtime_error(fmt::format("Mismatch in ENVR and active site positions: {} != {}", env2R.get_position() , active_sites.back()));
    auto mpo = get_multimpo();
    tools::log->trace("Starting contraction of hamiltonian subspace matrix squared");
    size_t log2chiL  = std::log2(dims[1]);
    size_t log2chiR  = std::log2(dims[2]);
    size_t log2spin  = std::log2(dims[0]);
    Eigen::Tensor<Scalar,2> H2;
    if(log2spin > log2chiL + log2chiR){
        if (log2chiL >= log2chiR){
            tools::log->trace("get_H2 path: log2spin > log2chiL + log2chiR  and  log2chiL >= log2chiR ");
            H2 =
                    eigvecs_tensor
                            .contract(env2L.block,                Textra::idx({1},{0}))
                            .contract(mpo  ,                      Textra::idx({0,4},{2,0}))
                            .contract(env2R.block,                Textra::idx({0,4},{0,2}))
                            .contract(mpo  ,                      Textra::idx({3,2,5},{2,0,1}))
                            .contract(eigvecs_tensor.conjugate(), Textra::idx({3,1,2},{0,1,2}));
        }
        else{
            tools::log->trace("get_H2 path: log2spin > log2chiL + log2chiR  and  log2chiL < log2chiR ");
            H2 =
                    eigvecs_tensor
                            .contract(env2R.block,                Textra::idx({2},{0}))
                            .contract(mpo  ,                      Textra::idx({0,4},{2,1}))
                            .contract(env2L.block,                Textra::idx({0,4},{0,2}))
                            .contract(mpo  ,                      Textra::idx({3,2,5},{2,1,0}))
                            .contract(eigvecs_tensor.conjugate(), Textra::idx({3,2,1},{0,1,2}));
        }
    }else{
        tools::log->trace("get_H2 path: log2spin <= log2chiL + log2chiR");

        H2 =
                eigvecs_tensor
                        .contract(env2L.block,                Textra::idx({1},{0}))
                        .contract(mpo  ,                      Textra::idx({0,4},{2,0}))
                        .contract(mpo  ,                      Textra::idx({5,3},{2,0}))
                        .contract(eigvecs_tensor.conjugate(), Textra::idx({5,2},{0,1}))
                        .contract(env2R.block,                Textra::idx({0,4,2,3},{0,1,2,3}));
    }
//    H2 =
//            eigvecs_tensor
//                    .contract(env2L.block,                Textra::idx({1},{0}))
//                    .contract(mpo  ,                      Textra::idx({0,4},{2,0}))
//                    .contract(mpo  ,                      Textra::idx({5,3},{2,0}))
//                    .contract(eigvecs_tensor.conjugate(), Textra::idx({5,2},{0,1}))
//                    .contract(env2R.block,                Textra::idx({0,4,2,3},{0,1,2,3}));
    tools::log->trace("Finished contraction of hamiltonian subspace matrix squared");
    return Eigen::Map<MType>(H2.data(),H2.dimension(0),H2.dimension(1));

}


void class_finite_state::unset_measurements()const {
    measurements = Measurements();
}

void class_finite_state::clear_cache()const {
    cache = Cache();
}

void class_finite_state::do_all_measurements()const {
    measurements.length                           = tools::finite::measure::length                        (*this);
    measurements.bond_dimension_current           = tools::finite::measure::bond_dimension_current        (*this);
    measurements.bond_dimension_midchain          = tools::finite::measure::bond_dimension_midchain       (*this);
    measurements.bond_dimensions                  = tools::finite::measure::bond_dimensions               (*this);
    measurements.norm                             = tools::finite::measure::norm                          (*this);
    measurements.energy                           = tools::finite::measure::energy                        (*this);
    measurements.energy_per_site                  = tools::finite::measure::energy_per_site               (*this);
    measurements.energy_variance_mpo              = tools::finite::measure::energy_variance               (*this);
    measurements.energy_variance_per_site         = tools::finite::measure::energy_variance_per_site      (*this);
    measurements.spin_components                  = tools::finite::measure::spin_components               (*this); // This will automatically measure sx,sy and sz as well
    measurements.entanglement_entropy_current     = tools::finite::measure::entanglement_entropy_current  (*this);
    measurements.entanglement_entropy_midchain    = tools::finite::measure::entanglement_entropy_midchain (*this);
    measurements.entanglement_entropies           = tools::finite::measure::entanglement_entropies        (*this);
}


void class_finite_state::tag_active_sites_have_been_updated(bool tag)     const{
    if (site_update_tags.size() != get_length()) throw std::runtime_error("Cannot tag active sites, size mismatch in site list");
    for (auto & site: active_sites){
        site_update_tags[site] = tag;
    }
}

void class_finite_state::tag_all_sites_have_been_updated(bool tag) const{
    if (site_update_tags.size() != get_length()) throw std::runtime_error("Cannot untag all sites, size mismatch in site list");
    site_update_tags = std::vector<bool>(get_length(),tag);
}

bool class_finite_state::all_sites_updated() const {
    if (site_update_tags.size() != get_length()) throw std::runtime_error("Cannot check update status on all sites, size mismatch in site list");
    return  std::all_of(site_update_tags.begin(), site_update_tags.end(), [](bool v) { return v; });
}

