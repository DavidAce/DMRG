//
// Created by david on 2019-01-29.
//

#include "class_state_finite.h"
#include <tools/finite/measure.h>
#include <tools/finite/mpo.h>
#include <tools/finite/mps.h>
#include <tools/finite/multisite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
// We need to make a destructor manually for the enclosing class "class_state_finite"
// that encloses "class_model_base". Otherwise unique_ptr will forcibly inline its
// own default deleter.
// This allows us to forward declare the abstract base class "class_model_base"
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_state_finite::~class_state_finite() = default;

class_state_finite::class_state_finite(const class_state_finite &other) { *this = other; }

class_state_finite &class_state_finite::operator=(const class_state_finite &other) {
    // check for self-assignment
    if(&other == this) return *this;

    // Copy all data members
    this->num_sweeps = other.num_sweeps;
    this->num_moves  = other.num_moves;
    this->direction  = other.direction;
    this->chi_lim    = other.chi_lim;
    this->chi_max    = other.chi_max;
    this->MPS_L      = other.MPS_L;
    this->MPS_R      = other.MPS_R;
    this->ENV_L      = other.ENV_L;
    this->ENV_R      = other.ENV_R;
    this->ENV2_L     = other.ENV2_L;
    this->ENV2_R     = other.ENV2_R;

    this->active_sites     = other.active_sites;
    this->truncation_error = other.truncation_error;
    this->truncated_variance = other.truncated_variance;
    this->measurements     = other.measurements;
    this->site_update_tags = other.site_update_tags;
    this->cache            = other.cache;

    // The MPO's are special and the whole point of doing this manually
    this->MPO_L.clear();
    this->MPO_R.clear();
    for(auto &mpo : other.MPO_L) this->MPO_L.emplace_back(mpo->clone());
    for(auto &mpo : other.MPO_R) this->MPO_R.emplace_back(mpo->clone());
    return *this;
}

void class_state_finite::do_all_measurements() {
    using namespace tools::finite;
    measurements.length                        = measure::length(*this);
    measurements.bond_dimension_current        = measure::bond_dimension_current(*this);
    measurements.bond_dimension_midchain       = measure::bond_dimension_midchain(*this);
    measurements.bond_dimensions               = measure::bond_dimensions(*this);
    measurements.norm                          = measure::norm(*this);
    measurements.energy                        = measure::energy(*this); // This number is needed for variance calculation!
    measurements.energy_per_site               = measure::energy_per_site(*this);
    measurements.energy_variance               = measure::energy_variance(*this);
    measurements.energy_variance_per_site      = measure::energy_variance_per_site(*this);
    measurements.entanglement_entropy_current  = measure::entanglement_entropy_current(*this);
    measurements.entanglement_entropy_midchain = measure::entanglement_entropy_midchain(*this);
    measurements.entanglement_entropies        = measure::entanglement_entropies(*this);
    measurements.spin_components               = measure::spin_components(*this);
}

void class_state_finite::set_positions() {
    size_t pos = 0;
    for(auto &MPS : MPS_L) {
        MPS.set_position(pos++);
    }
    for(auto &MPS : MPS_R) {
        MPS.set_position(pos++);
    }
    pos = 0;
    for(auto &ENV : ENV_L) {
        ENV.set_position(pos++);
    }
    for(auto &ENV : ENV_R) {
        ENV.set_position(pos++);
    }
    pos = 0;
    for(auto &ENV2 : ENV2_L) {
        ENV2.set_position(pos++);
    }
    for(auto &ENV2 : ENV2_R) {
        ENV2.set_position(pos++);
    }
    pos = 0;
    for(auto &MPO : MPO_L) {
        MPO->set_position(pos++);
    }
    for(auto &MPO : MPO_R) {
        MPO->set_position(pos++);
    }
}

size_t class_state_finite::get_length() const { return MPS_L.size() + MPS_R.size(); }
size_t class_state_finite::get_position() const { return MPS_L.size() - 1u; }

size_t class_state_finite::get_sweeps() const { return num_sweeps; }
size_t class_state_finite::reset_sweeps() {
    num_sweeps = 0;
    return num_sweeps;
}
void class_state_finite::set_sweeps(int num_sweeps_) { num_sweeps = num_sweeps_; }
void class_state_finite::increment_sweeps() { num_sweeps++; }

size_t class_state_finite::get_moves() const { return num_moves; }
size_t class_state_finite::reset_moves() {
    num_moves = 0;
    return num_moves;
}
void class_state_finite::set_moves(int num_moves_) { num_moves = num_moves_; }
void class_state_finite::increment_moves() { num_moves++; }

long class_state_finite::get_chi_lim() const {
    // Should get the the current limit on allowed bond dimension
    return chi_lim.value();
}
void class_state_finite::set_chi_lim(long chi_lim_) {
    // Should set the the current limit on allowed bond dimension
    if(chi_lim_ == 0) throw std::runtime_error("Can't set chi limit to zero!");
    chi_lim = chi_lim_;
}

long class_state_finite::get_chi_max() const {
    // Should get the the current limit on allowed bond dimension for the duration of the simulation
    return chi_max.value();
}
void class_state_finite::set_chi_max(long chi_max_) {
    // Should set the the highest limit on allowed bond dimension for the duration of the simulation
    if(chi_max_ == 0) throw std::runtime_error("Can't set chi max to zero!");
    chi_max = chi_max_;
}
long class_state_finite::find_largest_chi() const {
    auto bond_dimensions    = tools::finite::measure::bond_dimensions(*this);
    return *max_element(std::begin(bond_dimensions), std::end(bond_dimensions));
}
int  class_state_finite::get_direction() const { return direction; }
void class_state_finite::flip_direction() { direction *= -1; }

Eigen::DSizes<long, 3> class_state_finite::dimensions_2site() const {
    Eigen::DSizes<long, 3> dimensions;
    dimensions[1] = MPS_L.back().get_chiL();
    dimensions[2] = MPS_R.front().get_chiR();
    dimensions[0] = MPS_L.back().get_spin_dim() * MPS_R.front().get_spin_dim();
    return dimensions;
}
size_t class_state_finite::size_2site() const {
    auto dims = dimensions_2site();
    return dims[0] * dims[1] * dims[2];
}

bool class_state_finite::position_is_the_middle() const { return (size_t) get_position() + 1 == (size_t)(get_length() / 2.0) and direction == 1; }
bool class_state_finite::position_is_the_middle_any_direction() const { return (size_t) get_position() + 1 == (size_t)(get_length() / 2.0); }

bool class_state_finite::position_is_left_edge() const { return get_position() == 0 and direction == -1; }

bool class_state_finite::position_is_right_edge() const { return get_position() == get_length() - 2 and direction == 1; }

bool class_state_finite::position_is_any_edge() const { return position_is_left_edge() or position_is_right_edge(); }

bool class_state_finite::position_is_at(size_t pos) const { return get_position() == pos; }

bool class_state_finite::isReal() const {
    bool mps_real = true;
    bool mpo_real = true;
    bool env_real = true;
    for(auto &mps : MPS_L) {
        mps_real = mps_real and mps.isReal();
    }
    for(auto &mps : MPS_R) {
        mps_real = mps_real and mps.isReal();
    }
    for(auto &mpo : MPO_L) {
        mpo_real = mpo_real and mpo->isReal();
    }
    for(auto &mpo : MPO_R) {
        mpo_real = mpo_real and mpo->isReal();
    }
    for(auto &env : ENV_L) {
        env_real = env_real and env.isReal();
    }
    for(auto &env : ENV_R) {
        env_real = env_real and env.isReal();
    }
    for(auto &env : ENV2_R) {
        env_real = env_real and env.isReal();
    }
    for(auto &env : ENV2_L) {
        env_real = env_real and env.isReal();
    }
    return mps_real and mpo_real and env_real;
}

bool class_state_finite::hasNaN() const {
    for(auto &mps : MPS_L)
        if(mps.hasNaN()) return true;
    for(auto &mps : MPS_R)
        if(mps.hasNaN()) return true;
    for(auto &mpo : MPO_L)
        if(mpo->hasNaN()) return true;
    for(auto &mpo : MPO_R)
        if(mpo->hasNaN()) return true;
    for(auto &env : ENV_L)
        if(env.hasNaN()) return true;
    for(auto &env : ENV_R)
        if(env.hasNaN()) return true;
    for(auto &env : ENV2_L)
        if(env.hasNaN()) return true;
    for(auto &env : ENV2_R)
        if(env.hasNaN()) return true;
    return false;
}

void class_state_finite::assertValidity() const {
    for(auto &mps : MPS_L) mps.assertValidity();
    for(auto &mps : MPS_R) mps.assertValidity();
    for(auto &mpo : MPO_L) mpo->assertValidity();
    for(auto &mpo : MPO_R) mpo->assertValidity();
    for(auto &env : ENV_L) env.assertValidity();
    for(auto &env : ENV_R) env.assertValidity();
    for(auto &env : ENV2_L) env.assertValidity();
    for(auto &env : ENV2_R) env.assertValidity();
}

const Eigen::Tensor<class_state_finite::Scalar, 1> &class_state_finite::midchain_bond() const {
    size_t center_pos = (get_length() - 1) / 2;
    if(get_position() < center_pos) return get_MPS(center_pos).get_L();
    if(get_position() > center_pos) return get_MPS(center_pos + 1).get_L();
    if(get_position() == center_pos)
        return get_MPS(center_pos).get_LC();
    else
        throw std::logic_error("No valid position to find midchain_bond");
}

const Eigen::Tensor<class_state_finite::Scalar, 1> &class_state_finite::current_bond() const {
    if(MPS_L.back().isCenter()) {
        if(MPS_L.back().get_chiR() != MPS_L.back().get_LC().dimension(0)) throw std::runtime_error("Current center bond dimension mismatch!");
        return MPS_L.back().get_LC();
    } else {
        throw std::runtime_error("Current center bond LC hasn't been set.");
    }
}

// Eigen::Tensor<class_state_finite::Scalar,3> class_state_finite::get_A() const{
//    return Textra::asDiagonal(MPS_L.back().get_L()).contract(MPS_L.back().get_G(), Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
//}
//
// Eigen::Tensor<class_state_finite::Scalar,3> class_state_finite::get_B() const{
//    return MPS_R.front().get_G().contract(Textra::asDiagonal(MPS_R.front().get_L()), Textra::idx({2},{0}));
//}
//
Eigen::Tensor<class_state_finite::Scalar, 4> class_state_finite::get_theta() const {
    if(cache.theta) return cache.theta.value();
    tools::log->trace("Contracting theta...");
    cache.theta = MPS_L.back().get_M().contract(MPS_R.front().get_M(), Textra::idx({2}, {1}));
    tools::log->trace("Contracting theta... OK");
    return cache.theta.value();
}

const class_mps_site &class_state_finite::get_MPS(size_t pos) const {
    if(pos >= get_length()) throw std::range_error(fmt::format("get_MPS(pos) pos out of range: {}", pos));
    if(pos <= MPS_L.back().get_position()) {
        auto mps_it = std::next(MPS_L.begin(), pos);
        if(mps_it->get_position() != pos)
            throw std::range_error(fmt::format("get_MPS(pos): Mismatch in mps L position and pos: {} != {}", mps_it->get_position(), pos));
        return *mps_it;
    } else {
        if(pos < MPS_R.front().get_position())
            throw std::range_error(fmt::format("get_MPS(pos): Mismatch in pos and MPSR front position: {} < {}", pos, MPS_R.front().get_position()));
        auto mps_it = std::next(MPS_R.begin(), pos - MPS_R.front().get_position());
        if(mps_it->get_position() != pos)
            throw std::range_error(fmt::format("get_MPS(pos): Mismatch in mps R position and pos: {} != {}", mps_it->get_position(), pos));
        return *mps_it;
    }
}

class_mps_site &class_state_finite::get_MPS(size_t pos) { return const_cast<class_mps_site &>(static_cast<const class_state_finite &>(*this).get_MPS(pos)); }

const class_model_base &class_state_finite::get_MPO(size_t pos) const {
    if(pos >= get_length()) throw std::range_error(fmt::format("get_MPO(pos) pos out of range: {}", pos));
    if(pos <= MPO_L.back()->get_position()) {
        auto mpo_it = std::next(MPO_L.begin(), pos)->get();
        if(mpo_it->get_position() != pos)
            throw std::range_error(fmt::format("get_MPO(pos): Mismatch in mpo position and pos: {} != {}", mpo_it->get_position(), pos));
        return *mpo_it;
    } else {
        if(pos < MPO_R.front()->get_position())
            throw std::range_error(fmt::format("get_MPS(pos): Mismatch in pos and MPOR front position: {} < {}", pos, MPO_R.front()->get_position()));
        auto mpo_it = std::next(MPO_R.begin(), pos - MPO_R.front()->get_position())->get();
        if(mpo_it->get_position() != pos)
            throw std::range_error(fmt::format("get_MPO(pos): Mismatch in mpo position and pos: {} != {}", mpo_it->get_position(), pos));
        return *mpo_it;
    }
}

class_model_base &class_state_finite::get_MPO(size_t pos) {
    return const_cast<class_model_base &>(static_cast<const class_state_finite &>(*this).get_MPO(pos));
}

const class_environment &class_state_finite::get_ENVL(size_t pos) const {
    if(pos > ENV_L.back().get_position()) throw std::range_error(fmt::format("get_ENVL(pos):  pos is not in left side: {}", pos));
    if(pos >= ENV_L.size()) throw std::range_error(fmt::format("get_ENVL(pos) pos out of range: {}", pos));
    auto env_it = std::next(ENV_L.begin(), pos);
    if(env_it->get_position() != pos)
        throw std::range_error(fmt::format("get_ENVL(pos): Mismatch in env position and pos: {} != {}", env_it->get_position(), pos));
    return *env_it;
}

const class_environment &class_state_finite::get_ENVR(size_t pos) const {
    if(pos < ENV_R.front().get_position()) {
        throw std::range_error(fmt::format("get_ENVR(pos):  pos is not in right side: {}", pos));
    }
    if(pos >= get_length()) {
        throw std::range_error(fmt::format("get_ENVR(pos):  pos out of range: {}", pos));
    }
    if(pos < ENV_R.front().get_position())
        throw std::range_error(fmt::format("get_ENVR(pos): Mismatch in pos and ENVR front position: {} < {}", pos, ENV_R.front().get_position()));
    auto env_it = std::next(ENV_R.begin(), pos - ENV_R.front().get_position());
    if(env_it->get_position() != pos)
        throw std::range_error(fmt::format("get_ENVR(pos): Mismatch in env position and pos: {} != {}", env_it->get_position(), pos));
    return *env_it;
}

const class_environment_var &class_state_finite::get_ENV2L(size_t pos) const {
    if(pos > ENV2_L.back().get_position()) throw std::range_error(fmt::format("get_ENV2L(pos):  pos is not in left side: {}", pos));
    if(pos >= ENV2_L.size()) throw std::range_error(fmt::format("get_ENV2L(pos) pos out of range: {}", pos));
    auto env2_it = std::next(ENV2_L.begin(), pos);
    if(env2_it->get_position() != pos)
        throw std::range_error(fmt::format("get_ENV2L(pos): Mismatch in env position and pos: {} != {}", env2_it->get_position(), pos));
    return *env2_it;
}

const class_environment_var &class_state_finite::get_ENV2R(size_t pos) const {
    if(pos < ENV2_R.front().get_position()) throw std::range_error(fmt::format("get_ENV2R(pos):  pos is not in right side: {}", pos));
    if(pos > ENV2_R.back().get_position()) throw std::range_error(fmt::format("get_ENV2R(pos):  pos is not in right side: {}", pos));
    if(pos >= get_length()) throw std::range_error(fmt::format("get_ENV2R(pos):  pos out of range: {}", pos));
    if(pos < ENV2_R.front().get_position())
        throw std::range_error(fmt::format("get_ENV2R(pos): Mismatch in pos and ENV2R front position: {} < {}", pos, ENV2_R.front().get_position()));
    auto env2_it = std::next(ENV2_R.begin(), pos - ENV2_R.front().get_position());
    if(env2_it->get_position() != pos)
        throw std::range_error(fmt::format("get_ENV2R(pos): Mismatch in env2 position and pos: {} != {}", env2_it->get_position(), pos));
    return *env2_it;
}

// For reduced energy MPO's

bool class_state_finite::isReduced() const {
    bool reduced = MPO_L.front()->is_reduced();
    for(auto &mpo : MPO_L)
        if(reduced != mpo->is_reduced()) {
            throw std::runtime_error(
                fmt::format("First MPO has is_reduced: {}, but MPO at pos {} has is_reduced: {}", reduced, mpo->get_position(), mpo->is_reduced()));
        }
    for(auto &mpo : MPO_R)
        if(reduced != mpo->is_reduced()) {
            throw std::runtime_error(
                fmt::format("First MPO has is_reduced: {}, but MPO at pos {} has is_reduced: {}", reduced, mpo->get_position(), mpo->is_reduced()));
        }
    return reduced;
}

double class_state_finite::get_energy_reduced() const { return get_energy_per_site_reduced() * get_length(); }

double class_state_finite::get_energy_per_site_reduced() const {
    // Check that all energies are the same
    double e_reduced = MPO_L.front()->get_reduced_energy();
    for(auto &mpo : MPO_L) {
        mpo->get_reduced_energy();
        if(mpo->get_reduced_energy() != e_reduced) {
            throw std::runtime_error("Reduced energy mismatch!");
        }
    }
    for(auto &mpo : MPO_R) {
        mpo->get_reduced_energy();
        if(mpo->get_reduced_energy() != e_reduced) {
            throw std::runtime_error("Reduced energy mismatch!");
        }
    }
    return e_reduced;
}

void class_state_finite::set_reduced_energy(double total_energy) { set_reduced_energy_per_site(total_energy / get_length()); }

void class_state_finite::set_reduced_energy_per_site(double site_energy) {
    if(get_energy_per_site_reduced() == site_energy) return;
    clear_measurements();
    cache.multimpo = {};
    for(auto &mpo : MPO_L) mpo->set_reduced_energy(site_energy);
    for(auto &mpo : MPO_R) mpo->set_reduced_energy(site_energy);
    tools::finite::mps::rebuild_environments(*this);
}

void class_state_finite::perturb_hamiltonian(double coupling_ptb, double field_ptb, PerturbMode perturbMode) {
    tools::finite::mpo::perturb_hamiltonian(*this, coupling_ptb, field_ptb, perturbMode);
}

void class_state_finite::damp_hamiltonian(double coupling_damp, double field_damp){
    tools::finite::mpo::damp_hamiltonian(*this, coupling_damp, field_damp);
}



bool class_state_finite::is_perturbed() const {
    for(size_t pos = 0; pos < get_length(); pos++) {
        if(get_MPO(pos).is_perturbed()) return true;
    }
    return false;
}

bool class_state_finite::is_damped() const {
    for(size_t pos = 0; pos < get_length(); pos++) {
        if(get_MPO(pos).is_damped()) return true;
    }
    return false;
}


std::list<size_t> class_state_finite::activate_sites(const long threshold, const size_t max_sites, const size_t min_sites) {
    clear_cache();
    return active_sites = tools::finite::multisite::generate_site_list(*this, threshold, max_sites, min_sites);
}

std::list<size_t> class_state_finite::activate_truncated_sites(const long threshold, const size_t chi_lim, const size_t max_sites, const size_t min_sites) {
    clear_cache();
    return active_sites = tools::finite::multisite::generate_truncated_site_list(*this, threshold,chi_lim, max_sites, min_sites);
}

Eigen::DSizes<long, 3> class_state_finite::active_dimensions() const { return tools::finite::multisite::get_dimensions(*this, active_sites); }

size_t class_state_finite::active_problem_size() const { return tools::finite::multisite::get_problem_size(*this, active_sites); }

const Eigen::Tensor<class_state_finite::Scalar, 3> &class_state_finite::get_multitheta() const {
    if(cache.multitheta) return cache.multitheta.value();
    tools::log->trace("Contracting multi theta...");
    if(active_sites.empty()) {
        throw std::runtime_error("No active sites on which to build multitheta");
    }
    Eigen::Tensor<Scalar, 3> multitheta;
    Eigen::Tensor<Scalar, 3> temp;
    bool                     first = true;
    for(auto &site : active_sites) {
        if(first) {
            multitheta = get_MPS(site).get_M();
            first      = false;
            continue;
        }
        if(site != get_MPS(site).get_position()) throw std::runtime_error("Site mismatch in get_multitheta");
        auto M     = get_MPS(site).get_M();
        long dim0  = multitheta.dimension(0) * M.dimension(0);
        long dim1  = multitheta.dimension(1);
        long dim2  = M.dimension(2);
        temp       = multitheta.contract(M, Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(Textra::array3{dim0, dim1, dim2});
        multitheta = temp;
    }
    //    auto & L = get_L(active_sites.back()+1);
    //    temp = multitheta.contract(Textra::asDiagonal(L), Textra::idx({2},{0}));
    tools::log->trace("Contracting multi theta... OK");
    cache.multitheta = temp;
    return cache.multitheta.value();
}

const Eigen::Tensor<class_state_finite::Scalar, 4> &class_state_finite::get_multimpo() const {
    if(cache.multimpo) return cache.multimpo.value();
    tools::common::profile::t_mpo->tic();
    tools::log->trace("Contracting multi mpo...");
    if(active_sites.empty()) {
        throw std::runtime_error("No active sites on which to build multimpo");
    }
    Eigen::Tensor<Scalar, 4> multimpo;
    bool                     first = true;
    for(auto &site : active_sites) {
        if(first) {
            multimpo = get_MPO(site).MPO();
            first    = false;
            continue;
        }
        auto &                   mpo  = get_MPO(site).MPO();
        long                     dim0 = multimpo.dimension(0);
        long                     dim1 = mpo.dimension(1);
        long                     dim2 = multimpo.dimension(2) * mpo.dimension(2);
        long                     dim3 = multimpo.dimension(3) * mpo.dimension(3);
        Eigen::Tensor<Scalar, 4> temp(dim0, dim1, dim2, dim3);
        temp     = multimpo.contract(mpo, Textra::idx({1}, {0})).shuffle(Textra::array6{0, 3, 1, 4, 2, 5}).reshape(Textra::array4{dim0, dim1, dim2, dim3});
        multimpo = temp;
    }
    tools::log->trace("Contracting multi mpo... OK");
    tools::common::profile::t_mpo->toc();
    cache.multimpo = multimpo;
    return cache.multimpo.value();
}

std::pair<std::reference_wrapper<const class_environment>, std::reference_wrapper<const class_environment>> class_state_finite::get_multienv() const {
    return std::make_pair(get_ENVL(active_sites.front()), get_ENVR(active_sites.back()));
}

std::pair<std::reference_wrapper<const class_environment_var>, std::reference_wrapper<const class_environment_var>> class_state_finite::get_multienv2() const {
    return std::make_pair(get_ENV2L(active_sites.front()), get_ENV2R(active_sites.back()));
}

void class_state_finite::set_truncation_error(size_t left_site, double error)
/*! The truncation error vector has length + 1 elements, as many as there are bond matrices.
 *  Obviously, the edge matrices are always 1, so we expect these to have truncation_error = 0.
 *  Any schmidt decomposition involves two neighboriing sites, so "left_site" is the site number of the
 *  site on the left of the decomposition.
 */
{
    if(truncation_error.empty()) {
        tools::log->debug("Resizing truncation_error container to size: {}", get_length() + 1);
        truncation_error = std::vector<double>(get_length() + 1, 0.0);
    }
    truncation_error[left_site + 1] = error;
}

void class_state_finite::set_truncation_error(double error) { set_truncation_error(get_position(), error); }

double class_state_finite::get_truncation_error(size_t left_site) const {
    if(truncation_error.empty() or truncation_error.size() < get_length() + 1) {
        tools::log->warn("Truncation error container hasn't been initialized! You called this function prematurely");
        return std::numeric_limits<double>::infinity();
    }
    return truncation_error[left_site + 1];
}

double class_state_finite::get_truncation_error() const { return get_truncation_error(get_position()); }
double class_state_finite::get_truncation_error_midchain() const {
    size_t center_pos = (get_length() - 1) / 2;
    return get_truncation_error(center_pos);}

const std::vector<double> &class_state_finite::get_truncation_errors() const { return truncation_error; }

void class_state_finite::set_truncated_variance(size_t left_site, double error)
{
/*! The truncated variance vector has length + 1 elements, as many as there are bond matrices.
*  Obviously, the edge matrices are always 1, so we expect these to have truncated variance = 0.
*  Any schmidt decomposition involves two neighboriing sites, so "left_site" is the site number of the
*  site on the left of the decomposition.
*/
        if(truncated_variance.empty()) {
            tools::log->debug("Resizing truncated_variance container to size: {}", get_length() + 1);
            truncated_variance = std::vector<double>(get_length() + 1, 0.0);
        }
        truncated_variance[left_site + 1] = error;
}
void class_state_finite::set_truncated_variance(double error){
    set_truncated_variance(get_position(), error);
}

double class_state_finite::get_truncated_variance(size_t left_site) const{
    if(truncated_variance.empty() or truncated_variance.size() < get_length() + 1) {
        tools::log->warn("Truncated variance container hasn't been initialized! You called this function prematurely");
        return std::numeric_limits<double>::infinity();
    }
    return truncated_variance[left_site + 1];
}
double class_state_finite::get_truncated_variance() const { return get_truncated_variance(get_position()); }
const std::vector<double> &class_state_finite::get_truncated_variances() const{return truncated_variance;}


size_t class_state_finite::num_sites_truncated(double threshold) const {
    auto truncation_errors = get_truncation_errors();
    size_t trunc_bond_count = (size_t)  std::count_if(truncation_errors.begin(), truncation_errors.end(),
                                                      [threshold](auto const& val){ return val > threshold ; });
    return trunc_bond_count;
}


size_t class_state_finite::num_bonds_at_limit() const {
    auto bond_dims = tools::finite::measure::bond_dimensions(*this);
    size_t bonds_at_lim = (size_t)  std::count_if(bond_dims.begin(), bond_dims.end(),
                                                       [this](auto const& val){ return val >= (size_t)get_chi_lim(); });
    return bonds_at_lim;
}

bool class_state_finite::is_bond_limited(double threshold) const {
    return num_sites_truncated(threshold) > 0 and num_bonds_at_limit() > 0;
}


void class_state_finite::clear_measurements() const { measurements = Measurements(); }

void class_state_finite::clear_cache() const { cache = Cache(); }

void class_state_finite::do_all_measurements() const {
    tools::log->trace("Doing all measurements on state...");
    measurements.length                        = tools::finite::measure::length(*this);
    measurements.bond_dimension_current        = tools::finite::measure::bond_dimension_current(*this);
    measurements.bond_dimension_midchain       = tools::finite::measure::bond_dimension_midchain(*this);
    measurements.bond_dimensions               = tools::finite::measure::bond_dimensions(*this);
    measurements.norm                          = tools::finite::measure::norm(*this);
    measurements.energy                        = tools::finite::measure::energy(*this);
    measurements.energy_per_site               = tools::finite::measure::energy_per_site(*this);
    measurements.energy_variance               = tools::finite::measure::energy_variance(*this);
    measurements.energy_variance_per_site      = tools::finite::measure::energy_variance_per_site(*this);
    measurements.spin_components               = tools::finite::measure::spin_components(*this); // This will automatically measure sx,sy and sz as well
    measurements.entanglement_entropy_current  = tools::finite::measure::entanglement_entropy_current(*this);
    measurements.entanglement_entropy_midchain = tools::finite::measure::entanglement_entropy_midchain(*this);
    measurements.entanglement_entropies        = tools::finite::measure::entanglement_entropies(*this);
    tools::log->trace("Doing all measurements on state... OK");
}

void class_state_finite::tag_active_sites_have_been_updated(bool tag) const {
    if(site_update_tags.size() != get_length()) throw std::runtime_error("Cannot tag active sites, size mismatch in site list");
    for(auto &site : active_sites) {
        site_update_tags[site] = tag;
    }
}

void class_state_finite::tag_all_sites_have_been_updated(bool tag) const {
    if(site_update_tags.size() != get_length()) throw std::runtime_error("Cannot untag all sites, size mismatch in site list");
    site_update_tags = std::vector<bool>(get_length(), tag);
}

bool class_state_finite::all_sites_updated() const {
    if(site_update_tags.size() != get_length()) throw std::runtime_error("Cannot check update status on all sites, size mismatch in site list");
    return std::all_of(site_update_tags.begin(), site_update_tags.end(), [](bool v) { return v; });
}

bool class_state_finite::any_sites_updated() const {
    if(site_update_tags.size() != get_length()) throw std::runtime_error("Cannot check update status on all sites, size mismatch in site list");
    return get_sweeps() > 0 and std::any_of(site_update_tags.begin(), site_update_tags.end(), [](bool v) { return v; });
}

bool class_state_finite::active_sites_updated() const {
    if(site_update_tags.size() != get_length()) throw std::runtime_error("Cannot check update status on all sites, size mismatch in site list");
    if(active_sites.empty()) return false;
    auto first_site_ptr = site_update_tags.begin() + active_sites.front();
    auto last_site_ptr  = first_site_ptr + active_sites.size() - 1;
    tools::log->trace("Checking update status on active sites: {}", active_sites);
    bool updated = std::all_of(first_site_ptr, last_site_ptr, [](bool v) { return v; });
    tools::log->trace("Active sites updated: {}", updated);
    return updated;
}
