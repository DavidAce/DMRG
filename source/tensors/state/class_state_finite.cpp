//
// Created by david on 2019-01-29.
//

#include <general/nmspc_tensor_extra.h>
// -- (textra first)
#include "class_state_finite.h"
#include <config/nmspc_settings.h>
#include <tensors/state/class_mps_site.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/mpo.h>
#include <tools/finite/mps.h>
#include <tools/finite/multisite.h>
#include <tools/finite/svd.h>

class_state_finite::class_state_finite() = default; // Can't initialize lists since we don't know the model size yet

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_state_finite::~class_state_finite()                                   = default;            // default dtor
class_state_finite::class_state_finite(class_state_finite &&other) noexcept = default;            // default move ctor
class_state_finite &class_state_finite::operator=(class_state_finite &&other) noexcept = default; // default move assign

/* clang-format off */
class_state_finite::class_state_finite(const class_state_finite &other):
    iter(other.iter),
    step(other.step),
    direction(other.direction),
//    chi_lim(other.chi_lim),
//    chi_lim_max(other.cfg_chi_lim_max),
    cache(other.cache),
    site_update_tags(other.site_update_tags),
    active_sites(other.active_sites),
    measurements(other.measurements)
{
    mps_sites.clear();
    for(const auto &mps : other.mps_sites) mps_sites.emplace_back(std::make_unique<class_mps_site>(*mps));
}
/* clang-format on */

class_state_finite &class_state_finite::operator=(const class_state_finite &other) {
    // check for self-assignment
    if(this != &other) {
        iter                     = other.iter;
        step                     = other.step;
        direction                = other.direction;
//        chi_lim                  = other.chi_lim;
//        chi_lim_max              = other.cfg_chi_lim_max;
        cache                    = other.cache;
        site_update_tags         = other.site_update_tags;
        active_sites             = other.active_sites;
        measurements             = other.measurements;
        mps_sites.clear();
        for(const auto &mps : other.mps_sites) mps_sites.emplace_back(std::make_unique<class_mps_site>(*mps));
    }
    return *this;
}

void class_state_finite::initialize(ModelType model_type, size_t model_size, size_t position) {
    tools::log->info("Initializing state with {} sites at position {}", model_size, position);
    if(model_size < 2) throw std::logic_error("Tried to initialize state with less than 2 sites");
    if(model_size > 2048) throw std::logic_error("Tried to initialize state with more than 2048 sites");
    if(position >= model_size) throw std::logic_error("Tried to initialize state at a position larger than the number of sites");

    long spin_dim = 2;
    switch(model_type) {
        case ModelType::ising_tf_rf: spin_dim = settings::model::ising_tf_rf::spin_dim; break;
        case ModelType::ising_sdual: spin_dim = settings::model::ising_sdual::spin_dim; break;
        default: spin_dim = 2;
    }

    mps_sites.clear();

    // Generate a simple state with all spins equal
    Eigen::Tensor<Scalar, 3> M(static_cast<long>(spin_dim), 1, 1);
    Eigen::Tensor<Scalar, 1> L(1);
    M(0, 0, 0) = 0;
    M(1, 0, 0) = 1;
    L(0)       = 1;
    for(size_t site = 0; site < model_size; site++) {
        mps_sites.emplace_back(std::make_unique<class_mps_site>(M, L, site));
        if(site == position) mps_sites.back()->set_LC(L);
    }
    if(mps_sites.size() != model_size) throw std::logic_error("Initialized state with wrong size");
    if(not get_mps_site(position).isCenter()) throw std::logic_error("Initialized state center bond at the wrong position");
    if(get_position() != position) throw std::logic_error("Initialized state at the wrong position");
    site_update_tags = std::vector<bool>(model_size, false);
}

void class_state_finite::set_positions() {
    size_t pos = 0;
    for(auto &mps : mps_sites) mps->set_position(pos++);
}

size_t class_state_finite::get_length() const { return mps_sites.size(); }
size_t class_state_finite::get_position() const {
    size_t pos          = 0;
    bool   found_center = false;
    for(const auto &mps : mps_sites)
        if(mps->isCenter() and not found_center) {
            pos          = mps->get_position();
            found_center = true;
        } else if(mps->isCenter() and found_center)
            throw std::logic_error(fmt::format("Found multiple centers: first center at {} and another at {}", pos, mps->get_position()));
    if(not found_center) throw std::logic_error("Could not find center");
    return pos;
}

size_t class_state_finite::get_iteration() const { return iter; }
size_t class_state_finite::reset_iter() { return iter = 0; }
void   class_state_finite::set_iter(size_t iter_) { iter = iter_; }
void   class_state_finite::increment_iter() { iter++; }

size_t class_state_finite::get_step() const { return step; }
size_t class_state_finite::reset_step() { return step = 0; }
void   class_state_finite::set_step(size_t step_) { step = step_; }
void   class_state_finite::increment_step() { step++; }


//long class_state_finite::get_chi_lim() const {
//    // Should get the the current limit on allowed bond dimension
//    if(not chi_lim) throw std::runtime_error("Chi limit has not been set yet");
//    return chi_lim.value();
//}
//void class_state_finite::set_chi_lim(long chi_lim_) {
//    // Should set the the current limit on allowed bond dimension
//    if(chi_lim_ == 0) throw std::runtime_error("Can't set chi limit to zero!");
//    chi_lim = chi_lim_;
//}
//
//
//long class_state_finite::get_chi_lim_init() const {
//    // Should get the the current limit on allowed bond dimension for the duration of the simulation
//    if(not cfg_chi_lim_max) throw std::runtime_error("Chi maximum has not been set yet");
//    return cfg_chi_lim_max.value();
//}
//
//void class_state_finite::set_chi_lim_init(long chi_max_) {
//    // Should set the the highest limit on allowed bond dimension for the duration of the simulation
//    if(chi_max_ == 0) throw std::runtime_error("Can't set chi max to zero!");
//    cfg_chi_lim_max = chi_max_;
//}
//
//long class_state_finite::get_chi_lim_max() const {
//    // Should get the the current limit on allowed bond dimension for the duration of the simulation
//    if(not cfg_chi_lim_max) throw std::runtime_error("Chi maximum has not been set yet");
//    return cfg_chi_lim_max.value();
//}
//
//void class_state_finite::set_chi_lim_max(long chi_max_) {
//    // Should set the the highest limit on allowed bond dimension for the duration of the simulation
//    if(chi_max_ == 0) throw std::runtime_error("Can't set chi max to zero!");
//    cfg_chi_lim_max = chi_max_;
//}

long class_state_finite::find_largest_chi() const {
    auto bond_dimensions = tools::finite::measure::bond_dimensions(*this);
    return *max_element(std::begin(bond_dimensions), std::end(bond_dimensions));
}

int  class_state_finite::get_direction() const { return direction; }
void class_state_finite::flip_direction() { direction *= -1; }

Eigen::DSizes<long, 3> class_state_finite::dimensions_2site() const {
    Eigen::DSizes<long, 3> dimensions;
    auto                   pos  = get_position();
    const auto &           mpsL = get_mps_site(pos);
    const auto &           mpsR = get_mps_site(pos + 1);
    dimensions[1]               = mpsL.get_chiL();
    dimensions[2]               = mpsR.get_chiR();
    dimensions[0]               = mpsL.spin_dim() * mpsR.spin_dim();
    return dimensions;
}

long class_state_finite::size_2site() const {
    auto dims = dimensions_2site();
    return dims[0] * dims[1] * dims[2];
}

bool class_state_finite::position_is_the_middle() const {
    return (size_t) get_position() + 1 == (size_t)(static_cast<double>(get_length()) / 2.0) and direction == 1;
}
bool class_state_finite::position_is_the_middle_any_direction() const {
    return (size_t) get_position() + 1 == (size_t)(static_cast<double>(get_length()) / 2.0);
}

bool class_state_finite::position_is_left_edge() const { return get_position() == 0 and direction == -1; }

bool class_state_finite::position_is_right_edge() const { return get_position() == get_length() - 2 and direction == 1; }

bool class_state_finite::position_is_any_edge() const { return position_is_left_edge() or position_is_right_edge(); }

bool class_state_finite::position_is_at(size_t pos) const { return get_position() == pos; }

bool class_state_finite::is_real() const {
    bool mps_real = true;
    for(const auto &mps : mps_sites) mps_real = mps_real and mps->is_real();
    return mps_real;
}

bool class_state_finite::has_nan() const {
    for(const auto &mps : mps_sites)
        if(mps->has_nan()) return true;
    return false;
}

void class_state_finite::assert_validity() const {
    size_t pos = 0;
    for(const auto &mps : mps_sites){
        if(pos != mps->get_position())
            throw std::runtime_error(fmt::format("State is corrupted: position mismatch: expected position {} != mps position {}", pos, mps->get_position()));
        pos++;
    }
    for(const auto &mps : mps_sites) mps->assert_validity();
}

const Eigen::Tensor<class_state_finite::Scalar, 1> &class_state_finite::midchain_bond() const {
    size_t center_pos = (get_length() - 1) / 2;
    if(get_position() < center_pos) return get_mps_site(center_pos).get_L();
    if(get_position() > center_pos) return get_mps_site(center_pos + 1).get_L();
    if(get_position() == center_pos)
        return get_mps_site(center_pos).get_LC();
    else
        throw std::logic_error("No valid position to find midchain_bond");
}

const Eigen::Tensor<class_state_finite::Scalar, 1> &class_state_finite::current_bond() const { return get_mps_site(get_position()).get_LC(); }

const class_mps_site &class_state_finite::get_mps_site(size_t pos) const {
    if(pos >= get_length()) throw std::range_error(fmt::format("get_mps_site(pos): pos out of range: {}", pos));
    const auto &mps_ptr = *std::next(mps_sites.begin(), static_cast<long>(pos));
    if(mps_ptr->get_position() != pos) throw std::range_error(fmt::format("get_mps_site(pos): mismatch pos {} != mps pos {}", pos, mps_ptr->get_position()));
    return *mps_ptr;
}

class_mps_site &class_state_finite::get_mps_site(size_t pos) { return const_cast<class_mps_site &>(std::as_const(*this).get_mps_site(pos)); }

const class_mps_site &class_state_finite::get_mps_site() const { return get_mps_site(get_position()); }

class_mps_site &class_state_finite::get_mps_site() { return get_mps_site(get_position()); }

//
// const class_mpo_site &class_state_finite::get_2site_tensor(size_t pos) const {
//    if(pos >= MPO_L.size() + MPO_R.size()) throw std::range_error(fmt::format("get_2site_tensor(pos) pos out of range: {}", pos));
//    if(pos <= MPO_L.back()->get_position()) {
//        auto mpo_it = std::next(MPO_L.begin(), pos)->get();
//        if(mpo_it->get_position() != pos)
//            throw std::range_error(fmt::format("get_2site_tensor(pos): Mismatch in mpo position and pos: {} != {}", mpo_it->get_position(), pos));
//        return *mpo_it;
//    } else {
//        if(pos < MPO_R.front()->get_position())
//            throw std::range_error(fmt::format("get_mps_site(pos): Mismatch in pos and MPOR front position: {} < {}", pos, MPO_R.front()->get_position()));
//        auto mpo_it = std::next(MPO_R.begin(), pos - MPO_R.front()->get_position())->get();
//        if(mpo_it->get_position() != pos)
//            throw std::range_error(fmt::format("get_2site_tensor(pos): Mismatch in mpo position and pos: {} != {}", mpo_it->get_position(), pos));
//        return *mpo_it;
//    }
//}
//
// class_mpo_site &class_state_finite::get_2site_tensor(size_t pos) {
//    return const_cast<class_mpo_site &>(static_cast<const class_state_finite &>(*this).get_2site_tensor(pos));
//}
//
// const class_environment &class_state_finite::get_ENVL(size_t pos) const {
//    if(pos > ENV_L.back().get_position()) throw std::range_error(fmt::format("get_ENVL(pos):  pos is not in left side: {}", pos));
//    if(pos >= ENV_L.size()) throw std::range_error(fmt::format("get_ENVL(pos) pos out of range: {}", pos));
//    auto env_it = std::next(ENV_L.begin(), static_cast<long>(pos));
//    if(env_it->get_position() != pos)
//        throw std::range_error(fmt::format("get_ENVL(pos): Mismatch in env position and pos: {} != {}", env_it->get_position(), pos));
//    return *env_it;
//}
//
// const class_environment &class_state_finite::get_ENVR(size_t pos) const {
//    if(pos < ENV_R.front().get_position()) throw std::range_error(fmt::format("get_ENVR(pos):  pos is not in right side: {}", pos));
//    if(pos >= ENV_L.size() + ENV_R.size()) throw std::range_error(fmt::format("get_ENVR(pos):  pos out of range: {}", pos));
//    auto env_it = std::next(ENV_R.begin(), static_cast<long>(pos - ENV_R.front().get_position()));
//    if(env_it->get_position() != pos)
//        throw std::range_error(fmt::format("get_ENVR(pos): Mismatch in env position and pos: {} != {}", env_it->get_position(), pos));
//    return *env_it;
//}
//
// const class_environment_var &class_state_finite::get_ENV2L(size_t pos) const {
//    if(pos > ENV2_L.back().get_position()) throw std::range_error(fmt::format("get_ENV2L(pos):  pos is not in left side: {}", pos));
//    if(pos >= ENV2_L.size()) throw std::range_error(fmt::format("get_ENV2L(pos) pos out of range: {}", pos));
//    auto env2_it = std::next(ENV2_L.begin(), static_cast<long>(pos));
//    if(env2_it->get_position() != pos)
//        throw std::range_error(fmt::format("get_ENV2L(pos): Mismatch in env position and pos: {} != {}", env2_it->get_position(), pos));
//    return *env2_it;
//}
//
// const class_environment_var &class_state_finite::get_ENV2R(size_t pos) const {
//    if(pos < ENV2_R.front().get_position()) throw std::range_error(fmt::format("get_ENV2R(pos):  pos is not in right side: {}", pos));
//    if(pos >= ENV2_L.size() + ENV2_R.size()) throw std::range_error(fmt::format("get_ENV2R(pos):  pos out of range: {}", pos));
//    auto env2_it = std::next(ENV2_R.begin(), static_cast<long>(pos - ENV2_R.front().get_position()));
//    if(env2_it->get_position() != pos)
//        throw std::range_error(fmt::format("get_ENV2R(pos): Mismatch in env2 position and pos: {} != {}", env2_it->get_position(), pos));
//    return *env2_it;
//}
//
//// For reduced energy MPO's
//
// bool class_state_finite::is_reduced() const {
//    bool reduced = MPO_L.front()->is_reduced();
//    for(auto &mpo : MPO_L)
//        if(reduced != mpo->is_reduced()) {
//            throw std::runtime_error(
//                fmt::format("First MPO has is_reduced: {}, but MPO at pos {} has is_reduced: {}", reduced, mpo->get_position(), mpo->is_reduced()));
//        }
//    for(auto &mpo : MPO_R)
//        if(reduced != mpo->is_reduced()) {
//            throw std::runtime_error(
//                fmt::format("First MPO has is_reduced: {}, but MPO at pos {} has is_reduced: {}", reduced, mpo->get_position(), mpo->is_reduced()));
//        }
//    return reduced;
//}
//
// double class_state_finite::get_energy_reduced() const { return get_energy_per_site_reduced() * static_cast<double>(get_length()); }
//
// double class_state_finite::get_energy_per_site_reduced() const {
//    // Check that all energies are the same
//    double e_reduced = MPO_L.front()->get_reduced_energy();
//    for(auto &mpo : MPO_L) {
//        if(mpo->get_reduced_energy() != e_reduced) {
//            throw std::runtime_error("Reduced energy mismatch!");
//        }
//    }
//    for(auto &mpo : MPO_R) {
//        if(mpo->get_reduced_energy() != e_reduced) {
//            throw std::runtime_error("Reduced energy mismatch!");
//        }
//    }
//    return e_reduced;
//}
//
// void class_state_finite::set_reduced_energy(double total_energy) { set_reduced_energy_per_site(total_energy / static_cast<double>(get_length())); }
//
// void class_state_finite::set_reduced_energy_per_site(double site_energy) {
//    if(get_energy_per_site_reduced() == site_energy) return;
//    clear_measurements();
//    cache.multimpo = {};
//    for(auto &mpo : MPO_L) mpo->set_reduced_energy(site_energy);
//    for(auto &mpo : MPO_R) mpo->set_reduced_energy(site_energy);
//    tools::finite::mps::rebuild_all_edges(*this);
//}
//
// void class_state_finite::perturb_hamiltonian(double coupling_ptb, double field_ptb, PerturbMode perturbMode) {
//    tools::finite::mpo::perturb_hamiltonian(*this, coupling_ptb, field_ptb, perturbMode);
//}
//
// void class_state_finite::damp_hamiltonian(double coupling_damp, double field_damp) { tools::finite::mpo::damp_hamiltonian(*this, coupling_damp, field_damp);
// }
//
// bool class_state_finite::is_perturbed() const {
//    for(size_t pos = 0; pos < get_length(); pos++) {
//        if(get_2site_tensor(pos).is_perturbed()) return true;
//    }
//    return false;
//}
//
// bool class_state_finite::is_damped() const {
//    for(size_t pos = 0; pos < get_length(); pos++) {
//        if(get_2site_tensor(pos).is_damped()) return true;
//    }
//    return false;
//}

// std::list<size_t> class_state_finite::activate_sites(long threshold, size_t max_sites, size_t min_sites) {
//    clear_cache();
//    return active_sites = tools::finite::multisite::generate_site_list(*this, threshold, max_sites, min_sites);
//}
//
// std::list<size_t> class_state_finite::activate_truncated_sites(long threshold, long chi_lim, size_t max_sites, size_t min_sites) {
//    clear_cache();
//    return active_sites = tools::finite::multisite::generate_truncated_site_list(*this, threshold, chi_lim, max_sites, min_sites);
//}
//
Eigen::DSizes<long, 3> class_state_finite::active_dimensions() const { return tools::finite::multisite::get_dimensions(*this, active_sites); }

long class_state_finite::active_problem_size() const { return tools::finite::multisite::get_problem_size(*this, active_sites); }

Eigen::Tensor<class_state_finite::Scalar, 3> class_state_finite::get_multisite_tensor(const std::list<size_t> &sites) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mps tensor");
    if(sites == active_sites and cache.multisite_tensor) return cache.multisite_tensor.value();
    tools::log->trace("Contracting multisite mps tensor with {} sites", sites.size());
    tools::common::profile::t_mps->tic();
    Eigen::Tensor<Scalar, 3> multisite_tensor;
    constexpr auto           shuffle_idx  = Textra::array4{0, 2, 1, 3};
    constexpr auto           contract_idx = Textra::idx({2}, {1});
    Textra::array3           new_dims;
    Eigen::Tensor<Scalar, 3> temp;
    OMP                      omp;
    bool                     first = true;
    for(auto &site : sites) {
        if(first) {
            multisite_tensor = get_mps_site(site).get_M();
            first            = false;
            continue;
        }
        const auto &M    = get_mps_site(site).get_M();
        long        dim0 = multisite_tensor.dimension(0) * M.dimension(0);
        long        dim1 = multisite_tensor.dimension(1);
        long        dim2 = M.dimension(2);
        new_dims         = {dim0, dim1, dim2};
        temp.resize(new_dims);
        temp.device(omp.dev) = multisite_tensor.contract(M, contract_idx).shuffle(shuffle_idx).reshape(new_dims);
        multisite_tensor     = temp;
    }

    if constexpr(settings::debug) {
        // Check the norm of the tensor on debug builds
        Eigen::Tensor<Scalar,0> norm_scalar = multisite_tensor.contract(multisite_tensor.conjugate(), Textra::idx({0,1,2},{0,1,2}));
        double norm = std::abs(norm_scalar(0));
        if(std::abs(norm - 1) > settings::precision::max_norm_error)
            throw std::runtime_error(fmt::format("Multisite tensor is not normalized. Norm = {:.16f} {:+.16f}i", std::real(norm_scalar(0)), std::imag(norm_scalar(0))));
    }

    tools::common::profile::t_mps->toc();
    return multisite_tensor;
}

const Eigen::Tensor<class_state_finite::Scalar, 3> &class_state_finite::get_multisite_tensor() const {
    if(cache.multisite_tensor) return cache.multisite_tensor.value();
    cache.multisite_tensor = get_multisite_tensor(active_sites);
    return cache.multisite_tensor.value();
}

// const Eigen::Tensor<class_state_finite::Scalar, 4> &class_state_finite::get_multimpo() const {
//    if(cache.multimpo) return cache.multimpo.value();
//    tools::common::profile::t_mpo->tic();
//    tools::log->trace("Contracting multi mpo");
//    if(active_sites.empty()) {
//        throw std::runtime_error("No active sites on which to build multimpo");
//    }
//    Eigen::Tensor<Scalar, 4> multimpo;
//    bool                     first = true;
//    for(auto &site : active_sites) {
//        if(first) {
//            multimpo = get_2site_tensor(site).MPO();
//            first    = false;
//            continue;
//        }
//        auto &                   mpo  = get_2site_tensor(site).MPO();
//        long                     dim0 = multimpo.dimension(0);
//        long                     dim1 = mpo.dimension(1);
//        long                     dim2 = multimpo.dimension(2) * mpo.dimension(2);
//        long                     dim3 = multimpo.dimension(3) * mpo.dimension(3);
//        Eigen::Tensor<Scalar, 4> temp(dim0, dim1, dim2, dim3);
//        temp     = multimpo.contract(mpo, Textra::idx({1}, {0})).shuffle(Textra::array6{0, 3, 1, 4, 2, 5}).reshape(Textra::array4{dim0, dim1, dim2, dim3});
//        multimpo = temp;
//    }
//    tools::common::profile::t_mpo->toc();
//    cache.multimpo = multimpo;
//    return cache.multimpo.value();
//}

// std::pair<std::reference_wrapper<const class_environment>, std::reference_wrapper<const class_environment>> class_state_finite::get_multienv() const {
//    return std::make_pair(get_ENVL(active_sites.front()), get_ENVR(active_sites.back()));
//}
//
// std::pair<std::reference_wrapper<const class_environment_var>, std::reference_wrapper<const class_environment_var>> class_state_finite::get_multienv2() const
// {
//    return std::make_pair(get_ENV2L(active_sites.front()), get_ENV2R(active_sites.back()));
//}

void class_state_finite::set_truncation_error(size_t pos, double error) { get_mps_site(pos).set_truncation_error(error); }
void class_state_finite::set_truncation_error(double error) { set_truncation_error(get_position(), error); }

void class_state_finite::set_truncation_error_LC(double error) {
    auto &mps = get_mps_site(get_position());
    if(not mps.isCenter()) throw std::runtime_error("mps at current position is not a center");
    mps.set_truncation_error_LC(error);
}

double class_state_finite::get_truncation_error(size_t pos) const { return get_mps_site(pos).get_truncation_error(); }

double class_state_finite::get_truncation_error() const { return get_mps_site(get_position()).get_truncation_error(); }
double class_state_finite::get_truncation_error_LC() const { return get_mps_site(get_position()).get_truncation_error_LC(); }
double class_state_finite::get_truncation_error_midchain() const {
    size_t center_pos = (get_length() - 1) / 2;
    if(center_pos > get_position())
        return get_mps_site(center_pos).get_truncation_error();
    else if(center_pos == get_position())
        return get_mps_site(center_pos).get_truncation_error_LC();
    else
        return get_mps_site(center_pos + 1).get_truncation_error();
}

std::vector<double> class_state_finite::get_truncation_errors() const { return tools::finite::measure::truncation_errors(*this); }

void class_state_finite::set_truncated_variance(size_t left_site, double error) {
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
void class_state_finite::set_truncated_variance(double error) { set_truncated_variance(get_position(), error); }

double class_state_finite::get_truncated_variance(size_t left_site) const {
    if(truncated_variance.empty() or truncated_variance.size() < get_length() + 1) {
        tools::log->warn("Truncated variance container hasn't been initialized! You called this function prematurely");
        return std::numeric_limits<double>::infinity();
    }
    return truncated_variance[left_site + 1];
}
double                     class_state_finite::get_truncated_variance() const { return get_truncated_variance(get_position()); }
const std::vector<double> &class_state_finite::get_truncated_variances() const { return truncated_variance; }

size_t class_state_finite::num_sites_truncated(double truncation_threshold) const {
    auto truncation_errors = get_truncation_errors();
    auto trunc_bond_count =
        static_cast<size_t>(std::count_if(truncation_errors.begin(), truncation_errors.end(), [truncation_threshold](auto const &val) { return val > truncation_threshold; }));
    return trunc_bond_count;
}

size_t class_state_finite::num_bonds_reached_chi(long chi_level) const {
    auto bond_dims    = tools::finite::measure::bond_dimensions(*this);
    auto bonds_at_lim = static_cast<size_t>(std::count_if(bond_dims.begin(), bond_dims.end(), [chi_level](auto const &val) { return val >= chi_level; }));
    return bonds_at_lim;
}

bool class_state_finite::is_bond_limited(long chi_lim, double truncation_threshold) const { return num_sites_truncated(truncation_threshold) > 0 and num_bonds_reached_chi(chi_lim) > 0; }

void class_state_finite::clear_measurements() const { measurements = state_measure_finite(); }

void class_state_finite::clear_cache() const { cache = Cache(); }

void class_state_finite::do_all_measurements() const { tools::finite::measure::do_all_measurements(*this); }

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
    return get_iteration() > 0 and std::any_of(site_update_tags.begin(), site_update_tags.end(), [](bool v) { return v; });
}

bool class_state_finite::active_sites_updated() const {
    if(site_update_tags.size() != get_length()) throw std::runtime_error("Cannot check update status on all sites, size mismatch in site list");
    if(active_sites.empty()) return false;
    auto first_site_ptr = std::next(site_update_tags.begin(), static_cast<long>(active_sites.front()));
    auto last_site_ptr  = std::next(site_update_tags.begin(), static_cast<long>(active_sites.back()));
    tools::log->trace("Checking update status on active sites: {}", active_sites);
    bool updated = std::all_of(first_site_ptr, last_site_ptr, [](bool v) { return v; });
    tools::log->trace("Active sites updated: {}", updated);
    return updated;
}
