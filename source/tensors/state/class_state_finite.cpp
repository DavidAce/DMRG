//
// Created by david on 2019-01-29.
//

#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>

// -- (textra first)
#include "class_state_finite.h"
#include <config/nmspc_settings.h>
#include <tensors/state/class_mps_site.h>
#include <tools/common/log.h>
#include <tools/common/fmt.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/multisite.h>

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
class_state_finite::~class_state_finite()                                       = default; // default dtor
class_state_finite::class_state_finite(class_state_finite &&other)              = default; // default move ctor
class_state_finite &class_state_finite::operator=(class_state_finite &&other)   = default; // default move assign

/* clang-format off */
class_state_finite::class_state_finite(const class_state_finite &other):
    direction(other.direction),
    cache(other.cache),
    tag_normalized_sites(other.tag_normalized_sites),
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
        direction                = other.direction;
        cache                    = other.cache;
        tag_normalized_sites     = other.tag_normalized_sites;
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
    tag_normalized_sites = std::vector<bool>(model_size, false);
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

//size_t class_state_finite::get_iteration() const { return iter; }
//size_t class_state_finite::reset_iter() { return iter = 0; }
//void   class_state_finite::set_iter(size_t iter_) { iter = iter_; }
//void   class_state_finite::increment_iter() { iter++; }
//
//size_t class_state_finite::get_step() const { return step; }
//size_t class_state_finite::reset_step() { return step = 0; }
//void   class_state_finite::set_step(size_t step_) { step = step_; }
//void   class_state_finite::increment_step() { step++; }

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

Eigen::DSizes<long, 3> class_state_finite::active_dimensions() const { return tools::finite::multisite::get_dimensions(*this, active_sites); }

long class_state_finite::active_problem_size() const { return tools::finite::multisite::get_problem_size(*this, active_sites); }

Eigen::Tensor<class_state_finite::Scalar, 3> class_state_finite::get_multisite_mps(const std::vector<size_t> &sites) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mps tensor");
    if(sites == active_sites and cache.multisite_mps) return cache.multisite_mps.value();
    tools::log->trace("Contracting multisite mps tensor with {} sites", sites.size());
    tools::common::profile::get_default_prof()["t_mps"]->tic();
    Eigen::Tensor<Scalar, 3> multisite_tensor;
    constexpr auto           shuffle_idx  = Textra::array4{0, 2, 1, 3};
    constexpr auto           contract_idx = Textra::idx({2}, {1});
    Textra::array3           new_dims;
    Eigen::Tensor<Scalar, 3> temp;
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
        temp.device(Textra::omp::getDevice()) = multisite_tensor.contract(M, contract_idx).shuffle(shuffle_idx).reshape(new_dims);
        multisite_tensor     = temp;
    }
    tools::common::profile::get_default_prof()["t_mps"]->toc();
    if constexpr(settings::debug) {
        // Check the norm of the tensor on debug builds
        tools::common::profile::get_default_prof()["t_chk"]->tic();
        Eigen::Tensor<Scalar,0> norm_scalar = multisite_tensor.contract(multisite_tensor.conjugate(), Textra::idx({0,1,2},{0,1,2}));
        double norm = std::abs(norm_scalar(0));
        if(std::abs(norm - 1) > settings::precision::max_norm_error)
            throw std::runtime_error(fmt::format("Multisite tensor is not normalized. Norm = {:.16f} {:+.16f}i", std::real(norm_scalar(0)), std::imag(norm_scalar(0))));
        tools::common::profile::get_default_prof()["t_chk"]->toc();
    }
    return multisite_tensor;
}

const Eigen::Tensor<class_state_finite::Scalar, 3> &class_state_finite::get_multisite_mps() const {
    if(cache.multisite_mps) return cache.multisite_mps.value();
    cache.multisite_mps = get_multisite_mps(active_sites);
    return cache.multisite_mps.value();
}


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
std::vector<double> class_state_finite::get_truncation_errors_active() const {
    std::vector<double> truncation_errors;
    truncation_errors.reserve(active_sites.size());
    for(auto && pos : active_sites){
        // We are only interested in the truncation on bonds that are updated
        // when operating on active_sites. This excludes the outer bonds.
        if(get_mps_site(pos).isCenter())
            truncation_errors.emplace_back(get_truncation_error_LC());
        if(pos == active_sites.front()) continue;
        if(pos == active_sites.back()) continue;
        truncation_errors.emplace_back(get_truncation_error(pos));
    }
    return truncation_errors;
}

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

void class_state_finite::clear_measurements(LogPolicy logPolicy) const {
    if (logPolicy == LogPolicy::NORMAL) tools::log->trace("Clearing state measurements");
    measurements = state_measure_finite();
}

void class_state_finite::clear_cache(LogPolicy logPolicy) const {
    if (logPolicy == LogPolicy::NORMAL) tools::log->trace("Clearing state cache");
    cache = Cache();
}

void class_state_finite::do_all_measurements() const { tools::finite::measure::do_all_measurements(*this); }

void class_state_finite::tag_active_sites_normalized(bool tag) const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot tag active sites, size mismatch in site list");
    for(auto &site : active_sites)  tag_normalized_sites[site] = tag;
}

void class_state_finite::tag_all_sites_normalized(bool tag) const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot untag all sites, size mismatch in site list");
    tag_normalized_sites = std::vector<bool>(get_length(), tag);
}

void class_state_finite::tag_site_normalized(size_t pos, bool tag) const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot untag all sites, size mismatch in site list");
    tag_normalized_sites[pos] = tag;
}

bool class_state_finite::is_normalized_on_all_sites() const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot check normalization status on all sites, size mismatch in site list");
    tools::log->trace("Checking update status on all sites", active_sites);
    return std::all_of(tag_normalized_sites.begin(), tag_normalized_sites.end(), [](bool v) { return v; });
}

bool class_state_finite::is_normalized_on_any_sites() const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot check normalization status on any sites, size mismatch in site list");
    return std::any_of(tag_normalized_sites.begin(), tag_normalized_sites.end(), [](bool v) { return v; });
}

bool class_state_finite::is_normalized_on_active_sites() const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot check normalization status on active sites, size mismatch in site list");
    if(active_sites.empty()) return false;
    auto first_site_ptr = std::next(tag_normalized_sites.begin(), static_cast<long>(active_sites.front()));
    auto last_site_ptr  = std::next(tag_normalized_sites.begin(), static_cast<long>(active_sites.back()));
    tools::log->trace("Checking update status on active sites: {}", active_sites);
    bool normalized = std::all_of(first_site_ptr, last_site_ptr, [](bool v) { return v; });
    tools::log->trace("Active sites normalized: {}", normalized);
    return normalized;
}

bool class_state_finite::is_normalized_on_non_active_sites() const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot check update status on all sites, size mismatch in site list");
    if(active_sites.empty()) return is_normalized_on_all_sites();
    tools::log->trace("Checking normalization status on non-active sites", active_sites);
    for(size_t idx = 0; idx < get_length(); idx++)
        if(std::find(active_sites.begin(),active_sites.end(),idx) == active_sites.end() and not tag_normalized_sites[idx]) return false;
    return true;
}