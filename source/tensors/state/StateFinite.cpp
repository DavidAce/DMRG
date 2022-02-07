#include <math/tenx.h>

// -- (textra first)
#include "StateFinite.h"
#include <config/settings.h>
#include <debug/exceptions.h>
#include <general/iter.h>
#include <tensors/site/mps/MpsSite.h>
#include <tid/tid.h>
#include <tools/common/contraction.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>
#include <tools/finite/multisite.h>

StateFinite::StateFinite() = default; // Can't initialize lists since we don't know the model size yet

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
StateFinite::~StateFinite() noexcept                   = default;            // default dtor
StateFinite::StateFinite(StateFinite &&other) noexcept = default;            // default move ctor
StateFinite &StateFinite::operator=(StateFinite &&other) noexcept = default; // default move assign

/* clang-format off */
StateFinite::StateFinite(const StateFinite &other):
    direction(other.direction),
    cache(other.cache),
    tag_normalized_sites(other.tag_normalized_sites),
    name(other.name),
    algo(other.algo),
    active_sites(other.active_sites),
    measurements(other.measurements)
{
    mps_sites.clear();
    mps_sites.reserve(other.mps_sites.size());
    for(const auto &mps : other.mps_sites) mps_sites.emplace_back(std::make_unique<MpsSite>(*mps));
}
/* clang-format on */

StateFinite &StateFinite::operator=(const StateFinite &other) {
    // check for self-assignment
    if(this != &other) {
        direction            = other.direction;
        cache                = other.cache;
        tag_normalized_sites = other.tag_normalized_sites;
        name                 = other.name;
        algo                 = other.algo;
        active_sites         = other.active_sites;
        measurements         = other.measurements;
        mps_sites.clear();
        for(const auto &mps : other.mps_sites) mps_sites.emplace_back(std::make_unique<MpsSite>(*mps));
    }
    return *this;
}

StateFinite::StateFinite(AlgorithmType algo_type, ModelType model_type, size_t model_size, size_t position) {
    initialize(algo_type, model_type, model_size, position);
}

void StateFinite::initialize(AlgorithmType algo_type, ModelType model_type, size_t model_size, size_t position) {
    tools::log->debug("Initializing state: algorithm [{}] | model [{}] | sites [{}] | position [{}]", enum2sv(algo_type), enum2sv(model_type), model_size,
                      position);
    set_algorithm(algo_type);
    if(model_size < 2) throw std::logic_error("Tried to initialize state with less than 2 sites");
    if(model_size > 2048) throw std::logic_error("Tried to initialize state with more than 2048 sites");
    if(position >= model_size) throw std::logic_error("Tried to initialize state at a position larger than the number of sites");

    long spin_dim = 2;
    switch(model_type) {
        case ModelType::ising_tf_rf: spin_dim = settings::model::ising_tf_rf::spin_dim; break;
        case ModelType::ising_sdual: spin_dim = settings::model::ising_sdual::spin_dim; break;
        case ModelType::lbit: spin_dim = settings::model::lbit::spin_dim; break;
        default: spin_dim = 2;
    }

    mps_sites.clear();

    // Generate a simple state with all spins equal
    Eigen::Tensor<Scalar, 3> M(static_cast<long>(spin_dim), 1, 1);
    Eigen::Tensor<Scalar, 1> L(1);
    M(0, 0, 0)        = 0;
    M(1, 0, 0)        = 1;
    L(0)              = 1;
    std::string label = "A";
    for(size_t site = 0; site < model_size; site++) {
        mps_sites.emplace_back(std::make_unique<MpsSite>(M, L, site, 0.0, label));
        if(site == position) {
            mps_sites.back()->set_LC(L);
            label = "B";
        }
    }
    if(mps_sites.size() != model_size) throw std::logic_error("Initialized state with wrong size");
    if(not get_mps_site(position).isCenter()) throw std::logic_error("Initialized state center bond at the wrong position");
    if(get_position() != position) throw std::logic_error("Initialized state at the wrong position");
    tag_normalized_sites = std::vector<bool>(model_size, false);
    //    tag_edge_ene_status  = std::vector<EdgeStatus> (model_size,EdgeStatus::STALE);
    //    tag_edge_var_status  = std::vector<EdgeStatus> (model_size,EdgeStatus::STALE);
}

void             StateFinite::set_name(std::string_view statename) { name = statename; }
std::string_view StateFinite::get_name() const { return name; }

void          StateFinite::set_algorithm(const AlgorithmType &algo_type) { algo = algo_type; }
AlgorithmType StateFinite::get_algorithm() const { return algo; }

void StateFinite::set_positions() {
    for(auto &&[pos, mps] : iter::enumerate(mps_sites)) mps->set_position(pos);
}

template<typename T>
T StateFinite::get_length() const {
    return static_cast<T>(mps_sites.size());
}
template double StateFinite::get_length<double>() const;
template size_t StateFinite::get_length<size_t>() const;
template long   StateFinite::get_length<long>() const;
template int    StateFinite::get_length<int>() const;

template<typename T>
T StateFinite::get_position() const {
    std::optional<T> pos;
    for(const auto &mps : mps_sites)
        if(mps->isCenter()) {
            if(pos) throw except::logic_error("Found multiple centers: first center at {} and another at {}", pos.value(), mps->get_position());
            pos = mps->get_position<T>();
        }
    // If no center position was found then all sites are "B" sites. In that case, return -1 if T is signed, otherwise throw.
    if(not pos) {
        if constexpr(std::is_signed_v<T>)
            return -1;
        else
            throw except::runtime_error("could not find center position in current state: {}\n"
                                        "hint: Call get_position<long>() to get -1",
                                        get_labels());
    } else
        return pos.value();
}
template size_t StateFinite::get_position<size_t>() const;
template long   StateFinite::get_position<long>() const;

long StateFinite::find_largest_bond() const {
    auto bond_dimensions = tools::finite::measure::bond_dimensions(*this);
    return *max_element(std::begin(bond_dimensions), std::end(bond_dimensions));
}

int                           StateFinite::get_direction() const { return direction; }
std::vector<std::string_view> StateFinite::get_labels() const {
    std::vector<std::string_view> labels;
    labels.reserve(get_length());
    for(const auto &mps : mps_sites) labels.emplace_back(mps->get_label());
    return labels;
}

void StateFinite::flip_direction() { direction *= -1; }

std::array<long, 3> StateFinite::dimensions_1site() const {
    auto pos = get_position<long>();
    if(pos >= 0)
        return get_mps_site(pos).dimensions();
    else
        return {0, 0, 0};
}

std::array<long, 3> StateFinite::dimensions_2site() const {
    std::array<long, 3> dimensions{};
    auto                pos  = get_position<long>();
    auto                posL = std::clamp<long>(pos, 0, get_length<long>() - 2);
    auto                posR = std::clamp<long>(pos + 1, 0, get_length<long>() - 1);
    const auto         &mpsL = get_mps_site(posL);
    const auto         &mpsR = get_mps_site(posR);
    dimensions[1]            = mpsL.get_chiL();
    dimensions[2]            = mpsR.get_chiR();
    dimensions[0]            = posL != posR ? mpsL.spin_dim() * mpsR.spin_dim() : mpsL.spin_dim();
    return dimensions;
}

std::array<long, 3> StateFinite::dimensions_nsite() const { return tools::finite::multisite::get_dimensions(*this, active_sites); }

long StateFinite::size_1site() const {
    auto dims = dimensions_1site();
    return dims[0] * dims[1] * dims[2];
}

long StateFinite::size_2site() const {
    auto dims = dimensions_2site();
    return dims[0] * dims[1] * dims[2];
}

long StateFinite::size_nsite() const {
    auto dims = dimensions_nsite();
    return dims[0] * dims[1] * dims[2];
}

bool StateFinite::position_is_the_middle() const { return get_position<long>() + 1 == get_length<long>() / 2 and direction == 1; }
bool StateFinite::position_is_the_middle_any_direction() const { return get_position<long>() + 1 == get_length<long>() / 2; }

bool StateFinite::position_is_outward_edge_left([[maybe_unused]] size_t nsite) const {
    if(nsite == 1) {
        return get_position<long>() <= -1 and direction == -1; // i.e. all sites are B's
    } else
        return get_position<long>() == 0 and direction == -1 and get_mps_site().isCenter(); // left-most site is a an AC
}

bool StateFinite::position_is_outward_edge_right(size_t nsite) const {
    return get_position<long>() >= get_length<long>() - static_cast<long>(nsite) and direction == 1;
}

bool StateFinite::position_is_outward_edge(size_t nsite) const { return position_is_outward_edge_left(nsite) or position_is_outward_edge_right(nsite); }

bool StateFinite::position_is_inward_edge_left([[maybe_unused]] size_t nsite) const {
    return get_position<long>() == 0 and direction == 1; // i.e. first site is an AC going to the right
}

bool StateFinite::position_is_inward_edge_right(size_t nsite) const {
    return get_position<long>() >= get_length<long>() - static_cast<long>(nsite) and direction == -1;
}

bool StateFinite::position_is_inward_edge(size_t nsite) const { return position_is_inward_edge_left(nsite) or position_is_inward_edge_right(nsite); }

bool StateFinite::position_is_at(long pos) const { return get_position<long>() == pos; }

bool StateFinite::position_is_at(long pos, int dir) const { return get_position<long>() == pos and get_direction() == dir; }

bool StateFinite::position_is_at(long pos, int dir, bool isCenter) const {
    return get_position<long>() == pos and get_direction() == dir and (pos >= 0) == isCenter;
}

bool StateFinite::has_center_point() const { return get_position<long>() >= 0; }

bool StateFinite::is_real() const {
    bool mps_real = true;
    for(const auto &mps : mps_sites) mps_real = mps_real and mps->is_real();
    return mps_real;
}

bool StateFinite::has_nan() const {
    for(const auto &mps : mps_sites)
        if(mps->has_nan()) return true;
    return false;
}

void StateFinite::assert_validity() const {
    size_t pos = 0;
    for(const auto &mps : mps_sites) {
        if(pos != mps->get_position<size_t>())
            throw std::runtime_error(fmt::format("State is corrupted: position mismatch: expected position {} != mps position {}", pos, mps->get_position()));
        pos++;
    }

    for(const auto &mps : mps_sites) mps->assert_validity();
    if(settings::model::model_type == ModelType::ising_sdual or settings::model::model_type == ModelType::ising_majorana) {
        for(const auto &mps : mps_sites)
            if(not mps->is_real()) throw except::runtime_error("state has imaginary part at mps position {}", mps->get_position());
    }
}

const Eigen::Tensor<StateFinite::Scalar, 1> &StateFinite::midchain_bond() const {
    auto pos = get_position<long>();
    auto cnt = (get_length<long>() - 1) / 2;
    if(pos < cnt) return get_mps_site(cnt).get_L();
    if(pos > cnt) return get_mps_site(cnt + 1).get_L();
    return get_mps_site(cnt).get_LC();
}

const Eigen::Tensor<StateFinite::Scalar, 1> &StateFinite::current_bond() const { return get_mps_site(get_position()).get_LC(); }

template<typename T>
const MpsSite &StateFinite::get_mps_site(T pos) const {
    if constexpr(std::is_signed_v<T>)
        if(pos < 0) throw std::range_error(fmt::format("get_mps_site(pos): pos out of range: {}", pos));
    if(pos >= get_length<T>()) throw std::range_error(fmt::format("get_mps_site(pos): pos out of range: {}", pos));
    const auto &mps_ptr = *std::next(mps_sites.begin(), static_cast<long>(pos));
    if(mps_ptr->get_position<T>() != pos)
        throw std::range_error(fmt::format("get_mps_site(pos): mismatch pos {} != mps pos {}", pos, mps_ptr->get_position<T>()));
    return *mps_ptr;

    //    if(algo == AlgorithmType::fLBIT){
    //        // During fLBIT we can't assume the mps positions are in order since we could have swap operators.
    //        for(const auto & mps_ptr :  mps_sites ){
    //            if(mps_ptr->get_position<T>() == pos) return *mps_ptr;
    //        }
    //        throw std::runtime_error(fmt::format("get_mps_site(pos): pos {} not found", pos));
    //    }else{
    //        // There shouldn't be any swap operator, we can safely assume the mps positions are sorted
    //        const auto &mps_ptr = *std::next(mps_sites.begin(), static_cast<long>(pos));
    //        if(mps_ptr->get_position<T>() != pos)
    //            throw std::range_error(fmt::format("get_mps_site(pos): mismatch pos {} != mps pos {}", pos, mps_ptr->get_position<T>()));
    //        return *mps_ptr;
    //    }
}
template const MpsSite &StateFinite::get_mps_site(size_t pos) const;
template const MpsSite &StateFinite::get_mps_site(long pos) const;

template<typename T>
MpsSite &StateFinite::get_mps_site(T pos) {
    return const_cast<MpsSite &>(std::as_const(*this).get_mps_site<T>(pos));
}
template MpsSite &StateFinite::get_mps_site(size_t pos);
template MpsSite &StateFinite::get_mps_site(long pos);

const MpsSite &StateFinite::get_mps_site() const { return get_mps_site(get_position()); }

MpsSite &StateFinite::get_mps_site() { return get_mps_site(get_position()); }

std::vector<MpsSite> StateFinite::get_mps_sites(const std::vector<size_t> &sites) const {
    std::vector<MpsSite> mps_at_sites;
    for(const auto &site : sites) mps_at_sites.emplace_back(get_mps_site(site));
    return mps_at_sites;
}
void StateFinite::set_mps_sites(const std::vector<MpsSite> &mps_list) {
    for(const auto &mps_new : mps_list) {
        auto  pos     = mps_new.get_position();
        auto &mps_old = get_mps_site(pos);
        if(mps_new.has_L() and mps_new.has_M() and mps_old.get_label() == mps_new.get_label())
            mps_old = mps_new;
        else {
            mps_old.set_label(mps_new.get_label());
            if(mps_new.has_M()) mps_old.set_M(mps_new.get_M_bare());
            if(mps_new.has_L()) mps_old.set_L(mps_new.get_L());
            if(mps_new.has_LC()) mps_old.set_LC(mps_new.get_LC());
        }
    }
}

std::array<long, 3> StateFinite::active_dimensions() const { return tools::finite::multisite::get_dimensions(*this, active_sites); }

long StateFinite::active_problem_size() const { return tools::finite::multisite::get_problem_size(*this, active_sites); }

std::vector<long> StateFinite::get_spin_dims(const std::vector<size_t> &sites) const {
    if(sites.empty()) throw std::runtime_error("No sites on which to collect spin dimensions");
    std::vector<long> dims;
    dims.reserve(sites.size());
    for(const auto &site : sites) { dims.emplace_back(get_mps_site(site).spin_dim()); }
    return dims;
}

std::vector<long> StateFinite::get_spin_dims() const { return get_spin_dims(active_sites); }

Eigen::Tensor<StateFinite::Scalar, 3> StateFinite::get_multisite_mps(const std::vector<size_t> &sites) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mps tensor");
    if(sites == active_sites and cache.multisite_mps) return cache.multisite_mps.value();
    if constexpr(settings::debug) tools::log->trace("get_multisite_mps: sites {}", sites);
    auto                     t_mps = tid::tic_scope("gen_mps");
    Eigen::Tensor<Scalar, 3> multisite_mps;
    Eigen::Tensor<Scalar, 3> temp;
    for(auto &site : sites) {
        if(&site == &sites.front()) { // First site
            multisite_mps = get_mps_site(site).get_M();
            continue;
        }
        const auto &M = get_mps_site(site).get_M();
        multisite_mps = tools::common::contraction::contract_mps_mps_temp(multisite_mps, M, temp);
    }
    if(sites.front() != 0 and get_mps_site(sites.front()).get_label() == "B") {
        // In this case all sites are "B" and we need to prepend the the "L" from the site on the left to make a normalized multisite mps
        auto &mps_left = get_mps_site(sites.front() - 1);
        auto &L_left   = mps_left.isCenter() ? mps_left.get_LC() : mps_left.get_L();
        if(L_left.dimension(0) != multisite_mps.dimension(1))
            throw std::logic_error(
                fmt::format("get_multisite_mps: mismatching dimensions: L_left {} | multisite_mps {}", L_left.dimensions(), multisite_mps.dimensions()));
        multisite_mps = tools::common::contraction::contract_bnd_mps_temp(L_left, multisite_mps, temp);
    } else if(sites.back() < get_length() - 1 and get_mps_site(sites.back()).get_label() == "A") {
        // In this case all sites are "A" and we need to append the the "L" from the site on the right to make a normalized multisite mps
        auto &mps_right = get_mps_site(sites.back() + 1);
        auto &L_right   = mps_right.get_L();
        if(L_right.dimension(0) != multisite_mps.dimension(2))
            throw std::logic_error(
                fmt::format("get_multisite_mps: mismatching dimensions: L_right {} | multisite_mps {}", L_right.dimensions(), multisite_mps.dimensions()));
        multisite_mps = tools::common::contraction::contract_mps_bnd_temp(multisite_mps, L_right, temp);
    }

    if constexpr(settings::debug) {
        // Check the norm of the tensor on debug builds
        auto   t_dbg = tid::tic_scope("debug");
        double norm  = tools::common::contraction::contract_mps_norm(multisite_mps);
        if(std::abs(norm - 1) > settings::precision::max_norm_error) {
            if(sites.front() != 0 and get_mps_site(sites.front()).get_label() == "B") {
                // In this case all sites are "B" and we need to prepend the the "L" from the site on the left
                auto &mps_left = get_mps_site(sites.front() - 1);
                auto &L_left   = mps_left.isCenter() ? mps_left.get_LC() : mps_left.get_L();
                multisite_mps  = tools::common::contraction::contract_bnd_mps_temp(L_left, multisite_mps, temp);
                norm           = tools::common::contraction::contract_mps_norm(multisite_mps);
                tools::log->critical("Norm after adding L to B from the left: {:.16f}", norm);
            }
            throw std::runtime_error(fmt::format("get_multisite_mps: not normalized: sites {} | norm ⟨ψ|ψ⟩ = {:.16f}", sites, norm));
        }
    }
    return multisite_mps;
}

const Eigen::Tensor<StateFinite::Scalar, 3> &StateFinite::get_multisite_mps() const {
    if(cache.multisite_mps) return cache.multisite_mps.value();
    cache.multisite_mps = get_multisite_mps(active_sites);
    return cache.multisite_mps.value();
}

void StateFinite::set_truncation_error(size_t pos, double error) { get_mps_site(pos).set_truncation_error(error); }
void StateFinite::set_truncation_error(double error) { set_truncation_error(get_position(), error); }

void StateFinite::set_truncation_error_LC(double error) {
    auto &mps = get_mps_site(get_position());
    if(not mps.isCenter()) throw std::runtime_error("mps at current position is not a center");
    mps.set_truncation_error_LC(error);
}

void StateFinite::keep_max_truncation_errors(std::vector<double> &other_errors) {
    auto errors = get_truncation_errors();
    if(other_errors.size() != errors.size()) throw except::runtime_error("keep_max_truncation_errors: size mismatch");
    std::vector<double> max_errors(errors.size(), 0);
    for(auto &&[i, e] : iter::enumerate(max_errors)) e = std::max({errors[i], other_errors[i]});
    other_errors = max_errors;

    // Now set the maximum errors back to each site
    size_t past_center = 0;
    for(const auto &mps : mps_sites) {
        auto pos = mps->get_position<size_t>();
        auto idx = pos + past_center;
        set_truncation_error(pos, max_errors[idx]);
        if(mps->isCenter()) {
            past_center = 1;
            idx         = pos + past_center;
            set_truncation_error_LC(max_errors[idx]);
        }
    }
}

double StateFinite::get_truncation_error(size_t pos) const { return get_mps_site(pos).get_truncation_error(); }

double StateFinite::get_truncation_error() const {
    auto pos = get_position<long>();
    if(pos >= 0)
        return get_mps_site(pos).get_truncation_error();
    else
        return 0;
}

double StateFinite::get_truncation_error_LC() const { return get_mps_site(get_position()).get_truncation_error_LC(); }
double StateFinite::get_truncation_error_midchain() const {
    auto pos = get_position<long>();
    auto cnt = (get_length<long>() - 1) / 2;
    if(pos < cnt) return get_mps_site(cnt).get_truncation_error();
    if(pos > cnt) return get_mps_site(cnt + 1).get_truncation_error();
    return get_mps_site(cnt).get_truncation_error_LC();
}

std::vector<double> StateFinite::get_truncation_errors() const { return tools::finite::measure::truncation_errors(*this); }
std::vector<double> StateFinite::get_truncation_errors_active() const { return tools::finite::measure::truncation_errors_active(*this); }
double              StateFinite::get_truncation_error_active_max() const {
    auto   truncation_errors_active = get_truncation_errors_active();
    double truncation_error         = 0;
    if(not truncation_errors_active.empty()) truncation_error = *std::max_element(truncation_errors_active.begin(), truncation_errors_active.end());
    return truncation_error;
}

size_t StateFinite::num_sites_truncated(double truncation_threshold) const {
    auto truncation_errors = get_truncation_errors();
    auto trunc_bond_count  = static_cast<size_t>(
        std::count_if(truncation_errors.begin(), truncation_errors.end(), [truncation_threshold](auto const &val) { return val > truncation_threshold; }));
    return trunc_bond_count;
}

size_t StateFinite::num_bonds_at_limit(long bond_level) const {
    auto bond_dims    = tools::finite::measure::bond_dimensions(*this);
    auto bonds_at_lim = static_cast<size_t>(std::count_if(bond_dims.begin(), bond_dims.end(), [bond_level](auto const &val) { return val >= bond_level; }));
    return bonds_at_lim;
}

bool StateFinite::is_limited_by_bond(long bond_limit, [[maybe_unused]] double truncation_threshold) const {
    return num_bonds_at_limit(bond_limit) > 0;
    //    return num_sites_truncated(truncation_threshold) > 0 or num_bonds_at_limit(bond_limit) > 0;
}

void StateFinite::clear_measurements(LogPolicy logPolicy) const {
    if(logPolicy == LogPolicy::NORMAL) tools::log->trace("Clearing state measurements");
    measurements = MeasurementsStateFinite();
}

void StateFinite::clear_cache(LogPolicy logPolicy) const {
    if(logPolicy == LogPolicy::NORMAL) tools::log->trace("Clearing state cache");
    cache = Cache();
}

void StateFinite::do_all_measurements() const { tools::finite::measure::do_all_measurements(*this); }

void StateFinite::tag_active_sites_normalized(bool tag) const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot tag active sites, size mismatch in site list");
    for(auto &site : active_sites) tag_normalized_sites[site] = tag;
}

void StateFinite::tag_all_sites_normalized(bool tag) const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot untag all sites, size mismatch in site list");
    tag_normalized_sites = std::vector<bool>(get_length(), tag);
}

void StateFinite::tag_site_normalized(size_t pos, bool tag) const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot untag all sites, size mismatch in site list");
    tag_normalized_sites[pos] = tag;
}

bool StateFinite::is_normalized_on_all_sites() const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot check normalization status on all sites, size mismatch in site list");
    return std::all_of(tag_normalized_sites.begin(), tag_normalized_sites.end(), [](bool v) { return v; });
}

bool StateFinite::is_normalized_on_any_sites() const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot check normalization status on any sites, size mismatch in site list");
    return std::any_of(tag_normalized_sites.begin(), tag_normalized_sites.end(), [](bool v) { return v; });
}

bool StateFinite::is_normalized_on_active_sites() const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot check normalization status on active sites, size mismatch in site list");
    if(active_sites.empty()) return false;
    auto first_site_ptr = std::next(tag_normalized_sites.begin(), static_cast<long>(active_sites.front()));
    auto last_site_ptr  = std::next(tag_normalized_sites.begin(), static_cast<long>(active_sites.back()));
    return std::all_of(first_site_ptr, last_site_ptr, [](bool v) { return v; });
}

bool StateFinite::is_normalized_on_non_active_sites() const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot check update status on all sites, size mismatch in site list");
    if(active_sites.empty()) return is_normalized_on_all_sites();
    for(size_t idx = 0; idx < get_length(); idx++)
        if(std::find(active_sites.begin(), active_sites.end(), idx) == active_sites.end() and not tag_normalized_sites[idx]) return false;
    return true;
}

std::vector<size_t> StateFinite::get_active_ids() const {
    std::vector<size_t> ids;
    ids.reserve(active_sites.size());
    for(const auto &pos : active_sites) ids.emplace_back(get_mps_site(pos).get_unique_id());
    return ids;
}
