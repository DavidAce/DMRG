//
// Created by david on 2019-01-29.
//

#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>

// -- (textra first)
#include "class_state_finite.h"
#include <config/nmspc_settings.h>
#include <tensors/state/class_mps_site.h>
#include <tools/common/fmt.h>
#include <tools/common/log.h>
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
class_state_finite::~class_state_finite()                          = default;            // default dtor
class_state_finite::class_state_finite(class_state_finite &&other) = default;            // default move ctor
class_state_finite &class_state_finite::operator=(class_state_finite &&other) = default; // default move assign

/* clang-format off */
class_state_finite::class_state_finite(const class_state_finite &other):
    direction(other.direction),
    cache(other.cache),
    tag_normalized_sites(other.tag_normalized_sites),
    active_sites(other.active_sites),
    measurements(other.measurements)
{
    mps_sites.clear();
    mps_sites.reserve(other.mps_sites.size());
    for(const auto &mps : other.mps_sites) mps_sites.emplace_back(std::make_unique<class_mps_site>(*mps));
}
/* clang-format on */

class_state_finite &class_state_finite::operator=(const class_state_finite &other) {
    // check for self-assignment
    if(this != &other) {
        direction            = other.direction;
        cache                = other.cache;
        tag_normalized_sites = other.tag_normalized_sites;
        active_sites = other.active_sites;
        measurements = other.measurements;
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
        mps_sites.emplace_back(std::make_unique<class_mps_site>(M, L, site, 0.0, label));
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

void        class_state_finite::set_name(const std::string &statename) { name = statename; }
std::string class_state_finite::get_name() const { return name; }

void class_state_finite::set_positions() {
    long pos = 0;
    for(auto &mps : mps_sites) mps->set_position(pos++);
}

template<typename T>
T class_state_finite::get_length() const {
    return static_cast<T>(mps_sites.size());
}
template size_t class_state_finite::get_length<size_t>() const;
template long   class_state_finite::get_length<long>() const;
template int    class_state_finite::get_length<int>() const;

template<typename T>
T class_state_finite::get_position() const {
    std::optional<T> pos;
    for(const auto &mps : mps_sites)
        if(mps->isCenter()) {
            if(pos) throw std::logic_error(fmt::format("Found multiple centers: first center at {} and another at {}", pos.value(), mps->get_position()));
            pos = mps->get_position<T>();
        }
    // If no center position was found then all sites are "B" sites. In that case, return -1 if T is signed, otherwise throw.
    if(not pos) {
        if constexpr(std::is_signed_v<T>)
            return -1;
        else
            throw std::runtime_error(fmt::format("Could not find center position in current state: {}", get_labels()));
    } else
        return pos.value();
}
template size_t class_state_finite::get_position<size_t>() const;
template long   class_state_finite::get_position<long>() const;

long class_state_finite::find_largest_chi() const {
    auto bond_dimensions = tools::finite::measure::bond_dimensions(*this);
    return *max_element(std::begin(bond_dimensions), std::end(bond_dimensions));
}

int                      class_state_finite::get_direction() const { return direction; }
std::vector<std::string> class_state_finite::get_labels() const {
    std::vector<std::string> labels;
    labels.reserve(get_length());
    for(const auto &mps : mps_sites) labels.emplace_back(mps->get_label());
    return labels;
}

void class_state_finite::flip_direction() { direction *= -1; }


std::array<long, 3> class_state_finite::dimensions_1site() const {
    auto pos = get_position<long>();
    if(pos >= 0)
        return get_mps_site(pos).dimensions();
    else
        return {0,0,0};
}

std::array<long, 3> class_state_finite::dimensions_2site() const {
    std::array<long, 3> dimensions{};
    auto                   pos  = get_position<long>();
    auto                   posL = std::clamp<long>(pos, 0, get_length<long>() - 2);
    auto                   posR = std::clamp<long>(pos + 1, 0, get_length<long>() - 1);
    const auto &           mpsL = get_mps_site(posL);
    const auto &           mpsR = get_mps_site(posR);
    dimensions[1]               = mpsL.get_chiL();
    dimensions[2]               = mpsR.get_chiR();
    dimensions[0]               = posL != posR ? mpsL.spin_dim() * mpsR.spin_dim() : mpsL.spin_dim();
    return dimensions;
}

std::array<long, 3> class_state_finite::dimensions_nsite() const {
    return tools::finite::multisite::get_dimensions(*this,active_sites);
}

long class_state_finite::size_1site() const {
    auto dims = dimensions_1site();
    return dims[0] * dims[1] * dims[2];
}

long class_state_finite::size_2site() const {
    auto dims = dimensions_2site();
    return dims[0] * dims[1] * dims[2];
}


long class_state_finite::size_nsite() const {
    auto dims = dimensions_nsite();
    return dims[0] * dims[1] * dims[2];
}


bool class_state_finite::position_is_the_middle() const { return get_position() + 1 == static_cast<size_t>(get_length() / 2) and direction == 1; }
bool class_state_finite::position_is_the_middle_any_direction() const { return get_position() + 1 == static_cast<size_t>(get_length() / 2); }

bool class_state_finite::position_is_outward_edge_left([[maybe_unused]] size_t nsite) const {
    if(nsite == 1) {
        return get_position<long>() <= -1 and direction == -1; // i.e. all sites are B's
    } else
        return get_position<long>() == 0 and direction == -1 and get_mps_site().isCenter(); // left-most site is a an AC
}

bool class_state_finite::position_is_outward_edge_right(size_t nsite) const {
    return get_position<long>() >= get_length<long>() - static_cast<long>(nsite) and direction == 1;
}

bool class_state_finite::position_is_outward_edge(size_t nsite) const { return position_is_outward_edge_left(nsite) or position_is_outward_edge_right(nsite); }

bool class_state_finite::position_is_inward_edge_left([[maybe_unused]] size_t nsite) const {
    return get_position<long>() == 0 and direction == 1; // i.e. first site is an AC going to the right
}

bool class_state_finite::position_is_inward_edge_right(size_t nsite) const {
    return get_position<long>() >= get_length<long>() - static_cast<long>(nsite) and direction == -1;
}

bool class_state_finite::position_is_inward_edge(size_t nsite) const { return position_is_inward_edge_left(nsite) or position_is_inward_edge_right(nsite); }

bool class_state_finite::position_is_at(long pos) const { return get_position<long>() == pos; }

bool class_state_finite::position_is_at(long pos, int dir) const { return get_position<long>() == pos and get_direction() == dir; }

bool class_state_finite::position_is_at(long pos, int dir, bool isCenter) const {
    return get_position<long>() == pos and get_direction() == dir and (pos >= 0) == isCenter;
}

bool class_state_finite::has_center_point() const { return get_position<long>() >= 0; }

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
    for(const auto &mps : mps_sites) {
        if(pos != mps->get_position<size_t>())
            throw std::runtime_error(fmt::format("State is corrupted: position mismatch: expected position {} != mps position {}", pos, mps->get_position()));
        pos++;
    }
    for(const auto &mps : mps_sites) mps->assert_validity();
}

const Eigen::Tensor<class_state_finite::Scalar, 1> &class_state_finite::midchain_bond() const {
    auto pos = get_position<long>();
    auto cnt = (get_length<long>() - 1) / 2;
    if(pos < cnt) return get_mps_site(cnt).get_L();
    if(pos > cnt) return get_mps_site(cnt + 1).get_L();
    return get_mps_site(cnt).get_LC();
}

const Eigen::Tensor<class_state_finite::Scalar, 1> &class_state_finite::current_bond() const { return get_mps_site(get_position()).get_LC(); }

template<typename T>
const class_mps_site &class_state_finite::get_mps_site(T pos) const {
    if constexpr(std::is_signed_v<T>)
        if(pos < 0) throw std::range_error(fmt::format("get_mps_site(pos): pos out of range: {}", pos));
    if(pos >= get_length<T>()) throw std::range_error(fmt::format("get_mps_site(pos): pos out of range: {}", pos));
    const auto &mps_ptr = *std::next(mps_sites.begin(), static_cast<long>(pos));
    if(mps_ptr->get_position<T>() != pos)
        throw std::range_error(fmt::format("get_mps_site(pos): mismatch pos {} != mps pos {}", pos, mps_ptr->get_position<T>()));
    return *mps_ptr;
}
template const class_mps_site &class_state_finite::get_mps_site(size_t pos) const;
template const class_mps_site &class_state_finite::get_mps_site(long pos) const;

template<typename T>
class_mps_site &class_state_finite::get_mps_site(T pos) {
    return const_cast<class_mps_site &>(std::as_const(*this).get_mps_site<T>(pos));
}
template class_mps_site &class_state_finite::get_mps_site(size_t pos);
template class_mps_site &class_state_finite::get_mps_site(long pos);

const class_mps_site &class_state_finite::get_mps_site() const { return get_mps_site(get_position()); }

class_mps_site &class_state_finite::get_mps_site() { return get_mps_site(get_position()); }

std::vector<class_mps_site> class_state_finite::get_mps_sites(const std::vector<size_t> & sites) const{
    std::vector<class_mps_site> mps_at_sites;
    for(const auto & site: sites) mps_at_sites.emplace_back(get_mps_site(site));
    return mps_at_sites;
}
void class_state_finite::set_mps_sites(const std::vector<class_mps_site> & new_mps){
    for(const auto & mps : new_mps)
        get_mps_site(mps.get_position()) = mps;
}


std::array<long, 3> class_state_finite::active_dimensions() const { return tools::finite::multisite::get_dimensions(*this, active_sites); }

long class_state_finite::active_problem_size() const { return tools::finite::multisite::get_problem_size(*this, active_sites); }

std::vector<long> class_state_finite::get_spin_dims(const std::vector<size_t> & sites) const{
    if(sites.empty()) throw std::runtime_error("No sites on which to collect spin dimensions");
    std::vector<long> dims;
    dims.reserve(sites.size());
    for(const auto & site : sites ){
        dims.emplace_back(get_mps_site(site).spin_dim());
    }
    return dims;
}

std::vector<long> class_state_finite::get_spin_dims() const{
    return get_spin_dims(active_sites);
}

Eigen::Tensor<class_state_finite::Scalar, 3> class_state_finite::get_multisite_mps(const std::vector<size_t> &sites) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mps tensor");
    if(sites == active_sites and cache.multisite_mps) return cache.multisite_mps.value();
    if constexpr(settings::debug) tools::log->trace("Contracting multisite mps tensor with sites {}", sites);
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
        multisite_tensor                      = temp;
    }
    if(sites.front() != 0 and get_mps_site(sites.front()).get_label() == "B") {
        // In this case all sites are "B" and we need to prepend the the "L" from the site on the left to make a normalized multisite tensor
        auto &mps_left = get_mps_site(sites.front() - 1);
        auto &L_left   = mps_left.isCenter() ? mps_left.get_LC() : mps_left.get_L();
        if(L_left.dimension(0) != multisite_tensor.dimension(1))
            throw std::logic_error(fmt::format("Mismatching dimensions: L_left {} | multisite_tensor {}", L_left.dimensions(), multisite_tensor.dimensions()));
        temp.resize(multisite_tensor.dimension(0), L_left.dimension(0), multisite_tensor.dimension(2));
        temp.device(Textra::omp::getDevice()) = Textra::asDiagonal(L_left).contract(multisite_tensor, Textra::idx({1}, {1})).shuffle(Textra::array3{1, 0, 2});
        multisite_tensor                      = temp;
    } else if(sites.back() < get_length() - 1 and get_mps_site(sites.back()).get_label() == "A") {
        // In this case all sites are "A" and we need to append the the "L" from the site on the right to make a normalized multisite tensor
        auto &mps_right = get_mps_site(sites.back() + 1);
        auto &L_right   = mps_right.get_L();
        if(L_right.dimension(0) != multisite_tensor.dimension(2))
            throw std::logic_error(
                fmt::format("Mismatching dimensions: L_right {} | multisite_tensor {}", L_right.dimensions(), multisite_tensor.dimensions()));
        temp.resize(multisite_tensor.dimension(0), multisite_tensor.dimension(1), L_right.dimension(0));
        temp.device(Textra::omp::getDevice()) = multisite_tensor.contract(Textra::asDiagonal(L_right), Textra::idx({2}, {1}));
        multisite_tensor                      = temp;
    }

    tools::common::profile::get_default_prof()["t_mps"]->toc();
    if constexpr(settings::debug) {
        // Check the norm of the tensor on debug builds
        tools::common::profile::get_default_prof()["t_dbg"]->tic();
        double norm = Textra::norm(multisite_tensor.contract(multisite_tensor.conjugate(), Textra::idx({0, 1, 2}, {0, 1, 2})));
        if(std::abs(norm - 1) > settings::precision::max_norm_error) {
            for(const auto &site : sites) {
                auto &mps = get_mps_site(site);
                auto &M   = mps.get_M();
                tools::log->critical("{}({}) norm: {:.16f}", mps.get_label(), site, Textra::VectorMap(M).norm());
            }
            if(sites.front() != 0 and get_mps_site(sites.front()).get_label() == "B") {
                // In this case all sites are "B" and we need to prepend the the "L" from the site on the left
                auto &mps_left   = get_mps_site(sites.front() - 1);
                auto &L_left     = mps_left.isCenter() ? mps_left.get_LC() : mps_left.get_L();
                temp             = Textra::asDiagonal(L_left).contract(multisite_tensor, Textra::idx({1}, {1})).shuffle(Textra::array3{1, 0, 2});
                multisite_tensor = temp;
                norm             = Textra::norm(multisite_tensor.contract(multisite_tensor.conjugate(), Textra::idx({0, 1, 2}, {0, 1, 2})));
                tools::log->critical("Norm after adding L to B from the left: {:.16f}", norm);
            }

            throw std::runtime_error(fmt::format("Multisite tensor for sites {} is not normalized. Norm = {:.16f}", sites, norm));
        }
        tools::common::profile::get_default_prof()["t_dbg"]->toc();
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

double class_state_finite::get_truncation_error() const {
    auto pos = get_position<long>();
    if(pos >= 0)
        return get_mps_site(pos).get_truncation_error();
    else
        return 0;
}

double class_state_finite::get_truncation_error_LC() const { return get_mps_site(get_position()).get_truncation_error_LC(); }
double class_state_finite::get_truncation_error_midchain() const {
    auto pos = get_position<long>();
    auto cnt = (get_length<long>() - 1) / 2;
    if(pos < cnt) return get_mps_site(cnt).get_truncation_error();
    if(pos > cnt) return get_mps_site(cnt + 1).get_truncation_error();
    return get_mps_site(cnt).get_truncation_error_LC();
}

std::vector<double> class_state_finite::get_truncation_errors() const { return tools::finite::measure::truncation_errors(*this); }
std::vector<double> class_state_finite::get_truncation_errors_active() const {
    std::vector<double> truncation_errors;
    truncation_errors.reserve(active_sites.size());
    for(const auto &pos : active_sites) {
        // We are only interested in the truncation on bonds that are updated
        // when operating on active_sites. This excludes the outer bonds.
        if(get_mps_site(pos).isCenter()) truncation_errors.emplace_back(get_truncation_error_LC());
        if(pos == active_sites.front()) continue;
        if(pos == active_sites.back()) continue;
        truncation_errors.emplace_back(get_truncation_error(pos));
    }
    return truncation_errors;
}

size_t class_state_finite::num_sites_truncated(double truncation_threshold) const {
    auto truncation_errors = get_truncation_errors();
    auto trunc_bond_count  = static_cast<size_t>(
        std::count_if(truncation_errors.begin(), truncation_errors.end(), [truncation_threshold](auto const &val) { return val > truncation_threshold; }));
    return trunc_bond_count;
}

size_t class_state_finite::num_bonds_reached_chi(long chi_level) const {
    auto bond_dims    = tools::finite::measure::bond_dimensions(*this);
    auto bonds_at_lim = static_cast<size_t>(std::count_if(bond_dims.begin(), bond_dims.end(), [chi_level](auto const &val) { return val >= chi_level; }));
    return bonds_at_lim;
}

bool class_state_finite::is_bond_limited(long chi_lim, double truncation_threshold) const {
    return num_sites_truncated(truncation_threshold) > 0 and num_bonds_reached_chi(chi_lim) > 0;
}

void class_state_finite::clear_measurements(LogPolicy logPolicy) const {
    if(logPolicy == LogPolicy::NORMAL) tools::log->trace("Clearing state measurements");
    measurements = state_measure_finite();
}

void class_state_finite::clear_cache(LogPolicy logPolicy) const {
    if(logPolicy == LogPolicy::NORMAL) tools::log->trace("Clearing state cache");
    cache = Cache();
}

void class_state_finite::do_all_measurements() const { tools::finite::measure::do_all_measurements(*this); }

void class_state_finite::tag_active_sites_normalized(bool tag) const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot tag active sites, size mismatch in site list");
    for(auto &site : active_sites) tag_normalized_sites[site] = tag;
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
    tools::log->trace("Checking normalization status on all sites", active_sites);
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
    tools::log->trace("Checking normalization status on active sites: {}", active_sites);
    bool normalized = std::all_of(first_site_ptr, last_site_ptr, [](bool v) { return v; });
    tools::log->trace("Active sites normalized: {}", normalized);
    return normalized;
}

bool class_state_finite::is_normalized_on_non_active_sites() const {
    if(tag_normalized_sites.size() != get_length()) throw std::runtime_error("Cannot check update status on all sites, size mismatch in site list");
    if(active_sites.empty()) return is_normalized_on_all_sites();
    tools::log->trace("Checking normalization status on non-active sites", active_sites);
    for(size_t idx = 0; idx < get_length(); idx++)
        if(std::find(active_sites.begin(), active_sites.end(), idx) == active_sites.end() and not tag_normalized_sites[idx]) return false;
    return true;
}

std::vector<size_t> class_state_finite::get_active_ids() const {
    std::vector<size_t> ids;
    ids.reserve(active_sites.size());
    for(const auto &pos : active_sites) ids.emplace_back(get_mps_site(pos).get_unique_id());
    return ids;
}
