#include "math/tenx.h"

// -- (textra first)
#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "math/cast.h"
#include "math/num.h"
#include "StateFinite.h"
#include "tensors/site/mps/MpsSite.h"
#include "tid/tid.h"
#include "tools/common/contraction.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include "tools/finite/multisite.h"
#include <fmt/ranges.h>

namespace settings {
    inline constexpr bool debug_state = false;
}

StateFinite::StateFinite() = default; // Can't initialize lists since we don't know the model size yet

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise, unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
StateFinite::~StateFinite() noexcept                               = default; // default dtor
StateFinite:: StateFinite(StateFinite &&other) noexcept            = default; // default move ctor
StateFinite  &StateFinite::operator=(StateFinite &&other) noexcept = default; // default move assign

/* clang-format off */
StateFinite::StateFinite(const StateFinite &other):
    direction(other.direction),
    cache(other.cache),
    tag_normalized_sites(other.tag_normalized_sites),
    name(other.name),
    algo(other.algo),
    active_sites(other.active_sites),
    measurements(other.measurements),
    popcount(other.popcount)
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
        popcount             = other.popcount;
        mps_sites.clear();
        mps_sites.reserve(other.mps_sites.size());
        for(const auto &mps : other.mps_sites) mps_sites.emplace_back(std::make_unique<MpsSite>(*mps));
    }
    return *this;
}

StateFinite::StateFinite(AlgorithmType algo_type, size_t model_size, long position, long spin_dim) { initialize(algo_type, model_size, position, spin_dim); }

void StateFinite::initialize(AlgorithmType algo_type, size_t model_size, long position, long spin_dim) {
    set_algorithm(algo_type);
    tools::log->debug("Initializing state: sites {} | position {} | spin_dim {}", model_size, position, spin_dim);
    if(model_size < 2) throw except::logic_error("Tried to initialize state with less than 2 sites");
    if(model_size > 2048) throw except::logic_error("Tried to initialize state with more than 2048 sites");
    if(position >= safe_cast<long>(model_size)) throw except::logic_error("Tried to initialize state at a position larger than the number of sites");

    mps_sites.clear();

    // Generate a simple state with all spins equal
    Eigen::Tensor<cplx, 3> M(safe_cast<long>(spin_dim), 1, 1);
    Eigen::Tensor<cplx, 1> L(1);
    M(0, 0, 0) = 0;
    M(1, 0, 0) = 1;
    L(0)       = 1;
    for(size_t site = 0; site < model_size; site++) {
        std::string label = safe_cast<long>(site) <= position ? "A" : "B";
        mps_sites.emplace_back(std::make_unique<MpsSite>(M, L, site, 0.0, label));
        if(safe_cast<long>(site) == position) { mps_sites.back()->set_LC(L); }
    }
    if(mps_sites.size() != model_size) throw except::logic_error("Initialized state with wrong size");
    if(not get_mps_site(position).isCenter()) throw except::logic_error("Initialized state center bond at the wrong position");
    if(get_position<long>() != position) throw except::logic_error("Initialized state at the wrong position");
    tag_normalized_sites = std::vector<bool>(model_size, false);
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
    return safe_cast<T>(mps_sites.size());
}
template double             StateFinite::get_length<double>() const;
template size_t             StateFinite::get_length<size_t>() const;
template long               StateFinite::get_length<long>() const;
template int                StateFinite::get_length<int>() const;
template unsigned long long StateFinite::get_length<unsigned long long>() const; // hsize_t from hdf5

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
template size_t             StateFinite::get_position<size_t>() const;
template long               StateFinite::get_position<long>() const;
template int                StateFinite::get_position<int>() const;
template unsigned long long StateFinite::get_position<unsigned long long>() const; // hsize_t from hdf5

long StateFinite::get_largest_bond() const {
    auto bond_dimensions = tools::finite::measure::bond_dimensions(*this);
    return *max_element(std::begin(bond_dimensions), std::end(bond_dimensions));
}

double StateFinite::get_smallest_schmidt_value() const {
    double schmidt_min = 1;
    for(const auto &mps : mps_sites) {
        const auto &L = mps->get_L();
        schmidt_min   = std::min(schmidt_min, L.coeff(L.size() - 1).real());
        if(mps->isCenter()) {
            const auto &LC = mps->get_LC();
            schmidt_min    = std::min(schmidt_min, LC.coeff(LC.size() - 1).real());
        }
    }
    return schmidt_min;
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
    return get_position<long>() >= get_length<long>() - safe_cast<long>(nsite) and direction == 1;
}

bool StateFinite::position_is_outward_edge(size_t nsite) const { return position_is_outward_edge_left(nsite) or position_is_outward_edge_right(nsite); }

bool StateFinite::position_is_inward_edge_left([[maybe_unused]] size_t nsite) const {
    return get_position<long>() == 0 and direction == 1; // i.e. first site is an AC going to the right
}

bool StateFinite::position_is_inward_edge_right(size_t nsite) const {
    return get_position<long>() >= get_length<long>() - safe_cast<long>(nsite) and direction == -1;
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
            throw except::runtime_error("State is corrupted: position mismatch: expected position {} != mps position {}", pos, mps->get_position());
        pos++;
    }

    for(const auto &mps : mps_sites) mps->assert_validity();
    if(settings::model::model_type == ModelType::ising_sdual or settings::model::model_type == ModelType::ising_majorana) {
        for(const auto &mps : mps_sites)
            if(not mps->is_real()) throw except::runtime_error("state has imaginary part at mps position {}", mps->get_position());
    }
}

const Eigen::Tensor<cplx, 1> &StateFinite::get_bond(long posL, long posR) const {
    if(posL + 1 != posR) throw except::runtime_error("Expected posL+1 == posR, got: posL {}, posR {}", posL, posR);
    auto pos = get_position<long>();
    if(pos < posL) return get_mps_site(posL).get_L(); // B.B
    if(pos > posL) return get_mps_site(posR).get_L(); // A.A or A.AC
    return get_mps_site(posL).get_LC();               // AC.B
}

const Eigen::Tensor<cplx, 1> &StateFinite::get_midchain_bond() const {
    auto pos = get_position<long>();
    auto cnt = (get_length<long>() - 1) / 2;
    if(pos < cnt) return get_mps_site(cnt).get_L();
    if(pos > cnt) return get_mps_site(cnt + 1).get_L();
    return get_mps_site(cnt).get_LC();
}

const Eigen::Tensor<cplx, 1> &StateFinite::current_bond() const { return get_mps_site(get_position()).get_LC(); }

template<typename T>
const MpsSite &StateFinite::get_mps_site(T pos) const {
    if constexpr(std::is_signed_v<T>)
        if(pos < 0) throw except::range_error("get_mps_site(pos): pos out of range: {}", pos);
    if(pos >= get_length<T>()) throw except::range_error("get_mps_site(pos): pos out of range: {}", pos);
    const auto &mps_ptr = *std::next(mps_sites.begin(), safe_cast<long>(pos));
    if(mps_ptr->template get_position<T>() != pos)
        throw except::range_error("get_mps_site(pos): mismatch pos {} != mps pos {}", pos, mps_ptr->template get_position<T>());
    return *mps_ptr;

    //    if(algo == AlgorithmType::fLBIT){
    //        // During fLBIT we can't assume the mps positions are in order since we could have swap operators.
    //        for(const auto & mps_ptr :  mps_sites ){
    //            if(mps_ptr->get_position<T>() == pos) return *mps_ptr;
    //        }
    //        throw except::runtime_error("get_mps_site(pos): pos {} not found", pos);
    //    }else{
    //        // There shouldn't be any swap operator, we can safely assume the mps positions are sorted
    //        const auto &mps_ptr = *std::next(mps_sites.begin(), safe_cast<long>(pos));
    //        if(mps_ptr->get_position<T>() != pos)
    //            throw except::range_error("get_mps_site(pos): mismatch pos {} != mps pos {}", pos, mps_ptr->get_position<T>());
    //        return *mps_ptr;
    //    }
}
template const MpsSite &StateFinite::get_mps_site(size_t pos) const;
template const MpsSite &StateFinite::get_mps_site(long pos) const;
template const MpsSite &StateFinite::get_mps_site(int pos) const;
template const MpsSite &StateFinite::get_mps_site(unsigned long long pos) const; // hsize_t

template<typename T>
MpsSite &StateFinite::get_mps_site(T pos) {
    return const_cast<MpsSite &>(std::as_const(*this).get_mps_site<T>(pos));
}
template MpsSite &StateFinite::get_mps_site(size_t pos);
template MpsSite &StateFinite::get_mps_site(long pos);
template MpsSite &StateFinite::get_mps_site(int pos);
template MpsSite &StateFinite::get_mps_site(unsigned long long pos); // hsize_t

const MpsSite &StateFinite::get_mps_site() const { return get_mps_site(get_position()); }

MpsSite &StateFinite::get_mps_site() { return get_mps_site(get_position()); }

void StateFinite::set_mps(const std::vector<MpsSite> &mps_list) {
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

std::vector<std::reference_wrapper<const MpsSite>> StateFinite::get_mps(const std::vector<size_t> &sites) const {
    std::vector<std::reference_wrapper<const MpsSite>> mps;
    mps.reserve(sites.size());
    for(auto &site : sites) mps.emplace_back(get_mps_site(site));
    return mps;
}

std::vector<std::reference_wrapper<MpsSite>> StateFinite::get_mps(const std::vector<size_t> &sites) {
    std::vector<std::reference_wrapper<MpsSite>> mps;
    mps.reserve(sites.size());
    for(auto &site : sites) mps.emplace_back(get_mps_site(site));
    return mps;
}
std::vector<MpsSite> StateFinite::get_mps_copy(const std::vector<size_t> &sites) {
    std::vector<MpsSite> mps;
    mps.reserve(sites.size());
    for(auto &site : sites) mps.emplace_back(get_mps_site(site));
    return mps;
}

std::vector<std::reference_wrapper<const MpsSite>> StateFinite::get_mps_active() const { return get_mps(active_sites); }
std::vector<std::reference_wrapper<MpsSite>>       StateFinite::get_mps_active() { return get_mps(active_sites); }

std::array<long, 3> StateFinite::active_dimensions() const { return tools::finite::multisite::get_dimensions(*this, active_sites); }

long StateFinite::active_problem_size() const { return tools::finite::multisite::get_problem_size(*this, active_sites); }

std::vector<long> StateFinite::get_bond_dims(const std::vector<size_t> &sites) const {
    // If the sites are {2,3,4,5,6} this returns the 4 bonds connecting {2,3}, {3,4}, {4,5} and {5,6}
    // If sites is just {4}, it returns the bond between {4,5} when going left or right.
    auto t_chi = tid::tic_scope("bond_merged", tid::level::highest);
    if(sites.empty()) return {};
    if(sites.size() == 1) {
        // In single-site DMRG the active site is a center "AC" site:
        //  * Going left-to-right, the forward (right) bond is expanded, and this same bond is truncated when merging
        //  * Going right-to-left, the forward (left) bond is expanded (L), but LC is still the one truncated when merging.
        return {get_mps_site(sites.front()).get_chiR()};
    }
    if(sites.size() == 2) return {get_mps_site(sites.front()).get_chiR()};
    std::vector<long> bond_dimensions;
    for(const auto &pos : sites) {
        if(&pos == &sites.front()) continue;
        const auto &mps = get_mps_site(pos);
        bond_dimensions.push_back(mps.get_chiL());
    }
    return bond_dimensions;
}
std::vector<long> StateFinite::get_bond_dims_active() const { return get_bond_dims(active_sites); }

std::vector<long> StateFinite::get_spin_dims(const std::vector<size_t> &sites) const {
    if(sites.empty()) throw except::runtime_error("No sites on which to collect spin dimensions");
    std::vector<long> dims;
    dims.reserve(sites.size());
    for(const auto &site : sites) { dims.emplace_back(get_mps_site(site).spin_dim()); }
    return dims;
}

std::vector<long> StateFinite::get_spin_dims() const { return get_spin_dims(active_sites); }
long              StateFinite::get_spin_dim() const { return get_mps_site(0).spin_dim(); }

std::vector<std::array<long, 3>> StateFinite::get_mps_dims(const std::vector<size_t> &sites) const {
    std::vector<std::array<long, 3>> dims;
    for(const auto &pos : sites) dims.emplace_back(get_mps_site(pos).dimensions());
    return dims;
}

std::vector<std::array<long, 3>> StateFinite::get_mps_dims_active() const { return get_mps_dims(active_sites); }

template<typename Scalar>
Eigen::Tensor<Scalar, 3> StateFinite::get_multisite_mps(const std::vector<size_t> &sites, bool use_cache) const {
    if(sites.empty()) throw except::runtime_error("No active sites on which to build a multisite mps tensor");
    if constexpr(std::is_same_v<Scalar, cplx>) {
        if(sites == active_sites and cache.multisite_mps) { return cache.multisite_mps.value(); }
        auto                   t_mps  = tid::tic_scope("gen_mps", tid::level::highest);
        auto                   csites = std::vector<size_t>{}; // Keeps track of the contracted sites
        auto                   length = get_length<size_t>();
        Eigen::Tensor<cplx, 3> multisite_mps;
        Eigen::Tensor<cplx, 3> temp;
        for(auto &site : sites) {
            csites.emplace_back(site);
            const auto mps_key = fmt::format("{}", csites);
            if(use_cache and cache.mps_cplx.find(mps_key) != cache.mps_cplx.end()) {
                tools::log->debug("multisite_mps: cache_hit: {}", csites);
                multisite_mps = cache.mps_cplx.at(mps_key);
            } else {
                const auto &mps       = get_mps_site(site);
                auto        M         = mps.get_M();
                bool        prepend_L = mps.get_label() == "B" and site > 0 and site == sites.front();
                bool        append_L  = mps.get_label() == "A" and site + 1 < length and site == sites.back();
                if(prepend_L) {
                    // In this case all sites are "B" and we need to prepend the "L" from the site on the left to make a normalized multisite mps
                    const auto &mps_left = get_mps_site(site - 1);
                    const auto &L        = mps_left.isCenter() ? mps_left.get_LC() : mps_left.get_L();
                    if(L.dimension(0) != M.dimension(1))
                        throw except::logic_error("get_multisite_mps<cplx>({}): mismatching dimensions: L (left) {} | M {}", sites, L.dimensions(),
                                                  M.dimensions());
                    M = tools::common::contraction::contract_bnd_mps_temp(L, M, temp);
                }
                if(append_L) {
                    // In this case all sites are "A" and we need to append the "L" from the site on the right to make a normalized multisite mps
                    const auto &mps_right = get_mps_site(site + 1);
                    const auto &L         = mps_right.get_L();
                    if(L.dimension(0) != M.dimension(2))
                        throw except::logic_error("get_multisite_mps<cplx>({}): mismatching dimensions: M {} | L (right) {}", sites, M.dimensions(),
                                                  L.dimensions());
                    M = tools::common::contraction::contract_mps_bnd_temp(M, L, temp);
                }

                if(&site == &sites.front()) { // First site
                    multisite_mps = std::move(M);
                } else { // Next sites
                    multisite_mps = tools::common::contraction::contract_mps_mps_temp(multisite_mps, M, temp);
                }

                if(use_cache and not append_L) // If it is the last site, we may have closed off by appending L
                    cache.mps_cplx[mps_key] = multisite_mps;
            }
        }
        if constexpr(settings::debug) {
            // Check the norm of the tensor on debug builds
            auto t_dbg = tid::tic_scope("debug");
            cplx norm  = tools::common::contraction::contract_mps_norm(multisite_mps);
            if constexpr(settings::debug_state) tools::log->trace("get_multisite_mps({}): norm ⟨ψ|ψ⟩ = {:.16f}", sites, norm);
            if(std::abs(norm - 1.0) > settings::precision::max_norm_error) {
                tools::log->warn("get_multisite_mps<cplx>({}): norm error |1-⟨ψ|ψ⟩| = {:.2e} > max_norm_error {:.2e}", sites, std::abs(norm - 1.0),
                                 settings::precision::max_norm_error);
                //                throw except::runtime_error("get_multisite_mps<cplx>({}): norm error |1-⟨ψ|ψ⟩| = {:.2e} > max_norm_error {:.2e}", sites,
                //                std::abs(norm - 1),
                //                                            settings::precision::max_norm_error);
            }
        }
        return multisite_mps;
    } else if constexpr(std::is_same_v<Scalar, real>) {
        auto                   t_mps  = tid::tic_scope("gen_mps", tid::level::highest);
        auto                   csites = std::vector<size_t>{}; // Keeps track of the contracted sites
        auto                   length = get_length<size_t>();
        Eigen::Tensor<real, 3> multisite_mps;
        Eigen::Tensor<real, 3> temp;
        for(auto &site : sites) {
            csites.emplace_back(site);
            const auto mps_key = fmt::format("{}", csites);
            if(use_cache and cache.mps_real.find(mps_key) != cache.mps_real.end()) {
                tools::log->trace("multisite_mps: cache_hit: {}", csites);
                multisite_mps = cache.mps_real.at(mps_key);
            } else {
                const auto &mps       = get_mps_site(site);
                auto        M         = Eigen::Tensor<real, 3>(mps.get_M().real());
                bool        prepend_L = mps.get_label() == "B" and site > 0 and site == sites.front();
                bool        append_L  = mps.get_label() == "A" and site + 1 < length and site == sites.back();
                if(prepend_L) {
                    // In this case all sites are "B" and we need to prepend the "L" from the site on the left to make a normalized multisite mps
                    tools::log->trace("Prepending L to B site {}", site);
                    auto        t_prepend = tid::tic_scope("prepend", tid::level::higher);
                    const auto &mps_left  = get_mps_site(site - 1);
                    const auto  L         = Eigen::Tensor<real, 1>(mps_left.isCenter() ? mps_left.get_LC().real() : mps_left.get_L().real());
                    if(L.dimension(0) != M.dimension(1))
                        throw except::logic_error("get_multisite_mps<real>: mismatching dimensions: L (left) {} | M {}", L.dimensions(), M.dimensions());
                    M = tools::common::contraction::contract_bnd_mps_temp(L, M, temp);
                }
                if(append_L) {
                    // In this case all sites are "A" and we need to append the "L" from the site on the right to make a normalized multisite mps
                    tools::log->trace("Appending L to A site {}", site);
                    auto        t_append  = tid::tic_scope("append", tid::level::higher);
                    const auto &mps_right = get_mps_site(site + 1);
                    const auto &L         = Eigen::Tensor<real, 1>(mps_right.get_L().real());
                    if(L.dimension(0) != M.dimension(2))
                        throw except::logic_error("get_multisite_mps<real>: mismatching dimensions: M {} | L (right) {}", M.dimensions(), L.dimensions());
                    M = tools::common::contraction::contract_mps_bnd_temp(M, L, temp);
                }

                if(&site == &sites.front()) { // First site
                    multisite_mps = std::move(M);
                } else { // Next sites
                    multisite_mps = tools::common::contraction::contract_mps_mps_temp(multisite_mps, M, temp);
                }
                if(use_cache and not append_L) // If it is the last site, we may have closed off by appending L
                    cache.mps_real[mps_key] = multisite_mps;
            }
        }
        if constexpr(settings::debug) {
            // Check the norm of the tensor on debug builds
            auto   t_dbg = tid::tic_scope("debug");
            double norm  = tools::common::contraction::contract_mps_norm(multisite_mps);
            if constexpr(settings::debug_state) tools::log->trace("get_multisite_mps({}): norm ⟨ψ|ψ⟩ = {:.16f}", sites, norm);
            if(std::abs(norm - 1) > settings::precision::max_norm_error) {
                tools::log->warn("get_multisite_mps<real>({}): norm error |1-⟨ψ|ψ⟩| = {:.2e} > max_norm_error {:.2e}", sites, std::abs(norm - 1),
                                 settings::precision::max_norm_error);
                //                throw except::runtime_error("get_multisite_mps<real>({}): norm error |1-⟨ψ|ψ⟩| = {:.2e} > max_norm_error {:.2e}", sites,
                //                std::abs(norm - 1),
                //                                            settings::precision::max_norm_error);
            }
        }
        return multisite_mps;
    }
}

template Eigen::Tensor<real, 3> StateFinite::get_multisite_mps<real>(const std::vector<size_t> &sites, bool use_cache) const;
template Eigen::Tensor<cplx, 3> StateFinite::get_multisite_mps<cplx>(const std::vector<size_t> &sites, bool use_cache) const;

const Eigen::Tensor<cplx, 3> &StateFinite::get_multisite_mps() const {
    if(cache.multisite_mps) return cache.multisite_mps.value();
    cache.multisite_mps = get_multisite_mps(active_sites);
    return cache.multisite_mps.value();
}

template<typename Scalar>
Eigen::Tensor<Scalar, 2> StateFinite::get_reduced_density_matrix(const std::vector<size_t> &sites) const {
    auto t_rho = tid::tic_scope("rho");
    auto asites = num::range<size_t>(sites.front(), sites.back() + 1); // Contiguous list of all sites
    if(sites == asites) {
        // We have a contiguous set
        auto mps = get_multisite_mps<Scalar>(sites, true);
        return tools::common::contraction::contract_mps_partial<std::array{1l, 2l}>(mps);
    } else {
        // Note that for discontiguous sets, we expect the front and back sites to be at the edges!
        auto &threads   = tenx::threads::get();
        auto  rho_temp  = Eigen::Tensor<Scalar, 4>(); // Will accumulate the sites
        auto  rho_temp2 = Eigen::Tensor<Scalar, 4>();
        auto  M         = Eigen::Tensor<Scalar, 3>();
        for(size_t i = sites.front(); i <= sites.back(); ++i) {
            const auto &mps = get_mps_site(i);
            if(i == sites.front()) {
                if constexpr(std::is_same_v<Scalar, real>) {
                    M = mps.get_label() == "B" ? get_multisite_mps({sites.front()}).real()
                                               : mps.get_M().real(); // Could be an A, AC or B. Either way we need it to include the left schmidt values
                } else {
                    M = mps.get_label() == "B" ? get_multisite_mps({sites.front()})
                                               : mps.get_M(); // Could be an A, AC or B. Either way we need it to include the left schmidt values
                }
                auto dim = M.dimensions();
                rho_temp.resize(std::array{dim[0], dim[0], dim[2], dim[2]});
                rho_temp.device(*threads->dev) = M.conjugate().contract(M, tenx::idx({1}, {1})).shuffle(std::array{0, 2, 1, 3});
            } else {
                if constexpr(std::is_same_v<Scalar, real>) {
                    M = (i == sites.back() and mps.get_label() == "A") ? get_multisite_mps({sites.back()}).real() : mps.get_M().real();
                    M = mps.get_M().real();
                } else {
                    M = (i == sites.back() and mps.get_label() == "A") ? get_multisite_mps({sites.back()}) : mps.get_M();
                }
                auto mps_dim  = M.dimensions();
                auto rho_dim  = rho_temp.dimensions();
                bool do_trace = std::find(sites.begin(), sites.end(), i) == sites.end();
                if(do_trace) {
                    auto new_dim = std::array{rho_dim[0], rho_dim[1], mps_dim[2], mps_dim[2]};
                    rho_temp2.resize(new_dim);
                    rho_temp2.device(*threads->dev) = rho_temp.contract(M.conjugate(), tenx::idx({2}, {1})).contract(M, tenx::idx({2, 3}, {1, 0}));
                } else {
                    auto new_dim = std::array{rho_dim[0] * mps_dim[0], rho_dim[1] * mps_dim[0], mps_dim[2], mps_dim[2]};
                    rho_temp2.resize(new_dim);
                    rho_temp2.device(*threads->dev) = rho_temp.contract(M.conjugate(), tenx::idx({2}, {1}))
                                                          .contract(M, tenx::idx({2}, {1}))
                                                          .shuffle(std::array{0, 2, 1, 4, 3, 5})
                                                          .reshape(new_dim);
                }
                rho_temp = std::move(rho_temp2);
            }
        }
        return rho_temp.trace(std::array{2, 3});
    }
}

template Eigen::Tensor<real, 2> StateFinite::get_reduced_density_matrix<real>(const std::vector<size_t> &sites) const;
template Eigen::Tensor<cplx, 2> StateFinite::get_reduced_density_matrix<cplx>(const std::vector<size_t> &sites) const;

void StateFinite::set_truncation_error(size_t pos, double error) { get_mps_site(pos).set_truncation_error(error); }
void StateFinite::set_truncation_error(double error) { set_truncation_error(get_position(), error); }

void StateFinite::set_truncation_error_LC(double error) {
    auto &mps = get_mps_site(get_position());
    if(not mps.isCenter()) throw except::runtime_error("mps at current position is not a center");
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
    auto trunc_bond_count  = safe_cast<size_t>(
        std::count_if(truncation_errors.begin(), truncation_errors.end(), [truncation_threshold](auto const &val) { return val > truncation_threshold; }));
    return trunc_bond_count;
}

size_t StateFinite::num_bonds_at_limit(long bond_lim) const {
    auto bond_dimensions = tools::finite::measure::bond_dimensions(*this);
    auto bonds_at_lim =
        safe_cast<size_t>(std::count_if(bond_dimensions.begin(), bond_dimensions.end(), [bond_lim](auto const &dim) { return dim >= bond_lim; }));
    return bonds_at_lim;
}

bool StateFinite::is_limited_by_bond(long bond_lim) const { return num_bonds_at_limit(bond_lim) > 0; }

bool StateFinite::is_truncated(double truncation_error_limit) const {
    auto truncation_errors = get_truncation_errors();
    auto num_above_lim     = static_cast<size_t>(
        std::count_if(truncation_errors.begin(), truncation_errors.end(), [truncation_error_limit](auto const &err) { return err >= truncation_error_limit; }));
    return num_above_lim > 0;
}

void StateFinite::clear_measurements(LogPolicy logPolicy) const {
    if(logPolicy == LogPolicy::NORMAL) tools::log->trace("Clearing state measurements");
    measurements = MeasurementsStateFinite();
}

void StateFinite::clear_cache(LogPolicy logPolicy) const {
    if(logPolicy == LogPolicy::NORMAL) tools::log->trace("Clearing state cache");
    cache = Cache();
}

void StateFinite::tag_active_sites_normalized(bool tag) const {
    if(tag_normalized_sites.size() != get_length()) throw except::runtime_error("Cannot tag active sites, size mismatch in site list");
    for(auto &site : active_sites) tag_normalized_sites[site] = tag;
}

void StateFinite::tag_all_sites_normalized(bool tag) const {
    if(tag_normalized_sites.size() != get_length()) throw except::runtime_error("Cannot untag all sites, size mismatch in site list");
    tag_normalized_sites = std::vector<bool>(get_length(), tag);
}

void StateFinite::tag_site_normalized(size_t pos, bool tag) const {
    if(tag_normalized_sites.size() != get_length()) throw except::runtime_error("Cannot untag all sites, size mismatch in site list");
    tag_normalized_sites[pos] = tag;
}

bool StateFinite::is_normalized_on_all_sites() const {
    if(tag_normalized_sites.size() != get_length()) throw except::runtime_error("Cannot check normalization status on all sites, size mismatch in site list");
    // If all tags are false then we should definitely normalize:
    auto normalized_none = std::none_of(tag_normalized_sites.begin(), tag_normalized_sites.end(), [](bool v) { return v; });
    if(normalized_none) {
        tools::log->debug("{} normalized: false (none)", get_name());
        return false;
    }

    if constexpr(settings::debug) {
        auto normalized_some = std::any_of(tag_normalized_sites.begin(), tag_normalized_sites.end(), [](bool v) { return v; });
        if(normalized_some) {
            // In debug mode we check if the tags are truthful
            for(const auto &mps : mps_sites) {
                auto pos = mps->get_position<size_t>();
                if(not tag_normalized_sites[pos]) {
                    if(mps->is_normalized(settings::precision::max_norm_error)) tag_normalized_sites[pos] = true;
                }
            }
        }
    }
    auto normalized_tags = std::all_of(tag_normalized_sites.begin(), tag_normalized_sites.end(), [](bool v) { return v; });
    auto normalized_fast = true;
    auto normalized_full = true;
    auto normalized_site = true;
    auto msg             = fmt::format("tags {}", normalized_tags);
    if(normalized_tags) {
        // We don't need this check if the tags already told us the state isn't normalized
        auto norm       = tools::finite::measure::norm(*this, false);
        normalized_fast = std::abs(norm - 1.0) <= settings::precision::max_norm_error;
        msg += fmt::format(" | fast {} {:.3e}", normalized_fast, norm);
    }
    if constexpr(settings::debug) {
        if(normalized_tags and normalized_fast) {
            auto norm       = tools::finite::measure::norm(*this, true);
            normalized_full = std::abs(norm - 1.0) <= settings::precision::max_norm_error;
            msg += fmt::format(" | full {} {:.3e}", normalized_full, norm);
        }
        if(normalized_tags and normalized_fast and normalized_full) {
            std::vector<long> site_list;
            for(const auto &mps : mps_sites) {
                if(not mps->is_normalized(settings::precision::max_norm_error)) { site_list.emplace_back(mps->get_position<long>()); }
            }
            if(not site_list.empty()) {
                normalized_site = false;
                msg += fmt::format(" | non-normalized sites {}", site_list);
            }
        }
    }
    tools::log->debug("{} normalized: {}", get_name(), msg);
    return normalized_tags and normalized_fast and normalized_full and normalized_site;
}

bool StateFinite::is_normalized_on_any_sites() const {
    if(tag_normalized_sites.size() != get_length()) throw except::runtime_error("Cannot check normalization status on any sites, size mismatch in site list");
    return std::any_of(tag_normalized_sites.begin(), tag_normalized_sites.end(), [](bool v) { return v; });
}

bool StateFinite::is_normalized_on_active_sites() const {
    if(tag_normalized_sites.size() != get_length())
        throw except::runtime_error("Cannot check normalization status on active sites, size mismatch in site list");
    if(active_sites.empty()) return false;
    auto first_site_ptr = std::next(tag_normalized_sites.begin(), safe_cast<long>(active_sites.front()));
    auto last_site_ptr  = std::next(tag_normalized_sites.begin(), safe_cast<long>(active_sites.back()));
    return std::all_of(first_site_ptr, last_site_ptr, [](bool v) { return v; });
}

bool StateFinite::is_normalized_on_non_active_sites() const {
    if(tag_normalized_sites.size() != get_length()) throw except::runtime_error("Cannot check update status on all sites, size mismatch in site list");
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
const std::vector<bool> &StateFinite::get_normalization_tags() const { return tag_normalized_sites; }
