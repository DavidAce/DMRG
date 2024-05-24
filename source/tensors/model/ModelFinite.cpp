#include "math/tenx.h"
// -- (textra first)
#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "math/cast.h"
#include "math/eig.h"
#include "math/linalg/tensor.h"
#include "math/svd.h"
#include "ModelFinite.h"
#include "ModelLocal.h"
#include "qm/spin.h"
#include "tensors/site/mpo/MpoFactory.h"
#include "tid/tid.h"
#include "tools/finite/mpo.h"
#include "tools/finite/multisite.h"

namespace settings {
    inline constexpr bool debug_nbody_ham   = false;
    inline constexpr bool verbose_nbody_ham = false;
}

ModelFinite::ModelFinite() = default; // Can't initialize lists since we don't know the model size yet

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise, unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
ModelFinite::~ModelFinite()                               = default; // default dtor
ModelFinite:: ModelFinite(ModelFinite &&other)            = default; // default move ctor
ModelFinite  &ModelFinite::operator=(ModelFinite &&other) = default; // default move assign

/* clang-format off */
ModelFinite::ModelFinite(const ModelFinite &other) :
    cache(other.cache),
    active_sites(other.active_sites),
    model_type(other.model_type)
{
    MPO.clear();
    MPO.reserve(other.MPO.size());
    for(const auto &other_mpo : other.MPO) MPO.emplace_back(other_mpo->clone());
    if constexpr (settings::debug)
        for(const auto &[idx,other_mpo] : iter::enumerate(other.MPO))
            if(MPO[idx]->get_unique_id() != other_mpo->get_unique_id()) throw std::runtime_error("ID mismatch after copying mpo");
}
/* clang-format on */

ModelFinite &ModelFinite::operator=(const ModelFinite &other) {
    // check for self-assignment
    if(this != &other) {
        cache = other.cache;
        MPO.clear();
        MPO.reserve(other.MPO.size());
        for(const auto &other_mpo : other.MPO) MPO.emplace_back(other_mpo->clone());
        active_sites = other.active_sites;
        model_type   = other.model_type;
        if constexpr(settings::debug)
            for(const auto &[idx, other_mpo] : iter::enumerate(other.MPO))
                if(MPO[idx]->get_unique_id() != other_mpo->get_unique_id()) throw std::runtime_error("ID mismatch after copying mpo");
    }
    return *this;
}

void ModelFinite::initialize(ModelType model_type_, size_t model_size) {
    tools::log->info("Initializing model {} with {} sites", enum2sv(model_type_), model_size);
    if(model_size < 2) throw except::logic_error("Tried to initialize model with less than 2 sites");
    if(model_size > 2048) throw except::logic_error("Tried to initialize model with more than 2048 sites");
    if(not MPO.empty()) throw except::logic_error("Tried to initialize over an existing model. This is usually not what you want!");
    // Generate MPO
    model_type = model_type_;
    for(size_t site = 0; site < model_size; site++) { MPO.emplace_back(MpoFactory::create_mpo(site, model_type)); }
    if(MPO.size() != model_size) throw except::logic_error("Initialized MPO with wrong size");
}

const MpoSite &ModelFinite::get_mpo(size_t pos) const {
    if(pos >= MPO.size()) throw except::range_error("get_mpo(pos) pos out of range: {}", pos);
    return **std::next(MPO.begin(), safe_cast<long>(pos));
}

MpoSite &ModelFinite::get_mpo(size_t pos) { return const_cast<MpoSite &>(static_cast<const ModelFinite &>(*this).get_mpo(pos)); }

size_t ModelFinite::get_length() const { return MPO.size(); }

bool ModelFinite::is_real() const {
    for(const auto &mpo : MPO)
        if(not mpo->is_real()) return false;
    ;
    return true;
}

bool ModelFinite::has_nan() const {
    for(const auto &mpo : MPO)
        if(mpo->has_nan()) return true;
    return false;
}

void ModelFinite::assert_validity() const {
    for(const auto &mpo : MPO) mpo->assert_validity();
    if(settings::model::model_type == ModelType::ising_sdual) {
        for(const auto &mpo : MPO)
            if(not mpo->is_real()) throw except::runtime_error("model has imaginary part at mpo position {}", mpo->get_position());
    }
}

// For energy-shifted MPO's
bool ModelFinite::has_energy_shifted_mpo() const {
    bool shifted = MPO.front()->has_energy_shifted_mpo();
    for(const auto &mpo : MPO)
        if(shifted != mpo->has_energy_shifted_mpo())
            throw except::logic_error(
                fmt::format("Mismatching has_energy_shifted_mpo: pos 0:{} | {}:{}", shifted, mpo->get_position(), mpo->has_energy_shifted_mpo()));
    if(shifted and not settings::precision::use_energy_shifted_mpo)
        throw except::logic_error("The MPO's are energy-shifted but settings::precision::use_energy_shifted_mpo is false");
    return shifted;
}

bool ModelFinite::has_compressed_mpo_squared() const {
    bool compressed = MPO.front()->has_compressed_mpo_squared();
    for(const auto &mpo : MPO)
        if(compressed != mpo->has_compressed_mpo_squared())
            throw except::runtime_error("MPO² at pos 0: has_compressed_mpo_squared == {}, but MPO² at pos {} has_compressed_mpo_squared == {}", compressed,
                                        mpo->get_position(), mpo->has_compressed_mpo_squared());
    return compressed;
}

double ModelFinite::get_energy_shift_mpo() const {
    // Check that all energies are the same
    double e_shift = MPO.front()->get_energy_shift_mpo().real();
    for(const auto &mpo : MPO)
        if(mpo->get_energy_shift_mpo() != e_shift) throw std::runtime_error("Shifted energy per site mismatch!");
    return e_shift * static_cast<double>(get_length());
}

std::vector<std::any> ModelFinite::get_parameter(std::string_view fieldname) {
    std::vector<std::any> fields;
    for(const auto &mpo : MPO) { fields.emplace_back(mpo->get_parameter(fieldname)); }
    return fields;
}

double ModelFinite::get_energy_upper_bound() const { return MPO.front()->get_global_energy_upper_bound(); }

void ModelFinite::randomize() {
    tools::log->info("Randomizing hamiltonian");
    std::vector<MpoSite::TableMap> all_params;
    for(const auto &mpo : MPO) {
        mpo->randomize_hamiltonian();
        all_params.emplace_back(mpo->get_parameters());
    }
    for(const auto &mpo : MPO) mpo->set_averages(all_params, false);
}

void ModelFinite::build_mpo() {
    tools::log->debug("Building MPO");
    cache.multisite_ham   = std::nullopt;
    cache.multisite_mpo   = std::nullopt;
    cache.multisite_ham_t = std::nullopt;
    cache.multisite_mpo_t = std::nullopt;
    for(const auto &mpo : MPO) mpo->build_mpo();
    for(const auto &mpo : MPO) mpo->build_mpo_t();
}

void ModelFinite::build_mpo_squared() {
    tools::log->debug("Building MPO²");
    cache.multisite_mpo_squared = std::nullopt;
    cache.multisite_ham_squared = std::nullopt;
    for(const auto &mpo : MPO) mpo->build_mpo_squared();
}

void ModelFinite::clear_mpo_squared() {
    tools::log->debug("Clearing MPO²");
    cache.multisite_mpo_squared = std::nullopt;
    cache.multisite_ham_squared = std::nullopt;
    for(const auto &mpo : MPO) mpo->clear_mpo_squared();
}

void ModelFinite::compress_mpo() {
    cache.multisite_mpo     = std::nullopt;
    cache.multisite_mpo_ids = std::nullopt;
    cache.multisite_ham     = std::nullopt;
    auto mpo_compressed     = get_compressed_mpos();
    for(const auto &[pos, mpo] : iter::enumerate(MPO)) mpo->set_mpo(mpo_compressed[pos]);
}

void ModelFinite::compress_mpo_squared() {
    cache.multisite_mpo_squared = std::nullopt;
    cache.multisite_ham_squared = std::nullopt;
    auto mpo_squared_compressed = get_compressed_mpos_squared();
    for(const auto &[pos, mpo] : iter::enumerate(MPO)) mpo->set_mpo_squared(mpo_squared_compressed[pos]);
}

bool ModelFinite::has_mpo() const {
    return std::all_of(MPO.begin(), MPO.end(), [](const auto &mpo) { return mpo->has_mpo(); });
}

bool ModelFinite::has_mpo_squared() const {
    return std::all_of(MPO.begin(), MPO.end(), [](const auto &mpo) { return mpo->has_mpo_squared(); });
}

std::vector<std::reference_wrapper<const MpoSite>> ModelFinite::get_mpo(const std::vector<size_t> &sites) const {
    std::vector<std::reference_wrapper<const MpoSite>> mpos;
    mpos.reserve(sites.size());
    for(auto &site : sites) { mpos.emplace_back(get_mpo(site)); }
    return mpos;
}

std::vector<std::reference_wrapper<MpoSite>> ModelFinite::get_mpo(const std::vector<size_t> &sites) {
    std::vector<std::reference_wrapper<MpoSite>> mpos;
    mpos.reserve(sites.size());
    for(auto &site : sites) { mpos.emplace_back(get_mpo(site)); }
    return mpos;
}

std::vector<std::reference_wrapper<const MpoSite>> ModelFinite::get_mpo_active() const {
    std::vector<std::reference_wrapper<const MpoSite>> mpos;
    mpos.reserve(active_sites.size());
    for(auto &site : active_sites) { mpos.emplace_back(get_mpo(site)); }
    return mpos;
}

std::vector<std::reference_wrapper<MpoSite>> ModelFinite::get_mpo_active() {
    std::vector<std::reference_wrapper<MpoSite>> mpos;
    mpos.reserve(active_sites.size());
    for(auto &site : active_sites) { mpos.emplace_back(get_mpo(site)); }
    return mpos;
}

std::vector<Eigen::Tensor<cplx, 4>> ModelFinite::get_all_mpo_tensors(MposWithEdges withEdges) {
    tools::log->trace("Collecting all MPO: {} sites | with edges {}", MPO.size(), static_cast<std::underlying_type_t<MposWithEdges>>(withEdges));
    // Collect all the mpo (doesn't matter if they are already compressed)
    std::vector<Eigen::Tensor<cplx, 4>> mpos;
    mpos.reserve(MPO.size());
    for(const auto &mpo : MPO) mpos.emplace_back(mpo->MPO());
    switch(withEdges) {
        case MposWithEdges::OFF: return mpos;
        case MposWithEdges::ON: {
            auto ledge = MPO.front()->get_MPO_edge_left();
            auto redge = MPO.back()->get_MPO_edge_right();
            return tools::finite::mpo::get_mpos_with_edges(mpos, ledge, redge);
        }
        default: throw except::runtime_error("Unrecognized enum value <CompressWithEdges>: {}", static_cast<std::underlying_type_t<MposWithEdges>>(withEdges));
    }
}

std::vector<Eigen::Tensor<cplx_t, 4>> ModelFinite::get_all_mpo_tensors_t(MposWithEdges withEdges) {
    tools::log->trace("Collecting all MPO_t: {} sites | with edges {}", MPO.size(), static_cast<std::underlying_type_t<MposWithEdges>>(withEdges));
    // Collect all the mpo (doesn't matter if they are already compressed)
    std::vector<Eigen::Tensor<cplx_t, 4>> mpos_t;
    mpos_t.reserve(MPO.size());

    for(const auto &mpo : MPO) mpos_t.emplace_back(mpo->MPO_t());
    switch(withEdges) {
        case MposWithEdges::OFF: return mpos_t;
        case MposWithEdges::ON: {
            Eigen::Tensor<cplx_t, 1> ledge = MPO.front()->get_MPO_edge_right().cast<cplx_t>();
            Eigen::Tensor<cplx_t, 1> redge = MPO.back()->get_MPO_edge_right().cast<cplx_t>();
            return tools::finite::mpo::get_mpos_with_edges_t(mpos_t, ledge, redge);
        }
        default: throw except::runtime_error("Unrecognized enum value <CompressWithEdges>: {}", static_cast<std::underlying_type_t<MposWithEdges>>(withEdges));
    }
}

std::vector<Eigen::Tensor<cplx, 4>> ModelFinite::get_compressed_mpos(MposWithEdges withEdges) {
    auto mpoComp = settings::precision::use_compressed_mpo;
    tools::log->trace("Compressing MPO: {} sites | with edges {} | compression {}", MPO.size(), enum2sv(withEdges), enum2sv(mpoComp));
    // Collect all the mpo (doesn't matter if they are already compressed)
    std::vector<Eigen::Tensor<cplx, 4>> mpos;
    mpos.reserve(MPO.size());
    for(const auto &mpo : MPO) mpos.emplace_back(mpo->MPO());
    switch(withEdges) {
        case MposWithEdges::OFF: return tools::finite::mpo::get_compressed_mpos(mpos, mpoComp);
        case MposWithEdges::ON: {
            auto ledge = MPO.front()->get_MPO_edge_left();
            auto redge = MPO.back()->get_MPO_edge_right();
            return tools::finite::mpo::get_compressed_mpos(mpos, ledge, redge, mpoComp);
        }
        default: throw except::runtime_error("Unrecognized enum value <CompressWithEdges>: {}", static_cast<std::underlying_type_t<MposWithEdges>>(withEdges));
    }
}

std::vector<Eigen::Tensor<cplx, 4>> ModelFinite::get_compressed_mpos_squared(MposWithEdges withEdges) {
    auto mpoComp = settings::precision::use_compressed_mpo_squared;
    tools::log->trace("Compressing MPO²: {} sites | with edges {} | compression type {}", MPO.size(), enum2sv(withEdges), enum2sv(mpoComp));
    if(not has_mpo_squared()) build_mpo_squared(); // Make sure they exist.
    // Collect all the mpo² (doesn't matter if they are already compressed)
    std::vector<Eigen::Tensor<cplx, 4>> mpos_sq;
    mpos_sq.reserve(MPO.size());
    for(const auto &mpo : MPO) mpos_sq.emplace_back(mpo->MPO2());
    switch(withEdges) {
        case MposWithEdges::OFF: return tools::finite::mpo::get_compressed_mpos(mpos_sq, mpoComp); break;
        case MposWithEdges::ON: {
            auto ledge = MPO.front()->get_MPO2_edge_left();
            auto redge = MPO.back()->get_MPO2_edge_right();
            return tools::finite::mpo::get_compressed_mpos(mpos_sq, ledge, redge, mpoComp);
            break;
        }
        default: throw except::runtime_error("Unrecognized enum value <CompressWithEdges>: {}", static_cast<std::underlying_type_t<MposWithEdges>>(withEdges));
    }
}

void ModelFinite::set_energy_shift_mpo(double energy_shift) {
    if(std::abs(get_energy_shift_mpo() - energy_shift) <= std::numeric_limits<double>::epsilon()) { return; }
    tools::log->trace("Shifting MPO energy: {:.16f}", energy_shift);
    double energy_shift_per_site = energy_shift / static_cast<double>(get_length());
    for(const auto &mpo : MPO) mpo->set_energy_shift_mpo(cplx(energy_shift_per_site, 0.0));
    clear_cache();
}

void ModelFinite::set_parity_shift_mpo(OptRitz ritz, int sign, std::string_view axis) {
    if(not qm::spin::half::is_valid_axis(axis)) return;
    auto axus = qm::spin::half::get_axis_unsigned(axis);
    if(get_parity_shift_mpo() == std::make_tuple(ritz, sign, axus)) {
        tools::log->debug("set_parity_shift_mpo: not needed -- parity shift is already [{} {} {}]", enum2sv(ritz), sign, axus);
        return;
    }
    tools::log->info("Setting MPO parity shift for target axis {}{}", sign == 0 ? "" : (sign < 0 ? "-" : "+"), qm::spin::half::get_axis_unsigned(axis));
    for(const auto &mpo : MPO) mpo->set_parity_shift_mpo(ritz, sign, axis);
    clear_cache();
}

std::tuple<OptRitz, int, std::string_view> ModelFinite::get_parity_shift_mpo() const {
    auto parity_shift_front = MPO.front()->get_parity_shift_mpo();
    for(const auto &mpo : MPO) {
        if(mpo->get_parity_shift_mpo() != parity_shift_front)
            throw except::logic_error("mpo parity shift at site {} differs from shift at site 0", mpo->get_position());
    }
    return parity_shift_front;
}

bool ModelFinite::has_parity_shifted_mpo() const {
    auto shifted = MPO.front()->has_parity_shifted_mpo();
    for(const auto &mpo : MPO)
        if(shifted != mpo->has_parity_shifted_mpo())
            throw except::logic_error(
                fmt::format("Mismatching has_parity_shifted_mpo: pos 0:{} | {}:{}", shifted, mpo->get_position(), mpo->has_parity_shifted_mpo()));
    if(shifted and not settings::precision::use_parity_shifted_mpo)
        throw except::logic_error("The MPO's are parity-shifted but settings::precision::use_parity_shifted_mpo is false");
    return shifted;
}

bool ModelFinite::has_parity_shifted_mpo_squared() const {
    auto shifted = MPO.front()->has_parity_shifted_mpo2();
    for(const auto &mpo : MPO)
        if(shifted != mpo->has_parity_shifted_mpo2())
            throw except::logic_error(
                fmt::format("Mismatching has_parity_shifted_mpo2: pos 0:{} | {}:{}", shifted, mpo->get_position(), mpo->has_parity_shifted_mpo2()));
    if(shifted and not settings::precision::use_parity_shifted_mpo_squared)
        throw except::logic_error("The MPO's are parity-shifted but settings::precision::use_parity_shifted_mpo_squared is false");
    return shifted;
}

void ModelFinite::set_parity_shift_mpo_squared(int sign, std::string_view axis) {
    if(not qm::spin::half::is_valid_axis(axis)) return;
    if(get_parity_shift_mpo_squared() == std::make_pair(sign, axis)) return;
    tools::log->info("Setting MPO² parity shift for target axis {}{}", sign == 0 ? "" : (sign < 0 ? "-" : "+"), qm::spin::half::get_axis_unsigned(axis));
    for(const auto &mpo : MPO) mpo->set_parity_shift_mpo_squared(sign, axis);
    clear_cache();
}

std::pair<int, std::string_view> ModelFinite::get_parity_shift_mpo_squared() const {
    auto parity_shift     = std::pair<int, std::string_view>{0, ""};
    bool parity_shift_set = false;
    for(const auto &mpo : MPO) {
        if(not parity_shift_set) {
            parity_shift     = mpo->get_parity_shift_mpo_squared();
            parity_shift_set = true;
        } else if(parity_shift != mpo->get_parity_shift_mpo_squared())
            throw except::logic_error("mpo² parity shift at site {} differs from shift at site 0", mpo->get_position());
    }
    return parity_shift;
}

ModelLocal ModelFinite::get_local(const std::vector<size_t> &sites) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite hamiltonian tensor");
    auto mlocal    = ModelLocal();
    auto positions = num::range<size_t>(sites.front(), sites.back() + 1);
    auto skip      = std::optional<std::vector<size_t>>();
    if(sites != positions) {
        skip = std::vector<size_t>{};
        for(const auto &pos : positions) {
            if(std::find(sites.begin(), sites.end(), pos) == sites.end()) skip->emplace_back(pos);
        }
    }
    for(const auto &pos : positions) {
        bool do_skip = std::find(skip->begin(), skip->end(), pos) != skip->end();
        if(do_skip) continue;
        mlocal.mpos.emplace_back(get_mpo(pos).clone());
    }
    return mlocal;
}

ModelLocal ModelFinite::get_local() const { return get_local(active_sites); }

std::array<long, 4> ModelFinite::active_dimensions() const { return tools::finite::multisite::get_dimensions(*this); }

Eigen::Tensor<cplx, 4> ModelFinite::get_multisite_mpo(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody, bool with_edgeL,
                                                      bool with_edgeR) const {
    // Observe that nbody empty/nullopt have very different meanings
    //      - empty means that no interactions should be taken into account, effectively setting all J(i,j...) = 0
    //      - nullopt means that we want the default mpo with (everything on)
    //      - otherwise nbody with values like {1,2} would imply we want 1 and 2-body interactions turned on
    //      - if nbody has a 0 value in it, it means we want to make an attempt to account for double-counting in multisite mpos.

    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mpo tensor");
    if(sites == active_sites and cache.multisite_mpo and not nbody) return cache.multisite_mpo.value();

    auto                   nbody_str    = fmt::format("{}", nbody.has_value() ? nbody.value() : std::vector<size_t>{});
    auto                   t_mpo        = tid::tic_scope("get_multisite_mpo", tid::level::highest);
    constexpr auto         shuffle_idx  = tenx::array6{0, 3, 1, 4, 2, 5};
    constexpr auto         contract_idx = tenx::idx({1}, {0});
    auto                   positions    = num::range<size_t>(sites.front(), sites.back() + 1);
    auto                   skip         = std::vector<size_t>{};
    auto                   keep_log     = std::vector<size_t>();
    auto                   skip_log     = std::vector<size_t>();
    bool                   do_cache     = !with_edgeL and !with_edgeR and nbody.has_value() and nbody->back() > 1; // Caching doesn't make sense for nbody == 1
    auto                  &threads      = tenx::threads::get();
    Eigen::Tensor<cplx, 4> multisite_mpo, mpoL, mpoR;
    Eigen::Tensor<cplx, 2> mpoR_traced;
    // The hamiltonian is the lower left corner he full system mpo chain, which we can extract using edgeL and edgeR
    Eigen::Tensor<cplx, 1> edgeL = get_mpo(sites.front()).get_MPO_edge_left();
    Eigen::Tensor<cplx, 1> edgeR = get_mpo(sites.back()).get_MPO_edge_right();

    tools::log->trace("Contracting multisite mpo tensor with sites {} | nbody {} ", sites, nbody_str);

    if(sites != positions) {
        for(const auto &pos : positions) {
            if(std::find(sites.begin(), sites.end(), pos) == sites.end()) skip.emplace_back(pos);
        }
    }

    for(const auto &pos : positions) {
        if constexpr(settings::verbose_nbody_ham) tools::log->trace("contracting position {}", pos);
        // sites needs to be sorted, but may skip sites.
        // For instance, sites == {3,9} is valid. Then sites 4,5,6,7,8 are skipped.
        // When a site is skipped, we set the contribution from its interaction terms to zero and trace over it so that
        // the physical dimension doesn't grow.
        bool do_trace = std::find(skip.begin(), skip.end(), pos) != skip.end();
        if(pos == positions.front()) {
            auto t_pre = tid::tic_scope("prepending", tid::level::highest);
            if(nbody or not skip.empty()) {
                multisite_mpo = get_mpo(pos).MPO_nbody_view(nbody, skip);
            } else {
                multisite_mpo = get_mpo(pos).MPO();
            }
            if(do_cache) {
                if(do_trace) {
                    skip_log.emplace_back(pos);
                } else {
                    keep_log.emplace_back(pos);
                }
            }
            if(with_edgeL and pos == positions.front()) {
                /* We can prepend edgeL to the first mpo to reduce the size of subsequent operations.
                 * Start by converting the edge from a rank1 to a rank2 with a dummy index of size 1:
                 *                        2               2
                 *                        |               |
                 *    0---[L]---(1)(0)---[M]---1 =  0---[LM]---1
                 *                        |               |
                 *                        3               3
                 */
                mpoL          = edgeL.reshape(tenx::array2{1, edgeL.size()}).contract(multisite_mpo, tenx::idx({1}, {0}));
                multisite_mpo = mpoL;
            }
            if(with_edgeR and pos == positions.back()) {
                /* This only happens when positions.size() == 1
                 * We can append edgeR to the last mpo to reduce the size of subsequent operations.
                 * Start by converting the edge from a rank1 to a rank2 with a dummy index of size 1:
                 *         2                              1                       2
                 *         |                              |                       |
                 *    0---[M]---1  0---[R]---1   =  0---[MR]---3  [shuffle]  0---[MR]---1
                 *         |                              |                       |
                 *         3                              2                       3
                 */
                auto mpoR_edgeR = Eigen::Tensor<cplx, 4>(multisite_mpo.contract(edgeR.reshape(tenx::array2{edgeR.size(), 1}), tenx::idx({1}, {0})));
                multisite_mpo   = mpoR_edgeR.shuffle(tenx::array4{0, 3, 1, 2});
            }
            continue;
        }

        mpoL = multisite_mpo;
        mpoR = nbody or not skip.empty() ? get_mpo(pos).MPO_nbody_view(nbody, skip) : get_mpo(pos).MPO();

        if(with_edgeR and pos == positions.back()) {
            /* We can append edgeL to the first mpo to reduce the size of subsequent operations.
             * Start by converting the edge from a rank1 to a rank2 with a dummy index of size 1:
             *         2                              1                       2
             *         |                              |                       |
             *    0---[M]---1  0---[R]---1   =  0---[MR]---3  [shuffle]  0---[MR]---1
             *         |                              |                       |
             *         3                              2                       3
             */
            auto mpoR_edgeR = Eigen::Tensor<cplx, 4>(mpoR.contract(edgeR.reshape(std::array<long, 2>{edgeR.size(), 1}), tenx::idx({1}, {0})));
            mpoR            = mpoR_edgeR.shuffle(tenx::array4{0, 3, 1, 2});
        }
        tools::log->trace("contracting position {} | mpoL {} | mpoR {}", pos, mpoL.dimensions(), mpoR.dimensions());

        // Determine if this position adds to the physical dimension or if it will get traced over
        long dim0     = mpoL.dimension(0);
        long dim1     = mpoR.dimension(1);
        long dim2     = mpoL.dimension(2) * (do_trace ? 1l : mpoR.dimension(2));
        long dim3     = mpoL.dimension(3) * (do_trace ? 1l : mpoR.dimension(3));
        auto new_dims = std::array<long, 4>{dim0, dim1, dim2, dim3};
        multisite_mpo.resize(new_dims);
        // Generate a unique cache string for the mpo that will be generated.
        // If there is a match for the string in cache, use the corresponding mpo, otherwise we make it.
        if(do_cache) {
            if(do_trace) {
                skip_log.emplace_back(pos);
            } else {
                keep_log.emplace_back(pos);
            }
        }
        auto new_cache_string = fmt::format("keep{}|skip{}|nbody{}|dims{}", keep_log, skip_log, nbody_str, new_dims);
        if(do_cache and cache.multisite_mpo_temps.find(new_cache_string) != cache.multisite_mpo_temps.end()) {
            if constexpr(settings::verbose_nbody_ham) tools::log->trace("cache hit: {}", new_cache_string);
            multisite_mpo = cache.multisite_mpo_temps.at(new_cache_string);
        } else {
            if constexpr(settings::verbose_nbody_ham) tools::log->trace("cache new: {}", new_cache_string);
            if(do_trace) {
                auto t_skip = tid::tic_scope("skipping", tid::level::highest);
                // Trace the physical indices of this skipped mpo (this should trace an identity)
                mpoR_traced = mpoR.trace(tenx::array2{2, 3});
                mpoR_traced *= mpoR_traced.constant(cplx(0.5, 0.0)); // divide by 2 (after tracing identity)
                // Append it to the multisite mpo
                multisite_mpo.device(*threads->dev) = mpoL.contract(mpoR_traced, tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2}).reshape(new_dims);
            } else {
                auto t_app                          = tid::tic_scope("appending", tid::level::highest);
                multisite_mpo.device(*threads->dev) = mpoL.contract(mpoR, contract_idx).shuffle(shuffle_idx).reshape(new_dims);
            }
            // This intermediate multisite_mpo_t could be the result we are looking for at a later time, so cache it!
            if(do_cache) cache.multisite_mpo_temps[new_cache_string] = multisite_mpo;
        }
    }
    if(with_edgeL) assert(multisite_mpo.dimension(0) == 1);
    if(with_edgeR) assert(multisite_mpo.dimension(1) == 1);
    return multisite_mpo;
}

Eigen::Tensor<cplx_t, 4> ModelFinite::get_multisite_mpo_t(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody, bool with_edgeL,
                                                          bool with_edgeR) const {
    // Observe that nbody empty/nullopt have very different meanings
    //      - empty means that no interactions should be taken into account, effectively setting all J(i,j...) = 0
    //      - nullopt means that we want the default mpo with (everything on)
    //      - otherwise nbody with values like {1,2} would imply we want 1 and 2-body interactions turned on
    //      - if nbody has a 0 value in it, it means we want to make an attempt to account for double-counting in multisite mpos.

    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mpo tensor");
    if(sites == active_sites and cache.multisite_mpo_t and not nbody) return cache.multisite_mpo_t.value();

    auto                     nbody_str    = fmt::format("{}", nbody.has_value() ? nbody.value() : std::vector<size_t>{});
    auto                     t_mpo        = tid::tic_scope("get_multisite_mpo_t", tid::level::highest);
    constexpr auto           shuffle_idx  = tenx::array6{0, 3, 1, 4, 2, 5};
    constexpr auto           contract_idx = tenx::idx({1}, {0});
    auto                     positions    = num::range<size_t>(sites.front(), sites.back() + 1);
    auto                     skip         = std::vector<size_t>{};
    auto                     keep_log     = std::vector<size_t>();
    auto                     skip_log     = std::vector<size_t>();
    bool                     do_cache = !with_edgeL and !with_edgeR and nbody.has_value() and nbody->back() > 1; // Caching doesn't make sense for nbody == 1
    auto                    &threads  = tenx::threads::get();
    Eigen::Tensor<cplx_t, 4> multisite_mpo_t, mpoL, mpoR;
    Eigen::Tensor<cplx_t, 2> mpoR_traced;
    // The hamiltonian is the lower left corner he full system mpo chain, which we can extract using edgeL and edgeR
    Eigen::Tensor<cplx_t, 1> edgeL = get_mpo(sites.front()).get_MPO_edge_left<cplx_t>();
    Eigen::Tensor<cplx_t, 1> edgeR = get_mpo(sites.back()).get_MPO_edge_right<cplx_t>();

    tools::log->trace("Contracting multisite mpo tensor with sites {} | nbody {} ", sites, nbody_str);

    if(sites != positions) {
        for(const auto &pos : positions) {
            if(std::find(sites.begin(), sites.end(), pos) == sites.end()) skip.emplace_back(pos);
        }
    }

    for(const auto &pos : positions) {
        if constexpr(settings::verbose_nbody_ham) tools::log->trace("contracting position {}", pos);
        // sites needs to be sorted, but may skip sites.
        // For instance, sites == {3,9} is valid. Then sites 4,5,6,7,8 are skipped.
        // When a site is skipped, we set the contribution from its interaction terms to zero and trace over it so that
        // the physical dimension doesn't grow.
        bool do_trace = std::find(skip.begin(), skip.end(), pos) != skip.end();
        if(pos == positions.front()) {
            auto t_pre = tid::tic_scope("prepending", tid::level::highest);
            if(nbody or not skip.empty()) {
                multisite_mpo_t = get_mpo(pos).MPO_nbody_view_t(nbody, skip);
            } else {
                multisite_mpo_t = get_mpo(pos).MPO_t();
            }
            if(do_cache) {
                if(do_trace) {
                    skip_log.emplace_back(pos);
                } else {
                    keep_log.emplace_back(pos);
                }
            }
            if(with_edgeL and pos == positions.front()) {
                /* We can prepend edgeL to the first mpo to reduce the size of subsequent operations.
                 * Start by converting the edge from a rank1 to a rank2 with a dummy index of size 1:
                 *                        2               2
                 *                        |               |
                 *    0---[L]---(1)(0)---[M]---1 =  0---[LM]---1
                 *                        |               |
                 *                        3               3
                 */
                mpoL            = edgeL.reshape(tenx::array2{1, edgeL.size()}).contract(multisite_mpo_t, tenx::idx({1}, {0}));
                multisite_mpo_t = mpoL;
            }
            if(with_edgeR and pos == positions.back()) {
                /* This only happens when positions.size() == 1
                 * We can append edgeR to the last mpo to reduce the size of subsequent operations.
                 * Start by converting the edge from a rank1 to a rank2 with a dummy index of size 1:
                 *        2                              1                       2
                 *        |                              |                       |
                 *   0---[M]---1  0---[R]---1   =  0---[MR]---3  [shuffle]  0---[MR]---1
                 *        |                              |                       |
                 *        3                              2                       3
                 */
                auto mpoR_edgeR = Eigen::Tensor<cplx_t, 4>(multisite_mpo_t.contract(edgeR.reshape(tenx::array2{edgeR.size(), 1}), tenx::idx({1}, {0})));
                multisite_mpo_t = mpoR_edgeR.shuffle(tenx::array4{0, 3, 1, 2});
            }
            continue;
        }

        mpoL = multisite_mpo_t;
        mpoR = nbody or not skip.empty() ? get_mpo(pos).MPO_nbody_view_t(nbody, skip) : get_mpo(pos).MPO_t();

        if(with_edgeR and pos == positions.back()) {
            /* We can append edgeL to the first mpo to reduce the size of subsequent operations.
             * Start by converting the edge from a rank1 to a rank2 with a dummy index of size 1:
             *        2                              1                       2
             *        |                              |                       |
             *   0---[M]---1  0---[R]---1   =  0---[MR]---3  [shuffle]  0---[MR]---1
             *        |                              |                       |
             *        3                              2                       3
             */
            auto mpoR_edgeR = Eigen::Tensor<cplx_t, 4>(mpoR.contract(edgeR.reshape(tenx::array2{edgeR.size(), 1}), tenx::idx({1}, {0})));
            mpoR            = mpoR_edgeR.shuffle(tenx::array4{0, 3, 1, 2});
        }

        // Determine if this position adds to the physical dimension or if it will get traced over
        long dim0     = mpoL.dimension(0);
        long dim1     = mpoR.dimension(1);
        long dim2     = mpoL.dimension(2) * (do_trace ? 1l : mpoR.dimension(2));
        long dim3     = mpoL.dimension(3) * (do_trace ? 1l : mpoR.dimension(3));
        auto new_dims = std::array<long, 4>{dim0, dim1, dim2, dim3};
        multisite_mpo_t.resize(new_dims);
        // Generate a unique cache string for the mpo that will be generated.
        // If there is a match for the string in cache, use the corresponding mpo, otherwise we make it.
        if(do_cache) {
            if(do_trace) {
                skip_log.emplace_back(pos);
            } else {
                keep_log.emplace_back(pos);
            }
        }
        auto new_cache_string = fmt::format("keep{}|skip{}|nbody{}|dims{}", keep_log, skip_log, nbody_str, new_dims);
        if(do_cache and cache.multisite_mpo_t_temps.find(new_cache_string) != cache.multisite_mpo_t_temps.end()) {
            if constexpr(settings::verbose_nbody_ham) tools::log->trace("cache hit: {}", new_cache_string);
            multisite_mpo_t = cache.multisite_mpo_t_temps.at(new_cache_string);
        } else {
            if constexpr(settings::verbose_nbody_ham) tools::log->trace("cache new: {}", new_cache_string);
            if(do_trace) {
                auto t_skip = tid::tic_scope("skipping", tid::level::highest);
                // Trace the physical indices of this skipped mpo (this should trace an identity)
                mpoR_traced = mpoR.trace(tenx::array2{2, 3});
                mpoR_traced *= mpoR_traced.constant(static_cast<real_t>(0.5)); // divide by 2 (after tracing identity)
                // Append it to the multisite mpo
                multisite_mpo_t.device(*threads->dev) = mpoL.contract(mpoR_traced, tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2}).reshape(new_dims);
            } else {
                auto t_app                            = tid::tic_scope("appending", tid::level::highest);
                multisite_mpo_t.device(*threads->dev) = mpoL.contract(mpoR, contract_idx).shuffle(shuffle_idx).reshape(new_dims);
            }
            // This intermediate multisite_mpo_t could be the result we are looking for at a later time, so cache it!
            if(do_cache) cache.multisite_mpo_t_temps[new_cache_string] = multisite_mpo_t;
        }
    }
    if(with_edgeL) assert(multisite_mpo_t.dimension(0) == 1);
    if(with_edgeR) assert(multisite_mpo_t.dimension(1) == 1);
    return multisite_mpo_t;
}

Eigen::Tensor<cplx, 2> ModelFinite::get_multisite_ham(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody) const {
    /*! We use edges without particle content to get a hamiltonian for a subsystem.
     *  To get an effective hamiltonian with particle-full edges, use TensorsFinite::get_effective_hamiltonian.
     */
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite hamiltonian tensor");
    if(sites == active_sites and cache.multisite_ham and not nbody) {
        tools::log->info("cache hit: sites{}|nbody{}", sites, nbody.has_value() ? nbody.value() : std::vector<size_t>{});
        return cache.multisite_ham.value();
    }
    tools::log->trace("get_multisite_ham(): Contracting effective Hamiltonian");
    long spin_dim = 1;
    for(const auto &pos : sites) { spin_dim *= get_mpo(pos).get_spin_dimension(); }
    auto dim2 = tenx::array2{spin_dim, spin_dim};
    if(sites.size() <= 4) {
        return get_multisite_mpo(sites, nbody, true, true).reshape(dim2);
    } else {
        // When there are many sites, it is better to split sites into two equal chunks and then merge them (because edgeL/edgeR makes them small)
        auto half   = static_cast<long>((sites.size() + 1) / 2); // Let the left side take one more site in odd cases, because we contract from the left
        auto sitesL = std::vector<size_t>(sites.begin(), sites.begin() + half);
        auto sitesR = std::vector<size_t>(sites.begin() + half, sites.end());
        auto mpoL   = get_multisite_mpo(sitesL, nbody, true, false); // Shuffle so we can use GEMM
        auto mpoR   = get_multisite_mpo(sitesR, nbody, false, true);
        auto mpoLR  = tenx::gemm_mpo(mpoL, mpoR);
        return mpoLR.reshape(tenx::array2{spin_dim, spin_dim});
    }
}

Eigen::Tensor<cplx_t, 2> ModelFinite::get_multisite_ham_t(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite hamiltonian tensor");
    if(sites == active_sites and cache.multisite_ham_t and not nbody) {
        tools::log->info("cache hit: sites{}|nbody{}", sites, nbody.has_value() ? nbody.value() : std::vector<size_t>{});
        return cache.multisite_ham_t.value();
    }
    long spin_dim = 1;
    for(const auto &pos : sites) { spin_dim *= get_mpo(pos).get_spin_dimension(); }
    auto dim2 = tenx::array2{spin_dim, spin_dim};
    if(sites.size() < 4) {
        return get_multisite_mpo_t(sites, nbody, true, true).reshape(dim2);
    } else {
        // When there are many sites, it's beneficial to split sites into two equal chunks and then merge them (because edgeL/edgeR makes them small)
        auto half   = static_cast<long>((sites.size() + 1) / 2); // Let the left side take one more site in odd cases, because we contract from the left
        auto sitesL = std::vector<size_t>(sites.begin(), sites.begin() + half);
        auto sitesR = std::vector<size_t>(sites.begin() + half, sites.end());
        auto mpoL   = get_multisite_mpo_t(sitesL, nbody, true, false); // Shuffle so we can use GEMM
        auto mpoR   = get_multisite_mpo_t(sitesR, nbody, false, true);
        auto mpoLR  = tenx::gemm_mpo(mpoL, mpoR);
        return mpoLR.reshape(tenx::array2{spin_dim, spin_dim});
    }
}

const Eigen::Tensor<cplx, 4> &ModelFinite::get_multisite_mpo() const {
    if(cache.multisite_mpo and cache.multisite_mpo_ids and not active_sites.empty()) {
        if constexpr(settings::debug) tools::log->trace("get_multisite_mpo: cache hit");
        // Check that the ids match
        auto active_ids = get_active_ids();
        if(active_ids != cache.multisite_mpo_ids)
            throw except::runtime_error("get_multisite_mpo: cache has mismatching ids: active ids: {} | cache ids {}", active_ids,
                                        cache.multisite_mpo_ids.value());
        return cache.multisite_mpo.value();
    }
    cache.multisite_mpo     = get_multisite_mpo(active_sites);
    cache.multisite_mpo_ids = get_active_ids();
    return cache.multisite_mpo.value();
}

const Eigen::Tensor<cplx_t, 4> &ModelFinite::get_multisite_mpo_t() const {
    if(cache.multisite_mpo_t and not active_sites.empty()) {
        if constexpr(settings::debug) tools::log->trace("multisite_mpo_t: cache hit");
        return cache.multisite_mpo_t.value();
    }
    cache.multisite_mpo_t = get_multisite_mpo_t(active_sites);
    return cache.multisite_mpo_t.value();
}

const Eigen::Tensor<cplx, 2> &ModelFinite::get_multisite_ham() const {
    if(cache.multisite_ham and not active_sites.empty()) return cache.multisite_ham.value();
    cache.multisite_ham = get_multisite_ham(active_sites);
    return cache.multisite_ham.value();
}

const Eigen::Tensor<cplx_t, 2> &ModelFinite::get_multisite_ham_t() const {
    if(cache.multisite_ham_t and not active_sites.empty()) return cache.multisite_ham_t.value();
    cache.multisite_ham_t = get_multisite_ham_t(active_sites);
    return cache.multisite_ham_t.value();
}

Eigen::Tensor<cplx, 4> ModelFinite::get_multisite_mpo_shifted_view(double energy_per_site) const {
    auto                   t_mpo = tid::tic_scope("mpo_shifted_view");
    Eigen::Tensor<cplx, 4> multisite_mpo, temp;
    constexpr auto         shuffle_idx  = tenx::array6{0, 3, 1, 4, 2, 5};
    constexpr auto         contract_idx = tenx::idx({1}, {0});
    auto                  &threads      = tenx::threads::get();
    for(const auto &site : active_sites) {
        if(multisite_mpo.size() == 0) {
            multisite_mpo = get_mpo(site).MPO_energy_shifted_view(energy_per_site);
            continue;
        }
        const auto         &mpo      = get_mpo(site);
        long                dim0     = multisite_mpo.dimension(0);
        long                dim1     = mpo.MPO().dimension(1);
        long                dim2     = multisite_mpo.dimension(2) * mpo.MPO().dimension(2);
        long                dim3     = multisite_mpo.dimension(3) * mpo.MPO().dimension(3);
        std::array<long, 4> new_dims = {dim0, dim1, dim2, dim3};
        temp.resize(new_dims);
        temp.device(*threads->dev) = multisite_mpo.contract(mpo.MPO_energy_shifted_view(energy_per_site), contract_idx).shuffle(shuffle_idx).reshape(new_dims);
        multisite_mpo              = std::move(temp);
    }
    return multisite_mpo;
}

Eigen::Tensor<cplx, 4> ModelFinite::get_multisite_mpo_squared_shifted_view(double energy_per_site) const {
    auto                   multisite_mpo_shifted = get_multisite_mpo_shifted_view(energy_per_site);
    long                   dim0                  = multisite_mpo_shifted.dimension(0) * multisite_mpo_shifted.dimension(0);
    long                   dim1                  = multisite_mpo_shifted.dimension(1) * multisite_mpo_shifted.dimension(1);
    long                   dim2                  = multisite_mpo_shifted.dimension(2);
    long                   dim3                  = multisite_mpo_shifted.dimension(3);
    std::array<long, 4>    mpo_squared_dims      = {dim0, dim1, dim2, dim3};
    Eigen::Tensor<cplx, 4> multisite_mpo_squared_shifted(mpo_squared_dims);
    auto                  &threads = tenx::threads::get();
    multisite_mpo_squared_shifted.device(*threads->dev) =
        multisite_mpo_shifted.contract(multisite_mpo_shifted, tenx::idx({3}, {2})).shuffle(tenx::array6{0, 3, 1, 4, 2, 5}).reshape(mpo_squared_dims);
    return multisite_mpo_squared_shifted;
}

std::array<long, 4> ModelFinite::active_dimensions_squared() const { return tools::finite::multisite::get_dimensions_squared(*this); }

// Eigen::Tensor<cplx, 4> ModelFinite::get_multisite_mpo_squared(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody, bool with_edgeL,
//                                                               bool with_edgeR) const {
//     // Observe that nbody empty/nullopt have very different meanings
//     //      - empty means that no interactions should be taken into account, effectively setting all J(i,j...) = 0
//     //      - nullopt means that we want the default mpo with (everything on)
//     //      - otherwise nbody with values like {1,2} would imply we want 1 and 2-body interactions turned on
//     //      - if nbody has a 0 value in it, it means we want to make an attempt to account for double-counting in multisite mpos.
//
//     if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mpo tensor");
//     if(sites == active_sites and cache.multisite_mpo_squared and not nbody) return cache.multisite_mpo_squared.value();
//
//     auto                   nbody_str    = fmt::format("{}", nbody.has_value() ? nbody.value() : std::vector<size_t>{});
//     auto                   t_mpo2       = tid::tic_scope("get_multisite_mpo_squared", tid::level::highest);
//     constexpr auto         shuffle_idx  = tenx::array6{0, 3, 1, 4, 2, 5};
//     constexpr auto         contract_idx = tenx::idx({1}, {0});
//     auto                   positions    = num::range<size_t>(sites.front(), sites.back() + 1);
//     auto                   skip         = std::vector<size_t>{};
//     auto                   keep_log     = std::vector<size_t>();
//     auto                   skip_log     = std::vector<size_t>();
//     bool                   do_cache     = !with_edgeL and !with_edgeR and nbody.has_value() and nbody->back() > 1; // Caching doesn't make sense for nbody == 1
//     auto                  &threads      = tenx::threads::get();
//     Eigen::Tensor<cplx, 4> multisite_mpo_squared, mpo2L, mpo2R;
//     Eigen::Tensor<cplx, 2> mpo2R_traced;
//     // The hamiltonian is the lower left corner he full system mpo chain, which we can extract using edgeL and edgeR
//     Eigen::Tensor<cplx, 1> edgeL = get_mpo(sites.front()).get_MPO2_edge_left();
//     Eigen::Tensor<cplx, 1> edgeR = get_mpo(sites.back()).get_MPO2_edge_right();
//
//     tools::log->trace("Contracting multisite mpo squared tensor with sites {} | nbody {} ", sites, nbody_str);
//
//     if(sites != positions) {
//         for(const auto &pos : positions) {
//             if(std::find(sites.begin(), sites.end(), pos) == sites.end()) skip.emplace_back(pos);
//         }
//     }
//
//     for(const auto &pos : positions) {
//         if constexpr(settings::verbose_nbody_ham) tools::log->trace("contracting position {}", pos);
//         // sites needs to be sorted, but may skip sites.
//         // For instance, sites == {3,9} is valid. Then sites 4,5,6,7,8 are skipped.
//         // When a site is skipped, we set the contribution from its interaction terms to zero and trace over it so that
//         // the physical dimension doesn't grow.
//         bool do_trace = std::find(skip.begin(), skip.end(), pos) != skip.end();
//         if(pos == positions.front()) {
//             auto t_pre = tid::tic_scope("prepending", tid::level::highest);
//             if(nbody or not skip.empty()) {
//                 multisite_mpo_squared = get_mpo(pos).MPO2_nbody_view(nbody, skip);
//             } else {
//                 multisite_mpo_squared = get_mpo(pos).MPO2();
//             }
//             if(do_cache) {
//                 if(do_trace) {
//                     skip_log.emplace_back(pos);
//                 } else {
//                     keep_log.emplace_back(pos);
//                 }
//             }
//             if(with_edgeL and pos == positions.front()) {
//                 /* We can prepend edgeL to the first mpo to reduce the size of subsequent operations.
//                  * Start by converting the edge from a rank1 to a rank2 with a dummy index of size 1:
//                  *                        2               2
//                  *                        |               |
//                  *    0---[L]---(1)(0)---[M]---1 =  0---[LM]---1
//                  *                        |               |
//                  *                        3               3
//                  */
//                 mpo2L                 = edgeL.reshape(tenx::array2{1, edgeL.size()}).contract(multisite_mpo_squared, tenx::idx({1}, {0}));
//                 multisite_mpo_squared = mpo2L;
//             }
//             if(with_edgeR and pos == positions.back()) {
//                 /* This only happens when positions.size() == 1
//                  * We can append edgeR to the last mpo to reduce the size of subsequent operations.
//                  * Start by converting the edge from a rank1 to a rank2 with a dummy index of size 1:
//                  *         2                              1                       2
//                  *         |                              |                       |
//                  *    0---[M]---1  0---[R]---1   =  0---[MR]---3  [shuffle]  0---[MR]---1
//                  *         |                              |                       |
//                  *         3                              2                       3
//                  */
//                 auto mpo2R_edgeR = Eigen::Tensor<cplx, 4>(multisite_mpo_squared.contract(edgeR.reshape(tenx::array2{edgeR.size(), 1}), tenx::idx({1}, {0})));
//                 multisite_mpo_squared = mpo2R_edgeR.shuffle(tenx::array4{0, 3, 1, 2});
//             }
//             continue;
//         }
//
//         mpo2L = multisite_mpo_squared;
//         mpo2R = nbody or not skip.empty() ? get_mpo(pos).MPO2_nbody_view(nbody, skip) : get_mpo(pos).MPO2();
//
//         if(with_edgeR and pos == positions.back()) {
//             /* We can append edgeL to the first mpo to reduce the size of subsequent operations.
//              * Start by converting the edge from a rank1 to a rank2 with a dummy index of size 1:
//              *         2                              1                       2
//              *         |                              |                       |
//              *    0---[M]---1  0---[R]---1   =  0---[MR]---3  [shuffle]  0---[MR]---1
//              *         |                              |                       |
//              *         3                              2                       3
//              */
//             auto mpoR_edgeR = Eigen::Tensor<cplx, 4>(mpo2R.contract(edgeR.reshape(std::array<long, 2>{edgeR.size(), 1}), tenx::idx({1}, {0})));
//             mpo2R           = mpoR_edgeR.shuffle(tenx::array4{0, 3, 1, 2});
//         }
//         tools::log->trace("contracting position {} | mpoL {} | mpoR {}", pos, mpo2L.dimensions(), mpo2R.dimensions());
//
//         // Determine if this position adds to the physical dimension or if it will get traced over
//         long dim0     = mpo2L.dimension(0);
//         long dim1     = mpo2R.dimension(1);
//         long dim2     = mpo2L.dimension(2) * (do_trace ? 1l : mpo2R.dimension(2));
//         long dim3     = mpo2L.dimension(3) * (do_trace ? 1l : mpo2R.dimension(3));
//         auto new_dims = std::array<long, 4>{dim0, dim1, dim2, dim3};
//         multisite_mpo_squared.resize(new_dims);
//         // Generate a unique cache string for the mpo that will be generated.
//         // If there is a match for the string in cache, use the corresponding mpo, otherwise we make it.
//         if(do_cache) {
//             if(do_trace) {
//                 skip_log.emplace_back(pos);
//             } else {
//                 keep_log.emplace_back(pos);
//             }
//         }
//         auto new_cache_string = fmt::format("keep{}|skip{}|nbody{}|dims{}", keep_log, skip_log, nbody_str, new_dims);
//         if(do_cache and cache.multisite_mpo_squared_temps.find(new_cache_string) != cache.multisite_mpo_squared_temps.end()) {
//             if constexpr(settings::verbose_nbody_ham) tools::log->trace("cache hit: {}", new_cache_string);
//             multisite_mpo_squared = cache.multisite_mpo_squared_temps.at(new_cache_string);
//         } else {
//             if constexpr(settings::verbose_nbody_ham) tools::log->trace("cache new: {}", new_cache_string);
//             if(do_trace) {
//                 auto t_skip = tid::tic_scope("skipping", tid::level::highest);
//                 // Trace the physical indices of this skipped mpo (this should trace an identity)
//                 mpo2R_traced = mpo2R.trace(tenx::array2{2, 3});
//                 mpo2R_traced *= mpo2R_traced.constant(cplx(0.5, 0.0)); // divide by 2 (after tracing identity)
//                 // Append it to the multisite mpo
//                 multisite_mpo_squared.device(*threads->dev) =
//                     mpo2L.contract(mpo2R_traced, tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2}).reshape(new_dims);
//             } else {
//                 auto t_app                                  = tid::tic_scope("appending", tid::level::highest);
//                 multisite_mpo_squared.device(*threads->dev) = mpo2L.contract(mpo2R, contract_idx).shuffle(shuffle_idx).reshape(new_dims);
//             }
//             // This intermediate multisite_mpo_t could be the result we are looking for at a later time, so cache it!
//             if(do_cache) cache.multisite_mpo_squared_temps[new_cache_string] = multisite_mpo_squared;
//         }
//     }
//     if(with_edgeL) assert(multisite_mpo_squared.dimension(0) == 1);
//     if(with_edgeR) assert(multisite_mpo_squared.dimension(1) == 1);
//     return multisite_mpo_squared;
// }

Eigen::Tensor<cplx, 4> ModelFinite::get_multisite_mpo_squared(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mpo squared tensor");
    if(sites == active_sites and cache.multisite_mpo_squared and not nbody) return cache.multisite_mpo_squared.value();
    tools::log->trace("Contracting multisite mpo² tensor with sites {}", sites);
    auto                   t_mpo     = tid::tic_scope("get_multisite_mpo_squared", tid::level::highest);
    auto                   positions = num::range<size_t>(sites.front(), sites.back() + 1);
    Eigen::Tensor<cplx, 4> multisite_mpo_squared;
    constexpr auto         shuffle_idx  = tenx::array6{0, 3, 1, 4, 2, 5};
    constexpr auto         contract_idx = tenx::idx({1}, {0});
    tenx::array4           new_dims;
    Eigen::Tensor<cplx, 4> temp;
    bool                   first   = true;
    auto                  &threads = tenx::threads::get();
    for(const auto &pos : positions) {
        // sites needs to be sorted, but may skip sites.
        // For instance, sites == {3,9} is valid. Then sites 4,5,6,7,8 are skipped.
        // When a site is skipped, we set the contribution from its interaction terms to zero and trace over it so that
        // the physical dimension doesn't grow.

        auto nbody_local = nbody;
        bool skip        = std::find(sites.begin(), sites.end(), pos) == sites.end();
        if(skip) nbody_local = std::vector<size_t>{};

        if(first) {
            if(nbody_local)
                multisite_mpo_squared = get_mpo(pos).MPO2_nbody_view(nbody_local);
            else
                multisite_mpo_squared = get_mpo(pos).MPO2();
            first = false;
            continue;
        }

        const auto &mpo  = get_mpo(pos);
        long        dim0 = multisite_mpo_squared.dimension(0);
        long        dim1 = mpo.MPO2().dimension(1);
        long        dim2 = multisite_mpo_squared.dimension(2) * mpo.MPO2().dimension(2);
        long        dim3 = multisite_mpo_squared.dimension(3) * mpo.MPO2().dimension(3);
        new_dims         = {dim0, dim1, dim2, dim3};
        temp.resize(new_dims);
        if(nbody_local) // Avoids creating a temporary
            temp.device(*threads->dev) = multisite_mpo_squared.contract(mpo.MPO2_nbody_view(nbody_local), contract_idx).shuffle(shuffle_idx).reshape(new_dims);
        else
            temp.device(*threads->dev) = multisite_mpo_squared.contract(mpo.MPO2(), contract_idx).shuffle(shuffle_idx).reshape(new_dims);

        if(skip) {
            /*! We just got handed a multisite-mpo created as
             *
             *       2   3                      2
             *       |   |                      |
             *  0 --[ mpo ]-- 1  --->    0 --[ mpo ]-- 1
             *       |   |                      |
             *       4   5                      3
             *
             *
             * In this step, a reshape brings back the 6 indices, and index 3 and 5 should be traced over.
             *
             */
            long                   d0    = dim0;
            long                   d1    = dim1;
            long                   d2    = multisite_mpo_squared.dimension(2);
            long                   d3    = mpo.MPO().dimension(2);
            long                   d4    = multisite_mpo_squared.dimension(3);
            long                   d5    = mpo.MPO().dimension(3);
            Eigen::Tensor<cplx, 4> temp2 = temp.reshape(tenx::array6{d0, d1, d2, d3, d4, d5}).trace(tenx::array2{3, 5});
            multisite_mpo_squared        = temp2 * temp2.constant(cplx(0.5, 0.0));
        } else {
            multisite_mpo_squared = std::move(temp);
        }
    }
    return multisite_mpo_squared;
}

Eigen::Tensor<cplx, 2> ModelFinite::get_multisite_ham_squared(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody) const {
    /*! We use edges without particle content to get a hamiltonian for a subsystem.
     *  To get an effective hamiltonian with particle-full edges, use TensorsFinite::get_effective_hamiltonian.
     */
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite hamiltonian squared tensor");
    if(sites == active_sites and cache.multisite_ham_squared and not nbody) {
        tools::log->info("cache hit: sites{}|nbody{}", sites, nbody.has_value() ? nbody.value() : std::vector<size_t>{});
        return cache.multisite_ham_squared.value();
    }
    tools::log->trace("get_multisite_ham_squared(): Contracting effective Hamiltonian squared");
    long spin_dim = 1;
    for(const auto &pos : sites) { spin_dim *= get_mpo(pos).get_spin_dimension(); }
    auto dim2 = tenx::array2{spin_dim, spin_dim};
    if(sites.size() <= 4) {
        return get_multisite_mpo_squared(sites, nbody).reshape(dim2);
    } else {
        // When there are many sites, it is better to split sites into two equal chunks and then merge them, especially if when taking all mpos.
        auto half   = static_cast<long>((sites.size() + 1) / 2); // Let the left side take one more site in odd cases, because we contract from the left
        auto sitesL = std::vector<size_t>(sites.begin(), sites.begin() + half);
        auto sitesR = std::vector<size_t>(sites.begin() + half, sites.end());
        auto mpo2L   = get_multisite_mpo_squared(sitesL, nbody);
        auto mpo2R   = get_multisite_mpo_squared(sitesR, nbody);
        auto mpo2LR  = tenx::gemm_mpo(mpo2L, mpo2R);
        return mpo2LR.reshape(tenx::array2{spin_dim, spin_dim});
    }
}

// Eigen::Tensor<cplx, 2> ModelFinite::get_multisite_ham_squared(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody) const {
/*! We use edges without particle content to get a hamiltonian for a subsystem.
 *  To get an effective hamiltonian with particle-full edges, use TensorsFinite::get_effective_hamiltonian_squared.
 */
// if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite hamiltonian tensor");
// if(sites == active_sites and cache.multisite_ham and not nbody) return cache.multisite_ham.value();
// tools::log->trace("ModelFinite::get_multisite_ham_squared(): contracting active sites {}", sites);
// // A multisite_ham is simply the corner of a multisite_mpo where the hamiltonian resides
// auto edgeL                 = get_mpo(sites.front()).get_MPO2_edge_left();
// auto edgeR                 = get_mpo(sites.back()).get_MPO2_edge_right();
// auto multisite_mpo_squared = get_multisite_mpo_squared(sites, nbody);
// auto rows                  = multisite_mpo_squared.dimension(2);
// auto cols                  = multisite_mpo_squared.dimension(3);
// tools::log->info("multisite_ham_squared: dims {} | edgeL {} | edgeR {}", multisite_mpo_squared.dimensions(), edgeL.dimensions(), edgeR.dimensions());
// tools::log->info("multisite_mpo_squared: edgeL \n{}\n", linalg::tensor::to_string(edgeL, 16));
// tools::log->info("multisite_mpo_squared: edgeR \n{}\n", linalg::tensor::to_string(edgeR, 16));
// tools::log->info("multisite_mpo_squared: ham2  \n{}\n", linalg::tensor::to_string(multisite_mpo_squared.real(), 16));

// Eigen::Tensor<cplx, 2> ham2 = edgeL.reshape(tenx::array3{1, 1, 1})
//                                   .contract(multisite_mpo_squared, tenx::idx({2}, {0}))
//                                   .contract(edgeR.reshape(tenx::array3{1, 1, 1}), tenx::idx({2}, {2}))
//                                   .shuffle(tenx::array6{3, 1, 5, 2, 0, 4})
//                                   .reshape(tenx::array2{rows, cols});
// return ham2;
// return edgeL.contract(multisite_mpo_squared, tenx::idx({0}, {0})).contract(edgeR, tenx::idx({0}, {0})).reshape(tenx::array2{rows, cols});
// }

const Eigen::Tensor<cplx, 4> &ModelFinite::get_multisite_mpo_squared() const {
    if(cache.multisite_mpo_squared and cache.multisite_mpo_squared_ids and not active_sites.empty()) {
        if constexpr(settings::debug) tools::log->trace("get_multisite_mpo_squared: cache hit");
        // Check that the ids match
        auto active_ids_sq = get_active_ids_sq();
        if(active_ids_sq != cache.multisite_mpo_squared_ids)
            throw except::runtime_error("get_multisite_mpo_squared: cache has mismatching ids: active ids: {} | cache ids {}", active_ids_sq,
                                        cache.multisite_mpo_squared_ids.value());
        return cache.multisite_mpo_squared.value();
    }
    cache.multisite_mpo_squared     = get_multisite_mpo_squared(active_sites);
    cache.multisite_mpo_squared_ids = get_active_ids_sq();
    return cache.multisite_mpo_squared.value();
}

const Eigen::Tensor<cplx, 2> &ModelFinite::get_multisite_ham_squared() const {
    if(cache.multisite_ham_squared and not active_sites.empty()) return cache.multisite_ham_squared.value();
    cache.multisite_ham_squared = get_multisite_ham_squared(active_sites);
    return cache.multisite_ham_squared.value();
}

void ModelFinite::clear_cache(LogPolicy logPolicy) const {
    if(logPolicy == LogPolicy::NORMAL) tools::log->trace("Clearing model cache");
    cache = Cache();
}

std::vector<size_t> ModelFinite::get_active_ids() const {
    std::vector<size_t> ids;
    ids.reserve(active_sites.size());
    for(const auto &pos : active_sites) ids.emplace_back(get_mpo(pos).get_unique_id());
    return ids;
}

std::vector<size_t> ModelFinite::get_active_ids_sq() const {
    std::vector<size_t> ids;
    ids.reserve(active_sites.size());
    for(const auto &pos : active_sites) ids.emplace_back(get_mpo(pos).get_unique_id_sq());
    return ids;
}