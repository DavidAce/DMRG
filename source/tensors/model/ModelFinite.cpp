#include "math/tenx.h"
// -- (textra first)
#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "math/eig.h"
#include "math/linalg/tensor.h"
#include "math/svd.h"
#include "ModelFinite.h"
#include "ModelLocal.h"
#include "qm/spin.h"
#include "tensors/site/mpo/MpoFactory.h"
#include "tid/tid.h"
#include "tools/finite/multisite.h"

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
ModelFinite::~ModelFinite()                              = default; // default dtor
ModelFinite::ModelFinite(ModelFinite &&other)            = default; // default move ctor
ModelFinite &ModelFinite::operator=(ModelFinite &&other) = default; // default move assign

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
    return **std::next(MPO.begin(), static_cast<long>(pos));
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
bool ModelFinite::is_shifted() const {
    bool shifted = MPO.front()->is_shifted();
    for(const auto &mpo : MPO)
        if(shifted != mpo->is_shifted())
            throw std::runtime_error(
                fmt::format("First MPO has is_shifted: {}, but MPO at pos {} has is_shifted: {}", shifted, mpo->get_position(), mpo->is_shifted()));
    return shifted;
}

bool ModelFinite::is_compressed_mpo_squared() const {
    bool compressed = MPO.front()->is_compressed_mpo_squared();
    for(const auto &mpo : MPO)
        if(compressed != mpo->is_compressed_mpo_squared())
            throw except::runtime_error("First MPO has is_compressed_mpo_squared: {}, but MPO at pos {} has is_compressed_mpo_squared: {}", compressed,
                                        mpo->get_position(), mpo->is_compressed_mpo_squared());
    return compressed;
}

double ModelFinite::get_energy_shift() const { return get_energy_shift_per_site() * static_cast<double>(get_length()); }

double ModelFinite::get_energy_shift_per_site() const {
    // Check that all energies are the same
    double e_shift = MPO.front()->get_energy_shift();
    for(const auto &mpo : MPO)
        if(mpo->get_energy_shift() != e_shift) throw std::runtime_error("Shifted energy mismatch!");
    return e_shift;
}

std::vector<std::any> ModelFinite::get_parameter(const std::string &fieldname) {
    std::vector<std::any> fields;
    for(const auto &mpo : MPO) { fields.emplace_back(mpo->get_parameter(fieldname)); }
    return fields;
}

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

void ModelFinite::compress_mpo_squared() {
    cache.multisite_mpo_squared = std::nullopt;
    cache.multisite_ham_squared = std::nullopt;
    auto mpo_squared_compressed = get_compressed_mpo_squared();
    for(const auto &[pos, mpo] : iter::enumerate(MPO)) mpo->set_mpo_squared(mpo_squared_compressed[pos]);
}

// void ModelFinite::sqrt_mpo_squared() {
//     cache.multisite_mpo_squared = std::nullopt;
//     cache.multisite_ham_squared = std::nullopt;
//     auto mpo_squared_sqrt       = get_sqrt_mpo_squared();
//     for(const auto &[pos, mpo] : iter::enumerate(MPO)) mpo->set_mpo_squared(mpo_squared_sqrt[pos]);
// }

bool ModelFinite::has_mpo_squared() const {
    return std::all_of(MPO.begin(), MPO.end(), [](const auto &mpo) { return mpo->has_mpo_squared(); });
}

std::vector<Eigen::Tensor<cplx, 4>> ModelFinite::get_compressed_mpo_squared() {
    tools::log->trace("Compressing MPO²: {} sites", MPO.size());
    if(not has_mpo_squared()) build_mpo_squared(); // Make sure they exist.
    // Collect all the mpo² (doesn't matter if they are already compressed)
    std::vector<Eigen::Tensor<cplx, 4>> mpos_sq;
    mpos_sq.reserve(MPO.size());
    for(const auto &mpo : MPO) mpos_sq.emplace_back(mpo->MPO2());

    // Setup SVD
    // Here we need a lot of precision:
    //  - Use very low svd threshold
    //  - Force the use of JacobiSVD by setting the switchsize_bdc to something large
    //  - Force the use of Lapacke -- it is more precise than Eigen (I don't know why)
    auto svd_cfg = svd::config();
    // Eigen Jacobi becomes ?gesvd (i.e. using QR) with the BLAS backend.
    // See here: https://eigen.tuxfamily.org/bz/show_bug.cgi?id=1732
    svd_cfg.svd_lib = svd::lib::lapacke;
    svd_cfg.svd_rtn = svd::rtn::gesvj;

    svd::solver svd(svd_cfg);

    // Print the results
    std::vector<std::string> report;
    if(tools::log->level() == spdlog::level::trace)
        for(const auto &mpo : mpos_sq) report.emplace_back(fmt::format("{}", mpo.dimensions()));

    for(size_t iter = 0; iter < 4; iter++) {
        // Next compress from left to right
        Eigen::Tensor<cplx, 2> T_l2r; // Transfer matrix
        Eigen::Tensor<cplx, 4> T_mpo_sq;
        for(const auto &[idx, mpo_sq] : iter::enumerate(mpos_sq)) {
            auto mpo_sq_dim_old = mpo_sq.dimensions();
            if(T_l2r.size() == 0)
                T_mpo_sq = mpo_sq;
            else
                T_mpo_sq = T_l2r.contract(mpo_sq, tenx::idx({1}, {0}));

            if(idx == mpos_sq.size() - 1) {
                mpo_sq = T_mpo_sq;
            } else {
                std::tie(mpo_sq, T_l2r) = svd.split_mpo_l2r(T_mpo_sq);
                if(idx + 1 == mpos_sq.size())
                    // The remaining transfer matrix T can be multiplied back into the last MPO from the right
                    mpo_sq = Eigen::Tensor<cplx, 4>(mpo_sq.contract(T_l2r, tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2}));
            }
            if constexpr(settings::debug) tools::log->trace("iter {} | idx {} | dim {} -> {}", iter, idx, mpo_sq_dim_old, mpo_sq.dimensions());
        }

        // Now we have done left to right. Next we do right to left
        Eigen::Tensor<cplx, 2> T_r2l;    // Transfer matrix
        Eigen::Tensor<cplx, 4> mpo_sq_T; // Absorbs transfer matrix
        for(const auto &[idx, mpo_sq] : iter::enumerate_reverse(mpos_sq)) {
            auto mpo_sq_dim_old = mpo_sq.dimensions();
            if(T_r2l.size() == 0)
                mpo_sq_T = mpo_sq;
            else
                mpo_sq_T = mpo_sq.contract(T_r2l, tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2});
            if(idx == 0) {
                mpo_sq = mpo_sq_T;
            } else {
                std::tie(T_r2l, mpo_sq) = svd.split_mpo_r2l(mpo_sq_T);
                if(idx == 0)
                    // The remaining transfer matrix T can be multiplied back into the first MPO from the left
                    mpo_sq = Eigen::Tensor<cplx, 4>(T_r2l.contract(mpo_sq, tenx::idx({1}, {0})));
            }
            if constexpr(settings::debug) tools::log->trace("iter {} | idx {} | dim {} -> {}", iter, idx, mpo_sq_dim_old, mpo_sq.dimensions());
        }
    }

    // Print the results
    if(tools::log->level() == spdlog::level::trace)
        for(const auto &[idx, msg] : iter::enumerate(report)) tools::log->debug("mpo² {}: {} -> {}", idx, msg, mpos_sq[idx].dimensions());

    return mpos_sq;
}

// std::vector<Eigen::Tensor<cplx, 4>> ModelFinite::get_sqrt_mpo_squared() {
//     if(is_compressed_mpo_squared()) throw except::logic_error("Can only take sqrt of non-compressed MPO²");
//     tools::log->trace("Taking square root of MPO²: {} sites", MPO.size());
//
//     // Setup SVD
//     // Here we need a lot of precision:
//     //  - Use very low svd threshold
//     //  - Force the use of JacobiSVD by setting the switchsize_bdc to something large
//     //  - Force the use of Lapacke -- it is more precise than Eigen (I don't know why)
//     auto svd_cfg = svd::config();
//     // Eigen Jacobi becomes ?gesvd (i.e. using QR) with the BLAS backend.
//     // See here: https://eigen.tuxfamily.org/bz/show_bug.cgi?id=1732
//     svd_cfg.svd_lib        = svd::lib::lapacke;
//     svd_cfg.switchsize_bdc = 4096;
//     svd_cfg.use_bdc        = false;
//
//     svd::solver svd(svd_cfg);
//
//     // Setup EIG
//     auto solver = eig::solver();
//
//     // Allocate resulting sqrt mpo²
//     std::vector<Eigen::Tensor<cplx, 4>> mpos_squared_sqrt;
//     mpos_squared_sqrt.reserve(MPO.size());
//
//     for(const auto &[pos, mpo] : iter::enumerate(MPO)) {
//         const auto &d = mpo->MPO().dimensions();
//         // Set up the re-shape and shuffles to convert mpo² to a rank2-tensor
//         auto shf6 = std::array<long, 6>{0, 2, 4, 1, 3, 5};
//         auto shp6 = std::array<long, 6>{d[0], d[0], d[1], d[1], d[2], d[3]};
//         auto shp2 = std::array<long, 2>{d[0] * d[1] * d[2], d[0] * d[1] * d[3]};
//         auto shp4 = std::array<long, 4>{d[0] * d[0], d[1] * d[1], d[2], d[3]};
//
//         Eigen::Tensor<cplx, 2> mpo_squared_matrix = mpo->MPO2().reshape(shp6).shuffle(shf6).reshape(shp2);
//
//         auto [U, S, VT] = svd.decompose(mpo_squared_matrix);
//         shf6            = {0, 3, 1, 4, 2, 5};
//         shp6            = {d[0], d[1], d[2], d[0], d[1], d[3]};
//
//         Eigen::Tensor<cplx, 4> mpo_squared_sqrt =
//             U.contract(tenx::asDiagonal(S.sqrt()), tenx::idx({1}, {0})).contract(VT, tenx::idx({1}, {0})).reshape(shp6).shuffle(shf6).reshape(shp4);
//         mpos_squared_sqrt.emplace_back(mpo_squared_sqrt);
//
//         solver.eig<eig::Form::SYMM>(mpo_squared_matrix.data(), mpo_squared_matrix.dimension(0), eig::Vecs::ON);
//         auto D           = tenx::TensorCast(eig::view::get_eigvals<real>(solver.result).cast<cplx>());
//         auto V           = tenx::TensorCast(eig::view::get_eigvecs<cplx>(solver.result));
//         mpo_squared_sqrt = V.contract(tenx::asDiagonal(D.sqrt()), tenx::idx({1}, {0}))
//                                .contract(V.shuffle(std::array<long, 2>{1, 0}), tenx::idx({1}, {0}))
//                                .reshape(shp6)
//                                .shuffle(shf6)
//                                .reshape(shp4);
//
//         Eigen::Tensor<double,1> S_real = S.real();
//         Eigen::Tensor<double,1> D_real = D.real();
//         tools::log->info("mpo²[{}] svds: {:8.4e}", static_cast<long>(pos), fmt::join(tenx::span(S_real),", "));
//         tools::log->info("mpo²[{}] eigv: {:8.4e}", static_cast<long>(pos), fmt::join(tenx::span(D_real),", "));
//         mpos_squared_sqrt.emplace_back(mpo_squared_sqrt);
//     }
//
//     return mpos_squared_sqrt;
// }

void ModelFinite::set_energy_shift(double total_energy) { set_energy_shift_per_site(total_energy / static_cast<double>(get_length())); }

void ModelFinite::set_energy_shift_per_site(double energy_shift_per_site) {
    if(get_energy_shift_per_site() == energy_shift_per_site) return;
    tools::log->debug("Shifting MPO energy per site: {:.16f}", energy_shift_per_site);
    for(const auto &mpo : MPO) mpo->set_energy_shift(energy_shift_per_site);
    clear_cache();
}

void ModelFinite::set_psfactor(double psfactor) {
    for(const auto &mpo : MPO) mpo->set_psfactor(psfactor);
    clear_cache();
}

bool ModelFinite::set_parity_shift_mpo_squared(int sign, std::string_view axis) {
    if(not qm::spin::half::is_valid_axis(axis)) return false;
    if(get_parity_shift_mpo_squared() == std::make_pair(sign, axis)) return false;
    tools::log->info("Setting MPO² parity shift for axis {}{}", sign == 0 ? "" : (sign < 0 ? "-" : "+"), qm::spin::half::get_axis_unsigned(axis));
    for(const auto &mpo : MPO) mpo->set_parity_shift_mpo_squared(sign, axis);
    clear_cache();
    return true;
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

bool ModelFinite::has_parity_shift_mpo_squared() const {
    auto parity_shift = get_parity_shift_mpo_squared();
    return parity_shift.first != 0 and not parity_shift.second.empty();
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

Eigen::Tensor<cplx, 4> ModelFinite::get_multisite_mpo(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody) const {
    // Observe that nbody empty/nullopt have very different meanings
    //      - empty means that no interactions should be taken into account, effectively setting all J(i,j...) = 0
    //      - nullopt means that we want the default mpo with (everything on)
    //      - otherwise nbody with values like {1,2} would imply we want 1 and 2-body interactions turned on
    //      - if nbody has a 0 value in it, it means we want to make an attempt to account for double-counting in multisite mpos.

    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mpo tensor");
    if(sites == active_sites and cache.multisite_mpo and not nbody) return cache.multisite_mpo.value();
    if(not nbody)
        tools::log->trace("Contracting multisite mpo tensor with sites {}", sites);
    else
        tools::log->trace("Contracting multisite mpo tensor with sites {} | nbody {} ", sites, nbody.value());

    auto                   t_mpo        = tid::tic_scope("get_multisite_mpo", tid::level::highest);
    constexpr auto         shuffle_idx  = tenx::array6{0, 3, 1, 4, 2, 5};
    constexpr auto         contract_idx = tenx::idx({1}, {0});
    auto                   positions    = num::range<size_t>(sites.front(), sites.back() + 1);
    auto                   skip         = std::optional<std::vector<size_t>>();
    auto                   skip_log     = std::vector<size_t>();
    bool                   do_cache     = nbody.has_value() and nbody->back() > 1; // Caching doesn't make sense for nbody == 1
    Eigen::Tensor<cplx, 4> multisite_mpo, mpoL, mpoR;
    Eigen::Tensor<cplx, 2> mpoR_traced;

    if(sites != positions) {
        skip = std::vector<size_t>{};
        for(const auto &pos : positions) {
            if(std::find(sites.begin(), sites.end(), pos) == sites.end()) skip->emplace_back(pos);
        }
    }

    for(const auto &pos : positions) {
        // sites needs to be sorted, but may skip sites.
        // For instance, sites == {3,9} is valid. Then sites 4,5,6,7,8 are skipped.
        // When a site is skipped, we set the contribution from its interaction terms to zero and trace over it so that
        // the physical dimension doesn't grow.
        bool do_trace = skip.has_value() and std::find(skip->begin(), skip->end(), pos) != skip->end();
        if(multisite_mpo.size() == 0) {
            auto t_pre = tid::tic_scope("prepending", tid::level::highest);
            if(nbody or skip)
                multisite_mpo = get_mpo(pos).MPO_nbody_view(nbody, skip);
            else
                multisite_mpo = get_mpo(pos).MPO();

            if(do_trace and do_cache) skip_log.emplace_back(pos);
            continue;
        }

        if(nbody or skip) {
            mpoL = multisite_mpo;
            mpoR = get_mpo(pos).MPO_nbody_view(nbody, skip);
        } else {
            mpoL = multisite_mpo;
            mpoR = get_mpo(pos).MPO();
        }

        // Determine if this position adds to the physical dimension or if it will get traced over
        long dim0     = mpoL.dimension(0);
        long dim1     = mpoR.dimension(1);
        long dim2     = mpoL.dimension(2) * (do_trace ? 1l : mpoR.dimension(2));
        long dim3     = mpoL.dimension(3) * (do_trace ? 1l : mpoR.dimension(3));
        auto new_dims = std::array<long, 4>{dim0, dim1, dim2, dim3};
        multisite_mpo.resize(new_dims);
        auto new_cache_string = fmt::format("siteL{}|pos{}|dims{}|skips{}|trace{}", sites.front(), pos, new_dims, skip_log, do_trace);
        if(do_cache and cache.multisite_mpo_temps.find(new_cache_string) != cache.multisite_mpo_temps.end()) {
            multisite_mpo = cache.multisite_mpo_temps.at(new_cache_string);
        } else {
            if(do_trace) {
                auto t_skip = tid::tic_scope("skipping", tid::level::highest);
                // Trace the physical indices of this skipped mpo (this should trace an identity)
                mpoR_traced = mpoR.trace(tenx::array2{2, 3});
                mpoR_traced *= mpoR_traced.constant(0.5); // divide by 2 (after tracing identity)
                // Append it to the multisite mpo
                multisite_mpo.device(tenx::threads::getDevice()) =
                    mpoL.contract(mpoR_traced, tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2}).reshape(new_dims);
            } else {
                auto t_app                                       = tid::tic_scope("appending", tid::level::highest);
                multisite_mpo.device(tenx::threads::getDevice()) = mpoL.contract(mpoR, contract_idx).shuffle(shuffle_idx).reshape(new_dims);
            }
            if(do_cache) cache.multisite_mpo_temps[new_cache_string] = multisite_mpo;
        }
        if(do_trace and do_cache) skip_log.emplace_back(pos);
    }
    return multisite_mpo;
}

Eigen::Tensor<cplx_t, 4> ModelFinite::get_multisite_mpo_t(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody) const {
    // Observe that nbody empty/nullopt have very different meanings
    //      - empty means that no interactions should be taken into account, effectively setting all J(i,j...) = 0
    //      - nullopt means that we want the default mpo with (everything on)
    //      - otherwise nbody with values like {1,2} would imply we want 1 and 2-body interactions turned on
    //      - if nbody has a 0 value in it, it means we want to make an attempt to account for double-counting in multisite mpos.

    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mpo tensor");
    if(sites == active_sites and cache.multisite_mpo_t and not nbody) return cache.multisite_mpo_t.value();
    if(not nbody)
        tools::log->trace("Contracting multisite mpo tensor with sites {}", sites);
    else
        tools::log->trace("Contracting multisite mpo tensor with sites {} | nbody {} ", sites, nbody.value());

    auto                     t_mpo        = tid::tic_scope("get_multisite_mpo_t", tid::level::highest);
    constexpr auto           shuffle_idx  = tenx::array6{0, 3, 1, 4, 2, 5};
    constexpr auto           contract_idx = tenx::idx({1}, {0});
    auto                     positions    = num::range<size_t>(sites.front(), sites.back() + 1);
    auto                     skip         = std::optional<std::vector<size_t>>();
    auto                     skip_log     = std::vector<size_t>();
    bool                     do_cache     = nbody.has_value() and nbody->back() > 1; // Caching doesn't make sense for nbody == 1
    Eigen::Tensor<cplx_t, 4> multisite_mpo_t, mpoL, mpoR;
    Eigen::Tensor<cplx_t, 2> mpoR_traced;

    if(sites != positions) {
        skip = std::vector<size_t>{};
        for(const auto &pos : positions) {
            if(std::find(sites.begin(), sites.end(), pos) == sites.end()) skip->emplace_back(pos);
        }
    }

    for(const auto &pos : positions) {
        // sites needs to be sorted, but may skip sites.
        // For instance, sites == {3,9} is valid. Then sites 4,5,6,7,8 are skipped.
        // When a site is skipped, we set the contribution from its interaction terms to zero and trace over it so that
        // the physical dimension doesn't grow.
        bool do_trace = skip.has_value() and std::find(skip->begin(), skip->end(), pos) != skip->end();
        if(multisite_mpo_t.size() == 0) {
            auto t_pre = tid::tic_scope("prepending", tid::level::highest);
            if(nbody or skip)
                multisite_mpo_t = get_mpo(pos).MPO_nbody_view_t(nbody, skip);
            else
                multisite_mpo_t = get_mpo(pos).MPO_t();
            if(do_trace and do_cache) skip_log.emplace_back(pos);
            continue;
        }

        if(nbody or skip) {
            mpoL = multisite_mpo_t;
            mpoR = get_mpo(pos).MPO_nbody_view_t(nbody, skip);
        } else {
            mpoL = multisite_mpo_t;
            mpoR = get_mpo(pos).MPO_t();
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
        auto new_cache_string = fmt::format("siteL{}|pos{}|dims{}|skips{}|trace{}", sites.front(), pos, new_dims, skip_log, do_trace);
        if(do_cache and cache.multisite_mpo_t_temps.find(new_cache_string) != cache.multisite_mpo_t_temps.end()) {
            multisite_mpo_t = cache.multisite_mpo_t_temps.at(new_cache_string);
        } else {
            if(do_trace) {
                auto t_skip = tid::tic_scope("skipping", tid::level::highest);
                // Trace the physical indices of this skipped mpo (this should trace an identity)
                mpoR_traced = mpoR.trace(tenx::array2{2, 3});
                mpoR_traced *= mpoR_traced.constant(0.5); // divide by 2 (after tracing identity)
                // Append it to the multisite mpo
                multisite_mpo_t.device(tenx::threads::getDevice()) =
                    mpoL.contract(mpoR_traced, tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2}).reshape(new_dims);
            } else {
                auto t_app                                         = tid::tic_scope("appending", tid::level::highest);
                multisite_mpo_t.device(tenx::threads::getDevice()) = mpoL.contract(mpoR, contract_idx).shuffle(shuffle_idx).reshape(new_dims);
            }
            if(do_cache) cache.multisite_mpo_t_temps[new_cache_string] = multisite_mpo_t;
        }
        if(do_trace and do_cache) skip_log.emplace_back(pos);
        // This intermediate multisite_mpo_t could be the result we are looking for at a later time, so cache it!
    }
    return multisite_mpo_t;
}

Eigen::Tensor<cplx, 2> ModelFinite::get_multisite_ham(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite hamiltonian tensor");
    if(sites == active_sites and cache.multisite_ham and not nbody) return cache.multisite_ham.value();
    // A multisite_ham is simply the corner of a multisite_mpo where the hamiltonian resides
    auto multisite_mpo = get_multisite_mpo(sites, nbody);
    auto edgeL         = get_mpo(sites.front()).get_MPO_edge_left();
    auto edgeR         = get_mpo(sites.back()).get_MPO_edge_right();
    return multisite_mpo.contract(edgeL, tenx::idx({0}, {0}))
        .contract(edgeR, tenx::idx({0}, {0}))
        .reshape(tenx::array2{multisite_mpo.dimension(2), multisite_mpo.dimension(3)});
}

Eigen::Tensor<cplx_t, 2> ModelFinite::get_multisite_ham_t(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite hamiltonian tensor");
    if(sites == active_sites and cache.multisite_ham_t and not nbody) return cache.multisite_ham_t.value();
    // A multisite_ham is simply the corner of a multisite_mpo where the hamiltonian resides
    auto multisite_mpo_t = get_multisite_mpo_t(sites, nbody);
    auto edgeL           = get_mpo(sites.front()).get_MPO_edge_left();
    auto edgeR           = get_mpo(sites.back()).get_MPO_edge_right();
    return multisite_mpo_t.contract(edgeL.cast<cplx_t>(), tenx::idx({0}, {0}))
        .contract(edgeR.cast<cplx_t>(), tenx::idx({0}, {0}))
        .reshape(tenx::array2{multisite_mpo_t.dimension(2), multisite_mpo_t.dimension(3)});
}

const Eigen::Tensor<cplx, 4> &ModelFinite::get_multisite_mpo() const {
    if(cache.multisite_mpo and not active_sites.empty()) return cache.multisite_mpo.value();
    cache.multisite_mpo = get_multisite_mpo(active_sites);
    return cache.multisite_mpo.value();
}
const Eigen::Tensor<cplx_t, 4> &ModelFinite::get_multisite_mpo_t() const {
    if(cache.multisite_mpo_t and not active_sites.empty()) return cache.multisite_mpo_t.value();
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
    for(const auto &site : active_sites) {
        if(multisite_mpo.size() == 0) {
            multisite_mpo = get_mpo(site).MPO_shifted_view(energy_per_site);
            continue;
        }
        const auto         &mpo      = get_mpo(site);
        long                dim0     = multisite_mpo.dimension(0);
        long                dim1     = mpo.MPO().dimension(1);
        long                dim2     = multisite_mpo.dimension(2) * mpo.MPO().dimension(2);
        long                dim3     = multisite_mpo.dimension(3) * mpo.MPO().dimension(3);
        std::array<long, 4> new_dims = {dim0, dim1, dim2, dim3};
        temp.resize(new_dims);
        temp.device(tenx::threads::getDevice()) =
            multisite_mpo.contract(mpo.MPO_shifted_view(energy_per_site), contract_idx).shuffle(shuffle_idx).reshape(new_dims);
        multisite_mpo = temp;
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
    multisite_mpo_squared_shifted.device(tenx::threads::getDevice()) =
        multisite_mpo_shifted.contract(multisite_mpo_shifted, tenx::idx({3}, {2})).shuffle(tenx::array6{0, 3, 1, 4, 2, 5}).reshape(mpo_squared_dims);
    return multisite_mpo_squared_shifted;
}

std::array<long, 4> ModelFinite::active_dimensions_squared() const { return tools::finite::multisite::get_dimensions_squared(*this); }

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
    bool                   first = true;
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
            temp.device(tenx::threads::getDevice()) =
                multisite_mpo_squared.contract(mpo.MPO2_nbody_view(nbody_local), contract_idx).shuffle(shuffle_idx).reshape(new_dims);
        else
            temp.device(tenx::threads::getDevice()) = multisite_mpo_squared.contract(mpo.MPO2(), contract_idx).shuffle(shuffle_idx).reshape(new_dims);

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
            long d0 = dim0;
            long d1 = dim1;
            long d2 = multisite_mpo_squared.dimension(2);
            long d3 = mpo.MPO().dimension(2);
            long d4 = multisite_mpo_squared.dimension(3);
            long d5 = mpo.MPO().dimension(3);

            Eigen::Tensor<cplx, 4> temp2 = temp.reshape(tenx::array6{d0, d1, d2, d3, d4, d5}).trace(tenx::array2{3, 5});
            multisite_mpo_squared        = temp2 * temp2.constant(0.5);
        } else {
            multisite_mpo_squared = temp;
        }
    }
    return multisite_mpo_squared;
}

Eigen::Tensor<cplx, 2> ModelFinite::get_multisite_ham_squared(const std::vector<size_t> &sites, std::optional<std::vector<size_t>> nbody) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite hamiltonian tensor");
    if(sites == active_sites and cache.multisite_ham and not nbody) return cache.multisite_ham.value();
    // A multisite_ham is simply the corner of a multisite_mpo where the hamiltonian resides
    auto edgeL = get_mpo(sites.front()).get_MPO2_edge_left();
    auto edgeR = get_mpo(sites.back()).get_MPO2_edge_right();
    // TODO Fix the code below

    //    if(sites == num::range<size_t>(0, get_length()) and not nbody) {
    //        // We want the full hamiltonian. We save time by not including the mpo edges.
    //        constexpr auto shuffle_idx  = tenx::array6{0, 3, 1, 4, 2, 5};
    //        constexpr auto contract_idx = tenx::idx({1}, {0});
    //        auto           bulk_sites   = num::range<size_t>(1, get_length() - 1);
    //        auto           mpo_front    = get_mpo(sites.front()).MPO2();
    //        auto           dim_front    = mpo_front.dimensions();
    //        auto           mpo_back     = get_mpo(sites.back()).MPO2();
    //        auto           dim_back     = mpo_back.dimensions();
    //        dim_front[0]                = edgeL.dimension(0) * edgeL.dimension(1);
    //        dim_back[0]                 = edgeR.dimension(0) * edgeR.dimension(1);
    //
    //        Eigen::Tensor<cplx, 4> mpo_edgeL = edgeL.contract(mpo_front, tenx::idx({2}, {0})).reshape(mpo_front.dimensions());
    //        Eigen::Tensor<cplx, 4> mpo_edgeR = edgeR.contract(mpo_back, tenx::idx({2}, {1})).shuffle(tenx::array5{2, 0, 1, 3,
    //        4}).reshape(mpo_back.dimensions()); Eigen::Tensor<cplx, 4> mpo_temp; for(const auto &pos : bulk_sites) {
    //            const auto  &mpo      = get_mpo(pos);
    //            long         dim0     = mpo_edgeL.dimension(0);
    //            long         dim1     = mpo.MPO2().dimension(1);
    //            long         dim2     = mpo_edgeL.dimension(2) * mpo.MPO2().dimension(2);
    //            long         dim3     = mpo_edgeL.dimension(3) * mpo.MPO2().dimension(3);
    //            tenx::array4 new_dims = {dim0, dim1, dim2, dim3};
    //            mpo_temp.resize(new_dims);
    //            mpo_temp.device(tenx::omp::getDevice()) = mpo_edgeL.contract(mpo.MPO2(), contract_idx).shuffle(shuffle_idx).reshape(new_dims);
    //            mpo_edgeL                               = mpo_temp;
    //        }
    //        long         dim0     = mpo_edgeL.dimension(0);
    //        long         dim1     = mpo_edgeR.dimension(1);
    //        long         dim2     = mpo_edgeL.dimension(2) * mpo_edgeR.dimension(2);
    //        long         dim3     = mpo_edgeL.dimension(3) * mpo_edgeR.dimension(3);
    //        tenx::array4 new_dims = {dim0, dim1, dim2, dim3};
    //        mpo_temp.resize(new_dims);
    //        mpo_temp.device(tenx::omp::getDevice()) = mpo_edgeL.contract(mpo_edgeR, contract_idx).shuffle(shuffle_idx).reshape(new_dims);
    //        return mpo_temp;
    //    }
    auto multisite_mpo_squared = get_multisite_mpo_squared(sites, nbody);
    return multisite_mpo_squared.contract(edgeL, tenx::idx({0}, {0}))
        .contract(edgeR, tenx::idx({0}, {0}))
        .reshape(tenx::array2{multisite_mpo_squared.dimension(2), multisite_mpo_squared.dimension(3)});
}

const Eigen::Tensor<cplx, 4> &ModelFinite::get_multisite_mpo_squared() const {
    if(cache.multisite_mpo_squared) return cache.multisite_mpo_squared.value();
    cache.multisite_mpo_squared = get_multisite_mpo_squared(active_sites);
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