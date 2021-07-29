#include <math/tenx.h>
// -- (textra first)
#include "ModelFinite.h"
#include <config/settings.h>
#include <general/iter.h>
#include <math/svd.h>
#include <tensors/site/mpo/MpoFactory.h>
#include <tid/tid.h>
#include <tools/finite/multisite.h>

ModelFinite::ModelFinite() = default; // Can't initialize lists since we don't know the model size yet

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
ModelFinite::~ModelFinite()                   = default;            // default dtor
ModelFinite::ModelFinite(ModelFinite &&other) = default;            // default move ctor
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
    tools::log->info("Initializing model {} with {} sites", enum2str(model_type_), model_size);
    if(model_size < 2) throw std::logic_error("Tried to initialize model with less than 2 sites");
    if(model_size > 2048) throw std::logic_error("Tried to initialize model with more than 2048 sites");
    if(not MPO.empty()) throw std::logic_error("Tried to initialize over an existing model. This is usually not what you want!");
    // Generate MPO
    model_type = model_type_;
    for(size_t site = 0; site < model_size; site++) { MPO.emplace_back(MpoFactory::create_mpo(site, model_type)); }
    if(MPO.size() != model_size) throw std::logic_error("Initialized MPO with wrong size");
}

const MpoSite &ModelFinite::get_mpo(size_t pos) const {
    if(pos >= MPO.size()) throw std::range_error(fmt::format("get_mpo(pos) pos out of range: {}", pos));
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
}

// For reduced energy MPO's

bool ModelFinite::is_reduced() const {
    bool reduced = MPO.front()->is_reduced();
    for(const auto &mpo : MPO)
        if(reduced != mpo->is_reduced())
            throw std::runtime_error(
                fmt::format("First MPO has is_reduced: {}, but MPO at pos {} has is_reduced: {}", reduced, mpo->get_position(), mpo->is_reduced()));
    return reduced;
}

bool ModelFinite::is_compressed_mpo_squared() const {
    bool compressed = MPO.front()->is_compressed_mpo_squared();
    for(const auto &mpo : MPO)
        if(compressed != mpo->is_compressed_mpo_squared())
            throw std::runtime_error(fmt::format("First MPO has is_compressed_mpo_squared: {}, but MPO at pos {} has is_compressed_mpo_squared: {}", compressed,
                                                 mpo->get_position(), mpo->is_compressed_mpo_squared()));
    return compressed;
}

double ModelFinite::get_energy_reduced() const { return get_energy_per_site_reduced() * static_cast<double>(get_length()); }

double ModelFinite::get_energy_per_site_reduced() const {
    // Check that all energies are the same
    double e_reduced = MPO.front()->get_reduced_energy();
    for(const auto &mpo : MPO)
        if(mpo->get_reduced_energy() != e_reduced) throw std::runtime_error("Reduced energy mismatch!");
    return e_reduced;
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

bool ModelFinite::has_mpo_squared() const {
    return std::all_of(MPO.begin(), MPO.end(), [](const auto &mpo) { return mpo->has_mpo_squared(); });
}

void ModelFinite::reset_mpo_squared() {
    tools::log->debug("Resetting MPO²");
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

void ModelFinite::compress_mpo_squared(std::optional<svd::settings> svd_settings) {
    tools::log->debug("Compressing MPO²");
    cache.multisite_mpo_squared = std::nullopt;
    cache.multisite_ham_squared = std::nullopt;
    auto mpo_compressed         = get_compressed_mpo_squared(svd_settings);
    for(const auto &[pos, mpo] : iter::enumerate(MPO)) mpo->set_mpo_squared(mpo_compressed[pos]);
}

std::vector<Eigen::Tensor<ModelFinite::Scalar, 4>> ModelFinite::get_compressed_mpo_squared(std::optional<svd::settings> svd_settings) {
    // First, rebuild the MPO's
    std::vector<Eigen::Tensor<Scalar, 4>> mpos_sq;
    for(const auto &mpo : MPO) mpos_sq.emplace_back(mpo->get_non_compressed_mpo_squared());
    tools::log->trace("Compressing MPO²");

    // Setup SVD
    // Here we need a lot of precision:
    //  - Use very low svd threshold
    //  - Force the use of JacobiSVD by setting the switchsize to something large
    //  - Force the use of Lapacke -- it is more precise than Eigen (I don't know why)
    if(not svd_settings) svd_settings = svd::settings();
    if(not svd_settings->threshold) svd_settings->threshold = std::numeric_limits<double>::epsilon();
    if(not svd_settings->switchsize) svd_settings->switchsize = 4096;
    if(not svd_settings->use_lapacke) svd_settings->use_lapacke = true;
    if(not svd_settings->use_bdc) svd_settings->use_bdc = false;
    if(not svd_settings->loglevel) svd_settings->loglevel = 2;

    svd::solver svd(svd_settings);

    // Print the results
    std::vector<std::string> report;
    if(tools::log->level() == spdlog::level::trace)
        for(const auto &mpo : mpos_sq) report.emplace_back(fmt::format("{}", mpo.dimensions()));

    for(size_t iter = 0; iter < 1; iter++) {
        // Next compress from left to right
        Eigen::Tensor<Scalar, 2> T_l2r; // Transfer matrix
        Eigen::Tensor<Scalar, 4> T_mpo_sq;
        for(const auto &[idx, mpo_sq] : iter::enumerate(mpos_sq)) {
            if(T_l2r.size() == 0)
                T_mpo_sq = mpo_sq;
            else
                T_mpo_sq = T_l2r.contract(mpo_sq, tenx::idx({1}, {0}));

            if(idx == mpos_sq.size() - 1) {
                mpo_sq = T_mpo_sq;
            } else {
                auto [U, S, V] = svd.split_mpo_l2r(T_mpo_sq);
                T_l2r          = tenx::asDiagonal(S).contract(V, tenx::idx({1}, {0}));
                if(idx < mpos_sq.size() - 1)
                    mpo_sq = U;
                else
                    // The remaining transfer matrix T can be multiplied back into the last MPO from the right
                    mpo_sq = U.contract(T_l2r, tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2});
            }
        }

        // Now we have done left to right. Next we do right to left
        Eigen::Tensor<Scalar, 2> T_r2l;    // Transfer matrix
        Eigen::Tensor<Scalar, 4> mpo_sq_T; // Absorbs transfer matrix
        for(const auto &[idx, mpo_sq] : iter::enumerate_reverse(mpos_sq)) {
            if(T_r2l.size() == 0)
                mpo_sq_T = mpo_sq;
            else
                mpo_sq_T = mpo_sq.contract(T_r2l, tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2});
            if(idx == 0) {
                mpo_sq = mpo_sq_T;
            } else {
                auto [U, S, V] = svd.split_mpo_r2l(mpo_sq_T);
                T_r2l          = U.contract(tenx::asDiagonal(S), tenx::idx({1}, {0}));
                if(idx > 0)
                    mpo_sq = V;
                else
                    // The remaining transfer matrix T can be multiplied back into the first MPO from the left
                    mpo_sq = T_r2l.contract(V, tenx::idx({1}, {0}));
            }
        }
    }

    // Print the results
    if(tools::log->level() == spdlog::level::trace)
        for(const auto &[idx, msg] : iter::enumerate(report)) tools::log->trace("mpo² {}: {} -> {}", idx, msg, mpos_sq[idx].dimensions());

    return mpos_sq;
}

void ModelFinite::set_reduced_energy(double total_energy) { set_reduced_energy_per_site(total_energy / static_cast<double>(get_length())); }

void ModelFinite::set_reduced_energy_per_site(double site_energy) {
    if(get_energy_per_site_reduced() == site_energy) return;
    tools::log->debug("Reducing MPO energy per site by {:.16f}", site_energy);
    for(const auto &mpo : MPO) mpo->set_reduced_energy(site_energy);
    clear_cache();
}

void ModelFinite::perturb_hamiltonian(double coupling_ptb, double field_ptb, PerturbMode perturbMode) {
    std::vector<MpoSite::TableMap> all_params;
    clear_cache();
    for(const auto &mpo : MPO) {
        mpo->set_perturbation(coupling_ptb, field_ptb, perturbMode);
        all_params.push_back(mpo->get_parameters());
    }
    for(const auto &mpo : MPO) mpo->set_averages(all_params, false);
    if(coupling_ptb == 0.0 and field_ptb == 0.0 and is_perturbed()) throw std::runtime_error("Model: Should have unperturbed!");
}

void ModelFinite::damp_model_disorder(double coupling_damp, double field_damp) {
    clear_cache();
    for(const auto &mpo : MPO) {
        mpo->set_coupling_damping(coupling_damp);
        mpo->set_field_damping(field_damp);
    }
    if(coupling_damp == 0.0 and field_damp == 0.0 and is_damped()) throw std::runtime_error("Model: Should have undamped!");
}

bool ModelFinite::is_perturbed() const {
    for(size_t pos = 0; pos < get_length(); pos++) {
        if(get_mpo(pos).is_perturbed()) return true;
    }
    return false;
}

bool ModelFinite::is_damped() const {
    for(size_t pos = 0; pos < get_length(); pos++) {
        if(get_mpo(pos).is_damped()) return true;
    }
    return false;
}

std::array<long, 4> ModelFinite::active_dimensions() const { return tools::finite::multisite::get_dimensions(*this); }

Eigen::Tensor<ModelFinite::Scalar, 4> ModelFinite::get_multisite_mpo(const std::vector<size_t> &sites, const std::vector<size_t> &nbody) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mpo tensor");
    if(sites == active_sites and cache.multisite_mpo and nbody.empty()) return cache.multisite_mpo.value();
    if(nbody.empty())
        tools::log->trace("Contracting multisite mpo tensor with sites {}", sites);
    else
        tools::log->trace("Contracting multisite mpo tensor with sites {} | nbody {} ", sites, nbody);

    auto                     t_mpo = tid::tic_scope("mpo");
    Eigen::Tensor<Scalar, 4> multisite_mpo, temp;
    constexpr auto           shuffle_idx  = tenx::array6{0, 3, 1, 4, 2, 5};
    constexpr auto           contract_idx = tenx::idx({1}, {0});
    for(const auto &site : sites) {
        if(multisite_mpo.size() == 0) {
            if(nbody.empty())
                multisite_mpo = get_mpo(site).MPO();
            else
                multisite_mpo = get_mpo(site).MPO_nbody_view(nbody);
            continue;
        }
        const auto         &mpo      = get_mpo(site);
        long                dim0     = multisite_mpo.dimension(0);
        long                dim1     = mpo.MPO().dimension(1);
        long                dim2     = multisite_mpo.dimension(2) * mpo.MPO().dimension(2);
        long                dim3     = multisite_mpo.dimension(3) * mpo.MPO().dimension(3);
        std::array<long, 4> new_dims = {dim0, dim1, dim2, dim3};
        temp.resize(new_dims);
        if(nbody.empty()) // Avoids creating a temporary
            temp.device(tenx::omp::getDevice()) = multisite_mpo.contract(mpo.MPO(), contract_idx).shuffle(shuffle_idx).reshape(new_dims);
        else
            temp.device(tenx::omp::getDevice()) = multisite_mpo.contract(mpo.MPO_nbody_view(nbody), contract_idx).shuffle(shuffle_idx).reshape(new_dims);
        multisite_mpo = temp;
    }
    return multisite_mpo;
}

Eigen::Tensor<ModelFinite::Scalar, 2> ModelFinite::get_multisite_ham(const std::vector<size_t> &sites, const std::vector<size_t> &nbody_terms) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite hamiltonian tensor");
    if(sites == active_sites and cache.multisite_ham and nbody_terms.empty()) return cache.multisite_ham.value();
    // A multisite_ham is simply the corner of a multisite_mpo where the hamiltonian resides
    auto multisite_mpo = get_multisite_mpo(sites, nbody_terms);
    auto edgeL         = get_mpo(sites.front()).get_MPO_edge_left();
    auto edgeR         = get_mpo(sites.back()).get_MPO_edge_right();
    return multisite_mpo.contract(edgeL, tenx::idx({0}, {0}))
        .contract(edgeR, tenx::idx({0}, {0}))
        .reshape(tenx::array2{multisite_mpo.dimension(2), multisite_mpo.dimension(3)});
}

const Eigen::Tensor<ModelFinite::Scalar, 4> &ModelFinite::get_multisite_mpo() const {
    if(cache.multisite_mpo and not active_sites.empty()) return cache.multisite_mpo.value();
    cache.multisite_mpo = get_multisite_mpo(active_sites);
    return cache.multisite_mpo.value();
}

const Eigen::Tensor<ModelFinite::Scalar, 2> &ModelFinite::get_multisite_ham() const {
    if(cache.multisite_ham and not active_sites.empty()) return cache.multisite_ham.value();
    cache.multisite_ham = get_multisite_ham(active_sites);
    return cache.multisite_ham.value();
}

Eigen::Tensor<ModelFinite::Scalar, 4> ModelFinite::get_multisite_mpo_reduced_view(double energy_per_site) const {
    auto                     t_mpo = tid::tic_scope("mpo");
    Eigen::Tensor<Scalar, 4> multisite_mpo, temp;
    constexpr auto           shuffle_idx  = tenx::array6{0, 3, 1, 4, 2, 5};
    constexpr auto           contract_idx = tenx::idx({1}, {0});
    for(const auto &site : active_sites) {
        if(multisite_mpo.size() == 0) {
            multisite_mpo = get_mpo(site).MPO_reduced_view(energy_per_site);
            continue;
        }
        const auto         &mpo      = get_mpo(site);
        long                dim0     = multisite_mpo.dimension(0);
        long                dim1     = mpo.MPO().dimension(1);
        long                dim2     = multisite_mpo.dimension(2) * mpo.MPO().dimension(2);
        long                dim3     = multisite_mpo.dimension(3) * mpo.MPO().dimension(3);
        std::array<long, 4> new_dims = {dim0, dim1, dim2, dim3};
        temp.resize(new_dims);
        temp.device(tenx::omp::getDevice()) =
            multisite_mpo.contract(mpo.MPO_reduced_view(energy_per_site), contract_idx).shuffle(shuffle_idx).reshape(new_dims);
        multisite_mpo = temp;
    }
    return multisite_mpo;
}

Eigen::Tensor<ModelFinite::Scalar, 4> ModelFinite::get_multisite_mpo_squared_reduced_view(double energy_per_site) const {
    auto                     multisite_mpo_reduced = get_multisite_mpo_reduced_view(energy_per_site);
    long                     dim0                  = multisite_mpo_reduced.dimension(0) * multisite_mpo_reduced.dimension(0);
    long                     dim1                  = multisite_mpo_reduced.dimension(1) * multisite_mpo_reduced.dimension(1);
    long                     dim2                  = multisite_mpo_reduced.dimension(2);
    long                     dim3                  = multisite_mpo_reduced.dimension(3);
    std::array<long, 4>      mpo_squared_dims      = {dim0, dim1, dim2, dim3};
    Eigen::Tensor<Scalar, 4> multisite_mpo_squared_reduced(mpo_squared_dims);
    multisite_mpo_squared_reduced.device(tenx::omp::getDevice()) =
        multisite_mpo_reduced.contract(multisite_mpo_reduced, tenx::idx({3}, {2})).shuffle(tenx::array6{0, 3, 1, 4, 2, 5}).reshape(mpo_squared_dims);
    return multisite_mpo_squared_reduced;
}

std::array<long, 4> ModelFinite::active_dimensions_squared() const { return tools::finite::multisite::get_dimensions_squared(*this); }

Eigen::Tensor<ModelFinite::Scalar, 4> ModelFinite::get_multisite_mpo_squared(const std::vector<size_t> &sites, const std::vector<size_t> &nbody) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mpo squared tensor");
    if(sites == active_sites and cache.multisite_mpo_squared) return cache.multisite_mpo_squared.value();
    tools::log->trace("Contracting multisite mpo squared tensor with {} sites", sites.size());
    auto                     t_mpo = tid::tic_scope("mpo");
    Eigen::Tensor<Scalar, 4> multisite_mpo_squared;
    constexpr auto           shuffle_idx  = tenx::array6{0, 3, 1, 4, 2, 5};
    constexpr auto           contract_idx = tenx::idx({1}, {0});
    tenx::array4             new_dims;
    Eigen::Tensor<Scalar, 4> temp;
    bool                     first = true;
    for(const auto &site : sites) {
        if(first) {
            if(nbody.empty())
                multisite_mpo_squared = get_mpo(site).MPO2();
            else
                multisite_mpo_squared = get_mpo(site).MPO2_nbody_view(nbody);
            first = false;
            continue;
        }

        const auto &mpo  = get_mpo(site);
        long        dim0 = multisite_mpo_squared.dimension(0);
        long        dim1 = mpo.MPO2().dimension(1);
        long        dim2 = multisite_mpo_squared.dimension(2) * mpo.MPO2().dimension(2);
        long        dim3 = multisite_mpo_squared.dimension(3) * mpo.MPO2().dimension(3);
        new_dims         = {dim0, dim1, dim2, dim3};
        temp.resize(new_dims);
        if(nbody.empty()) // Avoids creating a temporary
            temp.device(tenx::omp::getDevice()) = multisite_mpo_squared.contract(mpo.MPO2(), contract_idx).shuffle(shuffle_idx).reshape(new_dims);
        else
            temp.device(tenx::omp::getDevice()) =
                multisite_mpo_squared.contract(mpo.MPO2_nbody_view(nbody), contract_idx).shuffle(shuffle_idx).reshape(new_dims);
        multisite_mpo_squared = temp;
    }
    return multisite_mpo_squared;
}

Eigen::Tensor<ModelFinite::Scalar, 2> ModelFinite::get_multisite_ham_squared(const std::vector<size_t> &sites, const std::vector<size_t> &nbody_terms) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite hamiltonian tensor");
    if(sites == active_sites and cache.multisite_ham and nbody_terms.empty()) return cache.multisite_ham.value();
    // A multisite_ham is simply the corner of a multisite_mpo where the hamiltonian resides
    auto multisite_mpo_squared = get_multisite_mpo_squared(sites, nbody_terms);
    auto edgeL                 = get_mpo(sites.front()).get_MPO2_edge_left();
    auto edgeR                 = get_mpo(sites.back()).get_MPO2_edge_right();
    return multisite_mpo_squared.contract(edgeL, tenx::idx({0}, {0}))
        .contract(edgeR, tenx::idx({0}, {0}))
        .reshape(tenx::array2{multisite_mpo_squared.dimension(2), multisite_mpo_squared.dimension(3)});
}

const Eigen::Tensor<ModelFinite::Scalar, 4> &ModelFinite::get_multisite_mpo_squared() const {
    if(cache.multisite_mpo_squared) return cache.multisite_mpo_squared.value();
    cache.multisite_mpo_squared = get_multisite_mpo_squared(active_sites);
    return cache.multisite_mpo_squared.value();
}

const Eigen::Tensor<ModelFinite::Scalar, 2> &ModelFinite::get_multisite_ham_squared() const {
    if(cache.multisite_ham_squared and not active_sites.empty()) return cache.multisite_ham_squared.value();
    cache.multisite_ham_squared = get_multisite_ham_squared(active_sites);
    return cache.multisite_ham_squared.value();
}

void ModelFinite::clear_cache() const {
    tools::log->trace("Clearing model cache");
    cache = Cache();
}

std::vector<size_t> ModelFinite::get_active_ids() const {
    std::vector<size_t> ids;
    ids.reserve(active_sites.size());
    for(const auto &pos : active_sites) ids.emplace_back(get_mpo(pos).get_unique_id());
    return ids;
}