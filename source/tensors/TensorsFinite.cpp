
#include "TensorsFinite.h"
#include <debug/exceptions.h>
#include <config/debug.h>
#include <config/settings.h>
#include <math/num.h>
#include <tensors/edges/EdgesFinite.h>
#include <tensors/model/ModelFinite.h>
#include <tensors/site/mpo/MpoSite.h>
#include <tensors/site/mps/MpsSite.h>
#include <tensors/state/StateFinite.h>
#include <tools/common/log.h>
#include <tools/finite/env.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/multisite.h>
#include <tools/finite/ops.h>

TensorsFinite::TensorsFinite() : state(std::make_unique<StateFinite>()), model(std::make_unique<ModelFinite>()), edges(std::make_unique<EdgesFinite>()) {
    tools::log->trace("Constructing TensorsFinite");
}

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
TensorsFinite::~TensorsFinite()                     = default;            // default dtor
TensorsFinite::TensorsFinite(TensorsFinite &&other) = default;            // default move ctor
TensorsFinite &TensorsFinite::operator=(TensorsFinite &&other) = default; // default move assign

TensorsFinite::TensorsFinite(const TensorsFinite &other)
    : state(std::make_unique<StateFinite>(*other.state)), model(std::make_unique<ModelFinite>(*other.model)),
      edges(std::make_unique<EdgesFinite>(*other.edges)), active_sites(other.active_sites), measurements(other.measurements) {}
TensorsFinite &TensorsFinite::operator=(const TensorsFinite &other) {
    // check for self-assignment
    if(this != &other) {
        state        = std::make_unique<StateFinite>(*other.state);
        model        = std::make_unique<ModelFinite>(*other.model);
        edges        = std::make_unique<EdgesFinite>(*other.edges);
        active_sites = other.active_sites;
        measurements = other.measurements;
        sync_active_sites();
    }
    return *this;
}

TensorsFinite::TensorsFinite(AlgorithmType algo_type, ModelType model_type, size_t model_size, size_t position) : TensorsFinite() {
    initialize(algo_type, model_type, model_size, position);
}

void TensorsFinite::initialize(AlgorithmType algo_type, ModelType model_type, size_t model_size, size_t position) {
    tools::log->debug("Initializing tensors: algorithm [{}] | model [{}] | sites [{}] | position [{}]", enum2sv(algo_type), enum2sv(model_type), model_size,
                      position);
    state->initialize(algo_type, model_type, model_size, position);
    model->initialize(model_type, model_size);
    edges->initialize(model_size);
}

void TensorsFinite::randomize_model() {
    model->randomize();
    rebuild_mpo();
    rebuild_mpo_squared();
}

void TensorsFinite::randomize_state(StateInit state_init, std::string_view sector, long chi_lim, bool use_eigenspinors, std::optional<long> bitfield,
                                    std::optional<StateInitType> state_type, std::optional<svd::settings> svd_settings) {
    state->clear_measurements();
    if(not state_type) state_type = state->is_real() ? StateInitType::REAL : StateInitType::CPLX;
    tools::log->info("Randomizing state | norm {:.16f} | spins: {:+.16f}", tools::finite::measure::norm(*state),
                     fmt::join(tools::finite::measure::spin_components(*state), ", "));
    tools::finite::mps::randomize_state(*state, state_init, state_type.value(), sector, chi_lim, use_eigenspinors, bitfield);
    if(settings::strategy::project_initial_state and tools::finite::mps::init::is_valid_axis(sector)) {
        tools::log->info("Projecting state  | norm {:.16f} | spin components: {:+.16f}", tools::finite::measure::norm(*state),
                         fmt::join(tools::finite::measure::spin_components(*state), ", "));
        tools::finite::ops::project_to_nearest_sector(*state, sector, chi_lim, svd_settings); // Normalization happens after projection

        // Note! After running this function we should rebuild edges! However, there are usually no sites active at this point, so don't here.
    }
}

void TensorsFinite::normalize_state(long chi_lim, std::optional<svd::settings> svd_settings, NormPolicy norm_policy) {
    // Normalize if unity was lost for some reason (numerical error buildup)
    auto has_normalized = tools::finite::mps::normalize_state(*state, chi_lim, svd_settings, norm_policy);
    if(has_normalized) {
        state->clear_cache();
        clear_measurements();
        rebuild_edges();
    }
}

const Eigen::Tensor<TensorsFinite::Scalar, 3> &TensorsFinite::get_multisite_mps() const { return state->get_multisite_mps(); }

const Eigen::Tensor<TensorsFinite::Scalar, 4> &TensorsFinite::get_multisite_mpo() const { return model->get_multisite_mpo(); }

const Eigen::Tensor<TensorsFinite::Scalar, 4> &TensorsFinite::get_multisite_mpo_squared() const { return model->get_multisite_mpo_squared(); }

env_pair<const Eigen::Tensor<TensorsFinite::Scalar, 3>> TensorsFinite::get_multisite_env_ene_blk() const {
    return std::as_const(*edges).get_multisite_env_ene_blk();
}
env_pair<const Eigen::Tensor<TensorsFinite::Scalar, 3>> TensorsFinite::get_multisite_env_var_blk() const {
    return std::as_const(*edges).get_multisite_env_var_blk();
}

StateFinite TensorsFinite::get_state_projected_to_nearest_sector(std::string_view sector, std::optional<long> chi_lim,
                                                                 std::optional<svd::settings> svd_settings) const {
    auto state_projected = *state;
    tools::finite::ops::project_to_nearest_sector(state_projected, sector, chi_lim, svd_settings);
    return state_projected;
}

void TensorsFinite::project_to_nearest_sector(std::string_view sector, std::optional<long> chi_lim, std::optional<svd::settings> svd_settings) {
    tools::finite::ops::project_to_nearest_sector(*state, sector, chi_lim, svd_settings);
    rebuild_edges();
    if constexpr(settings::debug) assert_validity();
}

StateFinite TensorsFinite::get_state_with_hamiltonian_applied(std::optional<long> chi_lim, std::optional<svd::settings> svd_settings) const {
    auto state_projected = *state;
    state_projected.clear_measurements(LogPolicy::QUIET);
    state_projected.clear_cache(LogPolicy::QUIET);
    if(not chi_lim) chi_lim = state_projected.find_largest_chi();
    tools::finite::mps::normalize_state(state_projected, chi_lim.value(), svd_settings, NormPolicy::IFNEEDED);
    std::vector<Eigen::Tensor<Scalar, 4>> mpos;
    for(const auto &mpo : model->MPO) mpos.emplace_back(mpo->MPO());
    auto Ledge = model->MPO.front()->get_MPO_edge_left();
    auto Redge = model->MPO.back()->get_MPO_edge_right();
    tools::finite::ops::apply_mpos(state_projected, mpos, Ledge, Redge);
    tools::finite::mps::normalize_state(state_projected, chi_lim.value(), svd_settings, NormPolicy::ALWAYS);
    return state_projected;
}

void TensorsFinite::apply_hamiltonian_on_state(std::optional<long> chi_lim, std::optional<svd::settings> svd_settings) {
    clear_measurements();
    if(not chi_lim) chi_lim = state->find_largest_chi();
    tools::finite::mps::normalize_state(*state, chi_lim.value(), svd_settings, NormPolicy::IFNEEDED);
    std::vector<Eigen::Tensor<Scalar, 4>> mpos;
    for(const auto &mpo : model->MPO) mpos.emplace_back(mpo->MPO_reduced_view(0.0));
    auto                     chiL  = state->mps_sites.front()->get_chiL();
    auto                     chiR  = state->mps_sites.back()->get_chiR();
    auto                     mpoL  = model->MPO.front()->MPO().dimension(0);
    auto                     mpoR  = model->MPO.front()->MPO().dimension(1);
    Eigen::Tensor<Scalar, 3> Ledge = model->MPO.front()->get_MPO_edge_left().reshape(std::array<long, 3>{chiL, chiL, mpoL});
    Eigen::Tensor<Scalar, 3> Redge = model->MPO.back()->get_MPO_edge_right().reshape(std::array<long, 3>{chiR, chiR, mpoR});
    tools::finite::ops::apply_mpos(*state, mpos, Ledge, Redge);
    normalize_state(chi_lim.value(), svd_settings, NormPolicy::ALWAYS); // Has to be normalized ALWAYS, applying H|psi> ruins normalization!
}

void TensorsFinite::perturb_model_params(double coupling_ptb, double field_ptb, PerturbMode perturbMode) {
    measurements = MeasurementsTensorsFinite(); // State measurements can remain
    model->clear_cache();
    model->perturb_hamiltonian(coupling_ptb, field_ptb, perturbMode);
    model->build_mpo_squared();
    model->assert_validity();
}

struct DebugStatus {
    double                                              ene = std::numeric_limits<double>::quiet_NaN();
    double                                              red = std::numeric_limits<double>::quiet_NaN();
    double                                              var = std::numeric_limits<double>::quiet_NaN();
    std::string                                         tag;
    std::pair<std::vector<size_t>, std::vector<size_t>> env_ids;
    std::vector<size_t>                                 mpo_ids;

    [[nodiscard]] std::string msg() const {
        std::string msg;
        msg.append(fmt::format("Energy   [{:<20}] = {:>20.16f}\n", tag, ene));
        msg.append(fmt::format("Shift    [{:<20}] = {:>20.16f}\n", tag, red));
        msg.append(fmt::format("σ²H      [{:<20}] = {:<20.16f}\n", tag, var));
        return msg;
    }
    void print() const {
        tools::log->debug("Energy   [{:<20}] = {:>20.16f}", tag, ene);
        tools::log->debug("Shift    [{:<20}] = {:>20.16f}", tag, red);
        tools::log->debug("σ²H      [{:<20}] = {:>20.16f}", tag, var);
    }
};

std::optional<DebugStatus> get_status(TensorsFinite &tensors, std::string_view tag) {
    if constexpr(not settings::debug) return std::nullopt;
    tensors.model->clear_cache();
    tensors.measurements = MeasurementsTensorsFinite();
    tensors.rebuild_mpo_squared();
    tensors.rebuild_edges();
    DebugStatus deb;
    deb.ene     = tools::finite::measure::energy(tensors);
    deb.red     = tools::finite::measure::energy_reduced(tensors);
    deb.var     = tools::finite::measure::energy_variance(tensors);
    deb.tag     = tag;
    deb.env_ids = tensors.edges->get_active_ids();
    deb.mpo_ids = tensors.model->get_active_ids();
    return deb;
}

void TensorsFinite::reduce_mpo_energy(std::optional<double> energy_reduce_per_site) {
    std::vector<std::optional<DebugStatus>> debs;
    debs.emplace_back(get_status(*this, "Before shift"));

    if(not energy_reduce_per_site) energy_reduce_per_site = tools::finite::measure::energy_per_site(*this);
    measurements = MeasurementsTensorsFinite(); // Resets model-related measurements but not state measurements, which can remain
    model->clear_cache();

    tools::log->trace("Setting MPO energy shift {:.16f} (all edges should be rebuilt after this)", energy_reduce_per_site.value());
    model->set_reduced_energy_per_site(energy_reduce_per_site.value());
    model->clear_mpo_squared();
    model->assert_validity();

    debs.emplace_back(get_status(*this, "After shift"));

    if(energy_reduce_per_site.value() != 0) {
        auto &bef = debs.front();
        auto &aft = debs.back();
        if(bef and aft) {
            if(bef->red != aft->red) {
                if(bef->mpo_ids == aft->mpo_ids)
                    throw std::runtime_error(fmt::format("MPO id's are unchanged after energy reduction\n{}\n{}", bef->msg(), aft->msg()));
                if(bef->env_ids == aft->env_ids)
                    throw std::runtime_error(fmt::format("ENV id's are unchanged after energy reduction\n{}\n{}", bef->msg(), aft->msg()));
            }
        }
    }

    for(const auto &deb : debs)
        if(deb) deb->print();

    if constexpr(settings::debug) {
        auto &bef = debs.front();
        auto &aft = debs.back();
        if(bef and aft) {
            double delta_ene     = std::abs(bef->ene - aft->ene);
            double delta_var     = std::abs(bef->var - aft->var);
            double delta_ene_rel = delta_ene / std::abs(aft->ene) * 100;
            double delta_var_rel = delta_var / std::abs(aft->var) * 100;
            double critical_cancellation_max_decimals =
                std::numeric_limits<double>::digits10 - std::max(0.0, std::log10(std::pow(get_length<double>() * energy_reduce_per_site.value(), 2)));
            double critical_cancellation_error = std::pow(10, -critical_cancellation_max_decimals);
            tools::log->debug("Variance change              {:>20.16f}", delta_var);
            tools::log->debug("Variance change percent      {:>20.16f}", delta_var_rel);
            tools::log->debug("Critical cancellation decs   {:>20.16f}", critical_cancellation_max_decimals);
            tools::log->debug("Critical cancellation error  {:>20.16f}", critical_cancellation_error);
            if(delta_var > 1e3 * critical_cancellation_error) {
                tools::log->warn("Variance changed significantly after energy reduction+compression");
                if(delta_var > 1e6 * critical_cancellation_error) throw std::runtime_error("Energy reduction destroyed variance precision");
            }
            if(delta_ene_rel > 1e-8) {
                tools::log->warn("Energy changed significantly after energy reduction+compression");
                if(delta_ene_rel > 1e-6)
                    throw except::runtime_error("Energy reduction changed energy level {:.16f} -> {:.16f} (Δ/E = {:.3f} %)", bef->ene, aft->ene, delta_ene_rel);
            }
        }
    }
}

void TensorsFinite::rebuild_mpo() { model->build_mpo(); }

void TensorsFinite::rebuild_mpo_squared(std::optional<bool> compress, std::optional<svd::settings> svd_settings) {
    measurements = MeasurementsTensorsFinite(); // Resets model-related measurements but not state measurements, which can remain
    model->clear_cache();
    if(not compress) compress = settings::precision::use_compressed_mpo_squared_all;
    if(compress.value())
        model->compress_mpo_squared(svd_settings);
    else
        model->build_mpo_squared();
    model->assert_validity();
}

// Active sites
void TensorsFinite::sync_active_sites() { activate_sites(state->active_sites); }

void TensorsFinite::clear_active_sites() {
    tools::log->trace("Clearing active sites {}", active_sites);
    active_sites.clear();
    state->active_sites.clear();
    model->active_sites.clear();
    edges->active_sites.clear();
}

void TensorsFinite::activate_sites(const std::vector<size_t> &sites) {
    if(num::all_equal(sites, active_sites, state->active_sites, model->active_sites, edges->active_sites)) return;
    active_sites        = sites;
    state->active_sites = active_sites;
    model->active_sites = active_sites;
    edges->active_sites = active_sites;
    clear_cache(LogPolicy::QUIET);
    clear_measurements(LogPolicy::QUIET);
    rebuild_edges();
    if constexpr(settings::debug) assert_validity();
}

void TensorsFinite::activate_sites(long threshold, size_t max_sites, size_t min_sites) {
    activate_sites(tools::finite::multisite::generate_site_list(*state, threshold, max_sites, min_sites));
}

std::array<long, 3> TensorsFinite::active_problem_dims() const { return tools::finite::multisite::get_dimensions(*state, active_sites); }
long                TensorsFinite::active_problem_size() const { return tools::finite::multisite::get_problem_size(*state, active_sites); }
bool                TensorsFinite::is_real() const {
    if(settings::model::model_type == ModelType::ising_sdual) {
        if(not state->is_real()) tools::log->critical("state has imaginary part");
        if(not model->is_real()) tools::log->critical("model has imaginary part");
        if(not edges->is_real()) tools::log->critical("edges has imaginary part");
    }
    return state->is_real() and model->is_real() and edges->is_real();
}
bool TensorsFinite::has_nan() const {
    if(state->has_nan()) tools::log->critical("state has nan");
    if(model->has_nan()) tools::log->critical("model has nan");
    if(edges->has_nan()) tools::log->critical("edges has nan");
    return state->has_nan() or model->has_nan() or edges->has_nan();
}
void TensorsFinite::assert_validity() const {
    if(settings::model::model_type == ModelType::ising_sdual) {
        if(not state->is_real()) throw std::runtime_error("state has imaginary part");
        if(not model->is_real()) throw std::runtime_error("model has imaginary part");
        if(not edges->is_real()) throw std::runtime_error("edges has imaginary part");
    }

    state->assert_validity();
    model->assert_validity();
    edges->assert_validity();
    assert_edges();
}

bool TensorsFinite::has_center_point() const { return state->has_center_point(); }

template<typename T>
T TensorsFinite::get_length() const {
    if(not num::all_equal(state->get_length<size_t>(), model->get_length(), edges->get_length()))
        throw std::runtime_error(
            fmt::format("All lengths are not equal: state {} | model {} | edges {}", state->get_length<size_t>(), model->get_length(), edges->get_length()));
    return state->get_length<T>();
}

template double TensorsFinite::get_length<double>() const;
template size_t TensorsFinite::get_length<size_t>() const;
template long   TensorsFinite::get_length<long>() const;
template int    TensorsFinite::get_length<int>() const;

template<typename T>
T TensorsFinite::get_position() const {
    return state->get_position<T>();
}
template size_t TensorsFinite::get_position<size_t>() const;
template long   TensorsFinite::get_position<long>() const;

bool   TensorsFinite::position_is_the_middle() const { return state->position_is_the_middle(); }
bool   TensorsFinite::position_is_the_middle_any_direction() const { return state->position_is_the_middle_any_direction(); }
bool   TensorsFinite::position_is_outward_edge_left(size_t nsite) const { return state->position_is_outward_edge_left(nsite); }
bool   TensorsFinite::position_is_outward_edge_right(size_t nsite) const { return state->position_is_outward_edge_right(nsite); }
bool   TensorsFinite::position_is_outward_edge(size_t nsite) const { return state->position_is_outward_edge(nsite); }
bool   TensorsFinite::position_is_inward_edge_left(size_t nsite) const { return state->position_is_inward_edge_left(nsite); }
bool   TensorsFinite::position_is_inward_edge_right(size_t nsite) const { return state->position_is_inward_edge_right(nsite); }
bool   TensorsFinite::position_is_inward_edge(size_t nsite) const { return state->position_is_inward_edge(nsite); }
bool   TensorsFinite::position_is_at(long pos) const { return state->position_is_at(pos); }
bool   TensorsFinite::position_is_at(long pos, int dir) const { return state->position_is_at(pos, dir); }
bool   TensorsFinite::position_is_at(long pos, int dir, bool isCenter) const { return state->position_is_at(pos, dir, isCenter); }
size_t TensorsFinite::move_center_point(long chi_lim, std::optional<svd::settings> svd_settings) {
    return tools::finite::mps::move_center_point_single_site(*state, chi_lim, svd_settings);
}
size_t TensorsFinite::move_center_point_to_edge(long chi_lim, std::optional<svd::settings> svd_settings) {
    return tools::finite::mps::move_center_point_to_edge(*state, chi_lim, svd_settings);
}
size_t TensorsFinite::move_center_point_to_middle(long chi_lim, std::optional<svd::settings> svd_settings) {
    return tools::finite::mps::move_center_point_to_middle(*state, chi_lim, svd_settings);
}

void TensorsFinite::merge_multisite_mps(const Eigen::Tensor<Scalar, 3> &multisite_tensor, long chi_lim, std::optional<svd::settings> svd_settings,
                                        LogPolicy log_policy) {
    // Make sure the active sites are the same everywhere
    if(not num::all_equal(active_sites, state->active_sites, model->active_sites, edges->active_sites))
        throw std::runtime_error("All active sites are not equal: tensors {} | state {} | model {} | edges {}");
    clear_measurements(LogPolicy::QUIET);
    tools::finite::mps::merge_multisite_mps(*state, multisite_tensor, active_sites, get_position<long>(), chi_lim, svd_settings, log_policy);
    normalize_state(chi_lim, svd_settings, NormPolicy::IFNEEDED);
}

std::vector<size_t> TensorsFinite::expand_subspace(std::optional<double> alpha, long chi_lim, std::optional<svd::settings> svd_settings) {
    if(active_sites.empty()) throw std::runtime_error("No active sites for subspace expansion");
    // Follows the subspace expansion technique explained in https://link.aps.org/doi/10.1103/PhysRevB.91.155115
    auto pos_expanded = tools::finite::env::expand_subspace(*state, *model, *edges, alpha, chi_lim, svd_settings);
    if(alpha) clear_measurements(LogPolicy::QUIET); // No change if alpha == std::nullopt
    if constexpr(settings::debug) assert_validity();
    return pos_expanded;
}

void TensorsFinite::assert_edges() const { tools::finite::env::assert_edges(*state, *model, *edges); }
void TensorsFinite::assert_edges_ene() const { tools::finite::env::assert_edges_ene(*state, *model, *edges); }
void TensorsFinite::assert_edges_var() const { tools::finite::env::assert_edges_var(*state, *model, *edges); }
void TensorsFinite::rebuild_edges() { tools::finite::env::rebuild_edges(*state, *model, *edges); }
void TensorsFinite::rebuild_edges_ene() { tools::finite::env::rebuild_edges_ene(*state, *model, *edges); }
void TensorsFinite::rebuild_edges_var() { tools::finite::env::rebuild_edges_var(*state, *model, *edges); }
void TensorsFinite::do_all_measurements() const { tools::finite::measure::do_all_measurements(*this); }
void TensorsFinite::redo_all_measurements() const {
    clear_cache();
    clear_measurements();
    tools::finite::measure::do_all_measurements(*this);
}

void TensorsFinite::clear_measurements(LogPolicy logPolicy) const {
    measurements = MeasurementsTensorsFinite();
    state->clear_measurements(logPolicy);
}

void TensorsFinite::clear_cache(LogPolicy logPolicy) const {
    model->clear_cache(logPolicy);
    state->clear_cache(logPolicy);
}