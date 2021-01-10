
#include "class_tensors_finite.h"
#include <config/debug.h>
#include <config/nmspc_settings.h>
#include <iostream>
#include <math/num.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/finite/env.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/multisite.h>
#include <tools/finite/ops.h>

class_tensors_finite::class_tensors_finite()
    : state(std::make_unique<class_state_finite>()), model(std::make_unique<class_model_finite>()), edges(std::make_unique<class_edges_finite>()) {
    tools::log->trace("Constructing class_tensors_finite");
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
class_tensors_finite::~class_tensors_finite()                            = default;            // default dtor
class_tensors_finite::class_tensors_finite(class_tensors_finite &&other) = default;            // default move ctor
class_tensors_finite &class_tensors_finite::operator=(class_tensors_finite &&other) = default; // default move assign

class_tensors_finite::class_tensors_finite(const class_tensors_finite &other)
    : state(std::make_unique<class_state_finite>(*other.state)), model(std::make_unique<class_model_finite>(*other.model)),
      edges(std::make_unique<class_edges_finite>(*other.edges)), active_sites(other.active_sites), measurements(other.measurements) {}
class_tensors_finite &class_tensors_finite::operator=(const class_tensors_finite &other) {
    // check for self-assignment
    if(this != &other) {
        state        = std::make_unique<class_state_finite>(*other.state);
        model        = std::make_unique<class_model_finite>(*other.model);
        edges        = std::make_unique<class_edges_finite>(*other.edges);
        active_sites = other.active_sites;
        measurements = other.measurements;
    }
    return *this;
}

void class_tensors_finite::initialize(ModelType model_type, size_t model_size, size_t position) {
    tools::log->trace("Initializing tensors with {} sites at position {}", model_size,position);
    state->initialize(model_type, model_size, position);
    model->initialize(model_type, model_size);
    edges->initialize(model_size);
    tools::log->trace("Initializing tensors with {} sites at position {} ... OK", state->get_length(),state->get_position());
}

void class_tensors_finite::randomize_model() {
    model->randomize();
    model->rebuild_mpo_squared();
    rebuild_edges();
}

void class_tensors_finite::randomize_state(StateInit state_init,const std::string &sector, long chi_lim, bool use_eigenspinors, std::optional<long> bitfield,  std::optional<StateInitType> state_type,
                                           std::optional<double> svd_threshold) {
    state->clear_measurements();
    if(not state_type) state_type = state->is_real() ? StateInitType::REAL : StateInitType::CPLX;
    tools::log->info("Randomizing state | norm {:.16f} | spins: {:+.16f}", tools::finite::measure::norm(*state), fmt::join(tools::finite::measure::spin_components(*state), ", "));
    tools::finite::mps::randomize_state(*state, state_init, state_type.value(),sector, chi_lim, use_eigenspinors, bitfield);
    if(settings::strategy::project_initial_state){
        tools::log->info("Projecting state  | norm {:.16f} | spins: {:+.16f}", tools::finite::measure::norm(*state), fmt::join(tools::finite::measure::spin_components(*state), ", "));
        project_to_nearest_sector(sector, chi_lim, svd_threshold); // Normalization happens after projection
    }
}

void class_tensors_finite::normalize_state(long chi_lim, std::optional<double> svd_threshold, NormPolicy norm_policy) {
    // Normalize if unity was lost for some reason (numerical error buildup)
    auto has_normalized = tools::finite::mps::normalize_state(*state, chi_lim, svd_threshold, norm_policy);
    if(has_normalized) {
        state->clear_cache();
        clear_measurements();
        rebuild_edges();
        assert_validity();
    }
}


const Eigen::Tensor<class_tensors_finite::Scalar, 3> & class_tensors_finite::get_multisite_mps() const{return state->get_multisite_mps();}

const Eigen::Tensor<class_tensors_finite::Scalar, 4> & class_tensors_finite::get_multisite_mpo() const{return model->get_multisite_mpo();}

const Eigen::Tensor<class_tensors_finite::Scalar, 4> & class_tensors_finite::get_multisite_mpo_squared() const { return model->get_multisite_mpo_squared(); }

env_pair<const Eigen::Tensor<class_tensors_finite::Scalar, 3>> class_tensors_finite::get_multisite_ene_blk() const { return std::as_const(*edges).get_multisite_ene_blk(); }
env_pair<const Eigen::Tensor<class_tensors_finite::Scalar, 3>> class_tensors_finite::get_multisite_var_blk() const { return std::as_const(*edges).get_multisite_var_blk(); }


class_state_finite class_tensors_finite::get_state_projected_to_nearest_sector(const std::string &sector, std::optional<long> chi_lim, std::optional<double> svd_threshold){
    auto state_projected = *state;
    state_projected.clear_measurements();
    state_projected.clear_cache();
    if(not chi_lim) chi_lim = state_projected.find_largest_chi();
    tools::finite::mps::normalize_state(state_projected, chi_lim.value(), svd_threshold, NormPolicy::IFNEEDED);
    tools::finite::ops::project_to_nearest_sector(state_projected, sector);
    tools::finite::mps::normalize_state(state_projected, chi_lim.value(), svd_threshold, NormPolicy::ALWAYS);
    return state_projected;
}

void class_tensors_finite::project_to_nearest_sector(const std::string &sector, std::optional<long> chi_lim, std::optional<double> svd_threshold) {
    clear_measurements();
    if(not chi_lim) chi_lim = state->find_largest_chi();
    tools::finite::mps::normalize_state(*state, chi_lim.value(), svd_threshold, NormPolicy::IFNEEDED);
    tools::finite::ops::project_to_nearest_sector(*state, sector);
    normalize_state(chi_lim.value(), svd_threshold, NormPolicy::ALWAYS); // Has to be normalized ALWAYS, projection ruins normalization!
}

void class_tensors_finite::perturb_model_params(double coupling_ptb, double field_ptb, PerturbMode perturbMode) {
    measurements = tensors_measure_finite(); // State measurements can remain
    model->clear_cache();
    model->perturb_hamiltonian(coupling_ptb, field_ptb, perturbMode);
    model->rebuild_mpo_squared();
    model->assert_validity();
    rebuild_edges();
}

void class_tensors_finite::reduce_mpo_energy(std::optional<double> site_energy) {
    [[maybe_unused]] double var_bef_reduce    = 1.0;
    [[maybe_unused]] double var_aft_reduce    = 1.0;
    [[maybe_unused]] double var_aft_compr     = 1.0;
    [[maybe_unused]] double ene_bef_reduce    = 1.0;
    [[maybe_unused]] double ene_aft_reduce    = 1.0;
    [[maybe_unused]] double ene_aft_compr     = 1.0;

    if constexpr(settings::debug) {
        model->reset_mpo_squared();
        rebuild_edges();
        clear_cache();
        clear_measurements();
        tools::log->info("Measuring ene bef reduce using sites {} with ids: mps {} mpo {} env {}",
                         active_sites,
                         state->get_active_ids(),
                         model->get_active_ids(),
                         edges->get_active_ids());
        ene_bef_reduce = tools::finite::measure::energy(*this);
        var_bef_reduce = tools::finite::measure::energy_variance(*this);
    }
    if(not site_energy) site_energy = tools::finite::measure::energy_per_site(*this);
    if(site_energy.value() == model->get_energy_per_site_reduced()) return;

    measurements = tensors_measure_finite(); // Resets model-related measurements but not state measurements, which can remain
    model->clear_cache();
    if(not model->has_mpo_squared()) model->reset_mpo_squared();
    auto mpo_ids_bef = model->get_active_ids();
    auto env_ids_bef = edges->get_active_ids();
    tools::log->trace("Reducing MPO energy (all edges should be rebuilt after this)");
    model->set_reduced_energy_per_site(site_energy.value());
    if(not model->has_mpo_squared()) model->reset_mpo_squared();
    rebuild_edges();
    auto mpo_ids_aft = model->get_active_ids();
    auto env_ids_aft = edges->get_active_ids();
    if(site_energy.value() != 0 and mpo_ids_bef == mpo_ids_aft) throw std::runtime_error(fmt::format("MPO id's at sites {} are unchanged after energy reduction", active_sites));
    if(site_energy.value() != 0 and env_ids_bef == env_ids_aft) throw std::runtime_error(fmt::format("ENV id's at sites {} are unchanged after energy reduction", active_sites));


    if constexpr(settings::debug) {
        tools::log->trace("Resetting MPO energy (var edges should be rebuilt after this)");
        if(not model->has_mpo_squared()) model->reset_mpo_squared();
        rebuild_edges();
        clear_cache();
        clear_measurements();
        tools::log->info("Measuring ene aft reduce using sites {} with ids: mps {} mpo {} env {}",
                         active_sites,
                         state->get_active_ids(),
                         model->get_active_ids(),
                         edges->get_active_ids());

        ene_aft_reduce     = tools::finite::measure::energy(*this);
        var_aft_reduce     = tools::finite::measure::energy_variance(*this);
    }

    model->rebuild_mpo_squared();
    model->assert_validity();
    rebuild_edges();

    if constexpr(settings::debug) {
        clear_cache();
        clear_measurements();
    }

    if constexpr(settings::debug) {
        double energy_change                      = std::abs(ene_bef_reduce - ene_aft_compr);
        double energy_change_percent              = energy_change / std::abs(ene_bef_reduce) * 100;
        double critical_cancellation_max_decimals = 15 - std::max(0.0,std::log10(std::pow(site_energy.value()* static_cast<double>(get_length()), 2)));
        double critical_cancellation_error        = std::pow(10, -critical_cancellation_max_decimals);
        double variance_change                    = std::abs(var_bef_reduce - var_aft_compr);
        double variance_change_percent            = variance_change / std::abs(var_bef_reduce) * 100;
        tools::log->debug("Energy   before mpo reduce       {:>20.16f}", ene_bef_reduce);
        tools::log->debug("Energy   after  mpo reduce       {:>20.16f}", ene_aft_reduce);
        tools::log->debug("Energy   after  mpo compression  {:>20.16f}", ene_aft_compr);
        tools::log->debug("Variance before mpo reduce       {:>20.16f}", std::log10(var_bef_reduce));
        tools::log->debug("Variance after  mpo reduce       {:>20.16f}", std::log10(var_aft_reduce));
        tools::log->debug("Variance after  mpo compression  {:>20.16f}", std::log10(var_aft_compr));
        tools::log->debug("Variance change                  {:>20.16f}", variance_change);
        tools::log->debug("Variance change percent          {:>20.16f}", variance_change_percent);
        tools::log->debug("Critical cancellation decimals   {:>20.16f}", critical_cancellation_max_decimals);
        tools::log->debug("Critical cancellation error      {:>20.16f}", critical_cancellation_error);
        if(variance_change > 1e3 * critical_cancellation_error) {
            tools::log->warn("Variance changed significantly after energy reduction+compression");
            if(variance_change > 1e6 * critical_cancellation_error)
                throw std::runtime_error("Energy reduction destroyed variance precision");
        }
        if(energy_change_percent > 1e-8){
            tools::log->warn("Energy changed significantly after energy reduction+compression");
            if(energy_change_percent > 1e-6)
                throw std::runtime_error("Energy reduction changed energy level");
        }
    }
}

void class_tensors_finite::rebuild_mpo_squared(std::optional<SVDMode> svdMode) {
    measurements = tensors_measure_finite(); // Resets model-related measurements but not state measurements, which can remain
    model->clear_cache();
    model->rebuild_mpo_squared(svdMode);
    model->assert_validity();
    rebuild_edges();
}


void class_tensors_finite::damp_model_disorder(double coupling_damp, double field_damp) {
    measurements = tensors_measure_finite(); // Resets model-related measurements but not state measurements, which can remain
    model->clear_cache();
    model->damp_model_disorder(coupling_damp, field_damp);
    model->rebuild_mpo_squared();
    model->assert_validity();
    rebuild_edges();
}

// Active sites
void class_tensors_finite::sync_active_sites() {
    active_sites        = state->active_sites;
    model->active_sites = state->active_sites;
    edges->active_sites = state->active_sites;
    rebuild_edges();
}

void class_tensors_finite::activate_sites(long threshold, size_t max_sites, size_t min_sites) {
    active_sites        = tools::finite::multisite::generate_site_list(*state, threshold, max_sites, min_sites);
    state->active_sites = active_sites;
    model->active_sites = active_sites;
    edges->active_sites = active_sites;
    rebuild_edges();
    clear_cache();
    clear_measurements();
}

void class_tensors_finite::activate_sites(const std::vector<size_t> &sites) {
    if(active_sites == sites) {
        state->active_sites = active_sites;
        model->active_sites = active_sites;
        edges->active_sites = active_sites;
        return;
    }

    active_sites        = sites;
    state->active_sites = active_sites;
    model->active_sites = active_sites;
    edges->active_sites = active_sites;
    rebuild_edges();
    clear_cache();
    clear_measurements();
}

void class_tensors_finite::activate_truncated_sites(long threshold, long chi_lim, size_t max_sites, size_t min_sites) {
    active_sites        = tools::finite::multisite::generate_truncated_site_list(*state, threshold, chi_lim, max_sites, min_sites);
    state->active_sites = active_sites;
    model->active_sites = active_sites;
    edges->active_sites = active_sites;
    rebuild_edges(); // will only produce missing edges
    clear_cache();   // clear multisite tensors and mpos which become invalid
    clear_measurements();
    // TODO: Should measurements be cleared here? Probably not
    //    clear_measurements();
}

Eigen::DSizes<long, 3> class_tensors_finite::active_problem_dims() const { return tools::finite::multisite::get_dimensions(*state, active_sites); }
long class_tensors_finite::active_problem_size() const { return tools::finite::multisite::get_problem_size(*state, active_sites); }
bool class_tensors_finite::is_real() const { return state->is_real() and model->is_real() and edges->is_real(); }
bool class_tensors_finite::has_nan() const { return state->has_nan() or model->has_nan() or edges->has_nan(); }
void class_tensors_finite::assert_validity() const {
    state->assert_validity();
    model->assert_validity();
    edges->assert_validity();
}

bool class_tensors_finite::has_center_point() const { return state->has_center_point();}

template<typename T>
T class_tensors_finite::get_length() const {
    if(not num::all_equal(state->get_length<size_t>(), model->get_length(), edges->get_length()))
        throw std::runtime_error(fmt::format("All lengths are not equal: state {} | model {} | edges {}",
                                             state->get_length<size_t>(), model->get_length(), edges->get_length()));
    return state->get_length<T>();
}

template size_t class_tensors_finite::get_length<size_t>() const;
template long class_tensors_finite::get_length<long>() const;

template<typename T>
T class_tensors_finite::get_position() const { return state->get_position<T>(); }
template size_t class_tensors_finite::get_position<size_t>() const;
template long class_tensors_finite::get_position<long>() const;

bool class_tensors_finite::position_is_the_middle() const { return state->position_is_the_middle(); }
bool class_tensors_finite::position_is_the_middle_any_direction() const { return state->position_is_the_middle_any_direction(); }
bool class_tensors_finite::position_is_outward_edge_left(size_t nsite) const { return state->position_is_outward_edge_left(nsite); }
bool class_tensors_finite::position_is_outward_edge_right(size_t nsite) const { return state->position_is_outward_edge_right(nsite); }
bool class_tensors_finite::position_is_outward_edge(size_t nsite) const { return state->position_is_outward_edge(nsite); }
bool class_tensors_finite::position_is_inward_edge_left(size_t nsite) const { return state->position_is_inward_edge_left(nsite); }
bool class_tensors_finite::position_is_inward_edge_right(size_t nsite) const { return state->position_is_inward_edge_right(nsite); }
bool class_tensors_finite::position_is_inward_edge(size_t nsite) const { return state->position_is_inward_edge(nsite); }
bool class_tensors_finite::position_is_at(long pos) const { return state->position_is_at(pos); }
bool class_tensors_finite::position_is_at(long pos, int dir) const { return state->position_is_at(pos, dir); }
bool class_tensors_finite::position_is_at(long pos, int dir, bool isCenter) const { return state->position_is_at(pos, dir, isCenter); }
size_t class_tensors_finite::move_center_point(long chi_lim, std::optional<double> svd_threshold) {
    return tools::finite::mps::move_center_point_single_site(*state, chi_lim, svd_threshold);
}
size_t class_tensors_finite::move_center_point_to_edge(long chi_lim, std::optional<double> svd_threshold) {
    return tools::finite::mps::move_center_point_to_edge(*state, chi_lim, svd_threshold);
}
size_t class_tensors_finite::move_center_point_to_middle(long chi_lim, std::optional<double> svd_threshold) {
    return tools::finite::mps::move_center_point_to_middle(*state, chi_lim, svd_threshold);
}

void class_tensors_finite::merge_multisite_tensor(const Eigen::Tensor<Scalar, 3> &multisite_tensor, long chi_lim, std::optional<double> svd_threshold, LogPolicy log_policy) {
    // Make sure the active sites are the same everywhere
    if(not num::all_equal(active_sites, state->active_sites, model->active_sites, edges->active_sites))
        throw std::runtime_error("All active sites are not equal: tensors {} | state {} | model {} | edges {}");
    clear_measurements();
    tools::finite::mps::merge_multisite_tensor(*state, multisite_tensor, active_sites, get_position<long>(), chi_lim, svd_threshold,log_policy);
    normalize_state(chi_lim, svd_threshold, NormPolicy::IFNEEDED);
    rebuild_edges();
}

void class_tensors_finite::assert_edges() const { tools::finite::env::assert_edges(*state,*model, *edges); }
void class_tensors_finite::assert_edges_ene() const { tools::finite::env::assert_edges_ene(*state,*model, *edges); }
void class_tensors_finite::assert_edges_var() const { tools::finite::env::assert_edges_var(*state,*model, *edges); }
void class_tensors_finite::rebuild_edges() { tools::finite::env::rebuild_edges(*state, *model, *edges); }
void class_tensors_finite::rebuild_edges_ene() { tools::finite::env::rebuild_edges_ene(*state, *model, *edges); }
void class_tensors_finite::rebuild_edges_var() { tools::finite::env::rebuild_edges_var(*state, *model, *edges); }
void class_tensors_finite::do_all_measurements() const { tools::finite::measure::do_all_measurements(*this); }
void class_tensors_finite::clear_measurements() const {
    measurements = tensors_measure_finite();
    state->clear_measurements();
}

void class_tensors_finite::clear_cache() const {
    model->clear_cache();
    state->clear_cache();
}