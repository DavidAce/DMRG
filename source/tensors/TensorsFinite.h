#pragma once

#include <array>
#include <complex>
#include <config/enums.h>
#include <general/eigen_tensor_fwd_decl.h>
#include <math/svd/settings.h>
#include <measure/tensors_measure_finite.h>
#include <memory>
#include <tensors/site/env/EnvPair.h>
class StateFinite;
class ModelFinite;
class EdgesFinite;

class TensorsFinite {
    public:
    using Scalar = std::complex<double>;
    std::unique_ptr<StateFinite> state;
    std::unique_ptr<ModelFinite> model;
    std::unique_ptr<EdgesFinite> edges;

    std::vector<size_t>            active_sites;
    mutable tensors_measure_finite measurements;

    // This class should have these responsibilities:
    //  - Initialize/randomize the tensors
    //  - Move/manage center position
    //  - Rebuild edges
    //  - Activate sites
    //  - Manage caches

    TensorsFinite();
    ~TensorsFinite();                                     // Read comment on implementation
    TensorsFinite(TensorsFinite &&other);                 // default move ctor
    TensorsFinite &operator=(TensorsFinite &&other);      // default move assign
    TensorsFinite(const TensorsFinite &other);            // copy ctor
    TensorsFinite &operator=(const TensorsFinite &other); // copy assign

    void initialize(ModelType model_type, size_t model_size, size_t position);
    void randomize_model();
    void randomize_state(StateInit state_init, std::string_view sector, long chi_lim, bool use_eigenspinors, std::optional<long> bitfield = std::nullopt,
                         std::optional<StateInitType> state_type = std::nullopt, std::optional<svd::settings> svd_settings = std::nullopt);
    void normalize_state(long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt, NormPolicy policy = NormPolicy::IFNEEDED);

    [[nodiscard]] const Eigen::Tensor<Scalar, 3>          &get_multisite_mps() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 4>          &get_multisite_mpo() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 4>          &get_multisite_mpo_squared() const;
    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 3>> get_multisite_ene_blk() const;
    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 3>> get_multisite_var_blk() const;

    [[nodiscard]] StateFinite get_state_projected_to_nearest_sector(std::string_view sector, std::optional<long> chi_lim = std::nullopt,
                                                                    std::optional<svd::settings> svd_settings = std::nullopt);
    void                      project_to_nearest_sector(std::string_view sector, std::optional<long> chi_lim = std::nullopt,
                                                        std::optional<svd::settings> svd_settings = std::nullopt);
    [[nodiscard]] StateFinite get_state_with_hamiltonian_applied(std::optional<long>          chi_lim      = std::nullopt,
                                                                 std::optional<svd::settings> svd_settings = std::nullopt);
    void                      apply_hamiltonian_on_state(std::optional<long> chi_lim = std::nullopt, std::optional<svd::settings> svd_settings = std::nullopt);

    void perturb_model_params(double coupling_ptb, double field_ptb, PerturbMode perturbMode);
    void damp_model_disorder(double coupling_damp, double field_damp);
    void reduce_mpo_energy(std::optional<double> energy_reduce_per_site = std::nullopt);
    void rebuild_mpo_squared(std::optional<bool> compress = std::nullopt, std::optional<svd::settings> svd_settings = std::nullopt);
    void compress_mpo_squared(std::optional<svd::settings> svd_settings = std::nullopt);

    void assert_validity() const;

    template<typename T = size_t>
    [[nodiscard]] T get_position() const;
    template<typename T = size_t>
    [[nodiscard]] T get_length() const;

    [[nodiscard]] bool is_real() const;
    [[nodiscard]] bool has_nan() const;
    [[nodiscard]] bool has_center_point() const;
    [[nodiscard]] bool position_is_the_middle() const;
    [[nodiscard]] bool position_is_the_middle_any_direction() const;
    [[nodiscard]] bool position_is_outward_edge_left(size_t nsite = 1) const;
    [[nodiscard]] bool position_is_outward_edge_right(size_t nsite = 1) const;
    [[nodiscard]] bool position_is_outward_edge(size_t nsite = 1) const;
    [[nodiscard]] bool position_is_inward_edge_left(size_t nsite = 1) const;
    [[nodiscard]] bool position_is_inward_edge_right(size_t nsite = 1) const;
    [[nodiscard]] bool position_is_inward_edge(size_t nsite = 1) const;
    [[nodiscard]] bool position_is_at(long pos) const;
    [[nodiscard]] bool position_is_at(long pos, int dir) const;
    [[nodiscard]] bool position_is_at(long pos, int dir, bool isCenter) const;

    void                sync_active_sites();
    void                activate_sites(long threshold, size_t max_sites, size_t min_sites = 1);
    void                activate_sites(const std::vector<size_t> &sites);
    void                activate_truncated_sites(long threshold, long chi_lim, size_t max_sites, size_t min_sites = 1);
    std::array<long, 3> active_problem_dims() const;
    long                active_problem_size() const;
    void                do_all_measurements() const;
    size_t              move_center_point(long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    size_t              move_center_point_to_edge(long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    size_t              move_center_point_to_middle(long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    void merge_multisite_tensor(const Eigen::Tensor<Scalar, 3> &multisite_tensor, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt,
                                LogPolicy log_policy = LogPolicy::QUIET);

    std::vector<size_t> expand_subspace(std::optional<double> alpha, long chi_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    void                assert_edges() const;
    void                assert_edges_ene() const;
    void                assert_edges_var() const;
    void                rebuild_edges();
    void                rebuild_edges_ene();
    void                rebuild_edges_var();
    void                clear_measurements() const;
    void                clear_cache() const;
};
