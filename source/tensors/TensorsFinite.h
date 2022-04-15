#pragma once

#include <array>
#include <complex>
#include <config/enums.h>
#include <math/svd/settings.h>
#include <math/tenx/fwd_decl.h>
#include <measure/MeasurementsTensorsFinite.h>
#include <memory>
#include <tensors/site/env/EnvPair.h>
class StateFinite;
class ModelFinite;
class EdgesFinite;

class TensorsFinite {
    private:
    using cplx = std::complex<double>;
    using real = double;
    struct Cache {
        std::optional<std::vector<size_t>>    cached_sites_hamiltonian           = std::nullopt;
        std::optional<std::vector<size_t>>    cached_sites_hamiltonian_squared   = std::nullopt;
        std::optional<Eigen::Tensor<cplx, 2>> effective_hamiltonian_cplx         = std::nullopt;
        std::optional<Eigen::Tensor<cplx, 2>> effective_hamiltonian_squared_cplx = std::nullopt;
        std::optional<Eigen::Tensor<real, 2>> effective_hamiltonian_real         = std::nullopt;
        std::optional<Eigen::Tensor<real, 2>> effective_hamiltonian_squared_real = std::nullopt;
    };
    mutable Cache cache;

    public:
    std::unique_ptr<StateFinite> state;
    std::unique_ptr<ModelFinite> model;
    std::unique_ptr<EdgesFinite> edges;

    std::vector<size_t>               active_sites;
    mutable MeasurementsTensorsFinite measurements;

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
    TensorsFinite(AlgorithmType algo_type, ModelType model_type, size_t model_size, size_t position);

    void initialize(AlgorithmType algo_type, ModelType model_type, size_t model_size, size_t position);
    void randomize_model();
    void randomize_state(StateInit state_init, std::string_view sector, long bond_lim, bool use_eigenspinors, std::optional<long> bitfield = std::nullopt,
                         std::optional<StateInitType> state_type = std::nullopt);
    void normalize_state(long bond_lim, std::optional<svd::settings> svd_settings = std::nullopt, NormPolicy policy = NormPolicy::IFNEEDED);

    [[nodiscard]] const Eigen::Tensor<cplx, 3>          &get_multisite_mps() const;
    [[nodiscard]] const Eigen::Tensor<cplx, 4>          &get_multisite_mpo() const;
    [[nodiscard]] const Eigen::Tensor<cplx, 4>          &get_multisite_mpo_squared() const;
    [[nodiscard]] env_pair<const Eigen::Tensor<cplx, 3>> get_multisite_env_ene_blk() const;
    [[nodiscard]] env_pair<const Eigen::Tensor<cplx, 3>> get_multisite_env_var_blk() const;
    /* clang-format off */
    template<typename Scalar> [[nodiscard]] const Eigen::Tensor<Scalar, 2> &get_effective_hamiltonian() const;
    template<typename Scalar> [[nodiscard]] const Eigen::Tensor<Scalar, 2> &get_effective_hamiltonian_squared() const;
    /* clang-format on */

    void project_to_nearest_sector(std::string_view sector, std::optional<long> bond_lim, std::optional<bool> use_mpo2_proj = std::nullopt,
                                   std::optional<svd::settings> svd_settings = std::nullopt);
    void shift_mpo_energy(std::optional<double> energy_shift_per_site = std::nullopt);
    void set_psfactor(double psfactor);
    void rebuild_mpo();
    void rebuild_mpo_squared(std::optional<bool> compress = std::nullopt, std::optional<svd::settings> svd_settings = std::nullopt);

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
    void                clear_active_sites();
    void                activate_sites(long threshold, size_t max_sites, size_t min_sites = 1);
    void                activate_sites(const std::vector<size_t> &sites);
    std::array<long, 3> active_problem_dims() const;
    long                active_problem_size() const;
    void                do_all_measurements() const;
    void                redo_all_measurements() const;
    size_t              move_center_point(long bond_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    size_t              move_center_point_to_edge(long bond_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    size_t              move_center_point_to_middle(long bond_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    void merge_multisite_mps(const Eigen::Tensor<cplx, 3> &multisite_tensor, long bond_lim, std::optional<svd::settings> svd_settings = std::nullopt,
                             LogPolicy log_policy = LogPolicy::QUIET);

    std::vector<size_t> expand_environment(std::optional<double> alpha, long bond_lim, std::optional<svd::settings> svd_settings = std::nullopt);
    void                assert_edges() const;
    void                assert_edges_ene() const;
    void                assert_edges_var() const;
    void                rebuild_edges();
    void                rebuild_edges_ene();
    void                rebuild_edges_var();
    void                clear_measurements(LogPolicy logPolicy = LogPolicy::QUIET) const;
    void                clear_cache(LogPolicy logPolicy = LogPolicy::QUIET) const;
};
