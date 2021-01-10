#pragma once
#include <complex>
#include <config/enums.h>
#include <general/eigen_tensor_fwd_decl.h>
#include <measure/tensors_measure_finite.h>
#include <memory>
#include <tensors/edges/class_env_pair.h>

class class_state_finite;
class class_model_finite;
class class_edges_finite;

class class_tensors_finite {
    public:
    using Scalar = std::complex<double>;
    std::unique_ptr<class_state_finite> state;
    std::unique_ptr<class_model_finite> model;
    std::unique_ptr<class_edges_finite> edges;

    std::vector<size_t>            active_sites;
    mutable tensors_measure_finite measurements;

    // This class should have these responsibilities:
    //  - Initialize/randomize the tensors
    //  - Move/manage center position
    //  - Rebuild edges
    //  - Activate sites
    //  - Manage caches

    class_tensors_finite();
    ~class_tensors_finite();                                            // Read comment on implementation
    class_tensors_finite(class_tensors_finite &&other);                 // default move ctor
    class_tensors_finite &operator=(class_tensors_finite &&other);      // default move assign
    class_tensors_finite(const class_tensors_finite &other);            // copy ctor
    class_tensors_finite &operator=(const class_tensors_finite &other); // copy assign

    void initialize(ModelType model_type, size_t model_size, size_t position);
    void randomize_model();
    void randomize_state(StateInit state_init, const std::string &sector, long chi_lim, bool use_eigenspinors, std::optional<long> bitfield = std::nullopt,
                         std::optional<StateInitType> state_type = std::nullopt, std::optional<double> svd_threshold = std::nullopt);
    void normalize_state(long chi_lim, std::optional<double> svd_threshold = std::nullopt, NormPolicy policy = NormPolicy::IFNEEDED);

    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &         get_multisite_mps() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 4> &         get_multisite_mpo() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 4> &         get_multisite_mpo_squared() const;
    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 3>> get_multisite_ene_blk() const;
    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 3>> get_multisite_var_blk() const;

    [[nodiscard]] class_state_finite get_state_projected_to_nearest_sector(const std::string &sector, std::optional<long> chi_lim = std::nullopt,
                                                                           std::optional<double> svd_threshold = std::nullopt);
    void project_to_nearest_sector(const std::string &sector, std::optional<long> chi_lim = std::nullopt, std::optional<double> svd_threshold = std::nullopt);
    void perturb_model_params(double coupling_ptb, double field_ptb, PerturbMode perturbMode);
    void damp_model_disorder(double coupling_damp, double field_damp);
    void reduce_mpo_energy(std::optional<double> site_energy = std::nullopt);
    void rebuild_mpo_squared(std::optional<SVDMode> svdMode = std::nullopt);

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

    void                   sync_active_sites();
    void                   activate_sites(long threshold, size_t max_sites, size_t min_sites = 1);
    void                   activate_sites(const std::vector<size_t> &sites);
    void                   activate_truncated_sites(long threshold, long chi_lim, size_t max_sites, size_t min_sites = 1);
    Eigen::DSizes<long, 3> active_problem_dims() const;
    long                   active_problem_size() const;
    void                   do_all_measurements() const;
    size_t                 move_center_point(long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    size_t                 move_center_point_to_edge(long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    size_t                 move_center_point_to_middle(long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    void merge_multisite_tensor(const Eigen::Tensor<Scalar, 3> &multisite_tensor, long chi_lim, std::optional<double> svd_threshold = std::nullopt,
                                LogPolicy log_policy = LogPolicy::QUIET);

    void assert_edges() const;
    void assert_edges_ene() const;
    void assert_edges_var() const;
    void rebuild_edges();
    void rebuild_edges_ene();
    void rebuild_edges_var();
    void clear_measurements() const;
    void clear_cache() const;
};
