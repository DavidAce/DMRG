#pragma once
#include <complex>
#include <config/enums.h>
#include <list>
#include <measure/tensors_measure_finite.h>
#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>
class class_state_finite;
class class_model_finite;
class class_edges_finite;

class class_tensors_finite {
    public:
    using Scalar = std::complex<double>;
    std::unique_ptr<class_state_finite> state;
    std::unique_ptr<class_model_finite> model;
    std::unique_ptr<class_edges_finite> edges;

    std::vector<size_t>              active_sites;
    mutable tensors_measure_finite measurements;

    // This class should have these responsibilities:
    //  - Initialize/randomize the tensors
    //  - Move/manage center position
    //  - Rebuild edges
    //  - Activate sites
    //  - Manage caches

    class_tensors_finite();
    ~class_tensors_finite();                                                // Read comment on implementation
    class_tensors_finite(class_tensors_finite &&other);                     // default move ctor
    class_tensors_finite &operator=(class_tensors_finite &&other);          // default move assign
    class_tensors_finite(const class_tensors_finite &other);                // copy ctor
    class_tensors_finite &operator=(const class_tensors_finite &other);     // copy assign

    void initialize(ModelType model_type, size_t model_size, size_t position);
    void randomize_model();
    void randomize_state(StateType state_type, const std::string &sector, long chi_lim, bool use_eigenspinors, std::optional<long> bitfield = std::nullopt, std::optional<double> svd_threshold = std::nullopt);
//    void randomize_from_current_state(const std::vector<std::string> &pauli_strings, const std::string &sector, long chi_lim, std::optional<double> svd_threshold = std::nullopt);
    void normalize_state(long chi_lim, std::optional<double> svd_threshold = std::nullopt, NormPolicy policy = NormPolicy::IFNEEDED);
//    void randomize_state(const std::string &sector, long bitfield, bool use_eigenspinors);
    void project_to_nearest_sector(const std::string & sector, std::optional<long> chi_lim = std::nullopt, std::optional<double> svd_threshold = std::nullopt);
    void perturb_model_params(double coupling_ptb, double field_ptb, PerturbMode perturbMode);
    void damp_model_disorder(double coupling_damp, double field_damp);
    void reduce_mpo_energy(std::optional<double> site_energy = std::nullopt);
    void rebuild_mpo_squared(std::optional<SVDMode> svdMode = std::nullopt);

    void                 assert_validity() const;
    [[nodiscard]] bool   is_real() const;
    [[nodiscard]] bool   has_nan() const;
    [[nodiscard]] size_t get_length() const;
    [[nodiscard]] size_t get_position() const;
    [[nodiscard]] bool   position_is_the_middle() const;
    [[nodiscard]] bool   position_is_the_middle_any_direction() const;
    [[nodiscard]] bool   position_is_left_edge() const;
    [[nodiscard]] bool   position_is_right_edge() const;
    [[nodiscard]] bool   position_is_any_edge() const;
    [[nodiscard]] bool   position_is_at(size_t pos) const;

    void sync_active_sites();
    void activate_sites(long threshold, size_t max_sites, size_t min_sites = 2);
    void activate_sites(const std::vector<size_t> & sites);
    void activate_truncated_sites(long threshold, long chi_lim, size_t max_sites, size_t min_sites = 2);
    long active_problem_size() const;
    void do_all_measurements() const;
    void move_center_point(long chi_lim,std::optional<double> svd_threshold=std::nullopt);
    void merge_multisite_tensor(const Eigen::Tensor<Scalar, 3> &multisite_tensor, long chi_lim,std::optional<double> svd_threshold = std::nullopt);

    void rebuild_edges();
    void eject_all_edges();
    void eject_inactive_edges();

    void clear_measurements() const;
    void clear_cache() const;
};
