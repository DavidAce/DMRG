#pragma once
#include <config/enums.h>
#include <list>
#include <measure/tensors_measure_finite.h>
#include <memory>
#include <complex>
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

    std::list<size_t>              active_sites;
    mutable tensors_measure_finite measurements;

    // This class should have these responsibilities:
    //  - Initialize the tensors
    //  - Move/manage center position
    //  - Rebuild edges
    //  - Activate sites
    //  - Manage measurements cache

    class_tensors_finite();
    ~class_tensors_finite();                                                // Read comment on implementation
    class_tensors_finite(class_tensors_finite &&other) noexcept;            // default move ctor
    class_tensors_finite &operator=(class_tensors_finite &&other) noexcept; // default move assign
    class_tensors_finite(const class_tensors_finite &other);                // copy ctor
    class_tensors_finite &operator=(const class_tensors_finite &other);     // copy assign

    void initialize(ModelType model_type, size_t model_size, size_t position);

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
    // Active sites
    void sync_active_sites();
    void activate_sites(long threshold, size_t max_sites, size_t min_sites = 2);
    void do_all_measurements() const;
    void move_center_point();
    void merge_multisite_tensor(const Eigen::Tensor<Scalar,3> & multisite_tensor);
    void rebuild_edges();
    void assert_validity() const;

    void clear_measurements() const;
    void clear_cache() const;
};
