#pragma once
#include <complex>
#include <config/enums.h>
#include <measure/tensors_measure_infinite.h>
#include <memory>
#include <general/eigen_tensor_fwd_decl.h>

class class_state_infinite;
class class_model_infinite;
class class_edges_infinite;

class class_tensors_infinite {
    public:
    using Scalar = std::complex<double>;
    std::unique_ptr<class_state_infinite> state;
    std::unique_ptr<class_model_infinite> model;
    std::unique_ptr<class_edges_infinite> edges;
    mutable tensors_measure_infinite      measurements;

    // This class should have these responsibilities:
    //  - Initialize the tensors
    //  - Move/manage center position
    //  - Rebuild edges
    //  - Activate sites
    //  - Manage measurements cache

    class_tensors_infinite();
    ~class_tensors_infinite();                                                  // Read comment on implementation
    class_tensors_infinite(class_tensors_infinite &&other);                     // default move ctor
    class_tensors_infinite &operator=(class_tensors_infinite &&other);          // default move assign
    class_tensors_infinite(const class_tensors_infinite &other);                // copy ctor
    class_tensors_infinite &operator=(const class_tensors_infinite &other);     // copy assign

    void                 initialize(ModelType model_type);
    void                 assert_validity() const;
    [[nodiscard]] size_t get_length() const;
    [[nodiscard]] size_t get_position() const;
    [[nodiscard]] bool   is_real() const;
    [[nodiscard]] bool   has_nan() const;

    /* clang-format off */
    void reset_to_random_product_state(const std::string & sector, long bitfield, bool use_eigenspinors);
    /* clang-format on */

    void reset_edges();
    void eject_edges();


    void merge_multisite_tensor(const Eigen::Tensor<Scalar, 3> &twosite_tensor);
    void enlarge();
    void do_all_measurements() const;
    void clear_measurements() const;
    void clear_cache() const;
};
