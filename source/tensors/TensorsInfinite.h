#pragma once
#include "config/enums.h"
#include "math/svd/config.h"
#include "math/tenx/fwd_decl.h"
#include "measure/MeasurementsTensorsInfinite.h"
#include <complex>
#include <memory>

class StateInfinite;
class ModelInfinite;
class EdgesInfinite;

class TensorsInfinite {
    public:
    using Scalar = std::complex<double>;
    std::unique_ptr<StateInfinite>      state;
    std::unique_ptr<ModelInfinite>      model;
    std::unique_ptr<EdgesInfinite>      edges;
    mutable MeasurementsTensorsInfinite measurements;

    // This class should have these responsibilities:
    //  - Initialize the tensors
    //  - Move/manage center position
    //  - Rebuild edges
    //  - Activate sites
    //  - Manage measurements cache

    TensorsInfinite();
    ~TensorsInfinite();                                       // Read comment on implementation
    TensorsInfinite(TensorsInfinite &&other);                 // default move ctor
    TensorsInfinite &operator=(TensorsInfinite &&other);      // default move assign
    TensorsInfinite(const TensorsInfinite &other);            // copy ctor
    TensorsInfinite &operator=(const TensorsInfinite &other); // copy assign

    void initialize(ModelType model_type);
    void initialize_model();
    void assert_validity() const;

    [[nodiscard]] size_t get_length() const;
    [[nodiscard]] size_t get_position() const;
    [[nodiscard]] bool   is_real() const;
    [[nodiscard]] bool   has_nan() const;

    /* clang-format off */
    void reset_to_random_product_state(std::string_view  sector, bool use_eigenspinors, std::string & pattern);
    /* clang-format on */

    void reset_edges();
    void eject_edges();

    void merge_twosite_tensor(const Eigen::Tensor<Scalar, 3> &twosite_tensor, std::optional<svd::config> svd_cfg = std::nullopt);
    void enlarge();
    void clear_measurements() const;
    void clear_cache() const;
};
