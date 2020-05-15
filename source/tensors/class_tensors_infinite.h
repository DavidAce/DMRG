#pragma once
#include <config/enums.h>
#include <list>
#include <measure/tensors_measure_infinite.h>
#include <memory>
class class_state_infinite;
class class_model_infinite;
class class_edges_infinite;

class class_tensors_infinite {
    public:
    std::unique_ptr<class_state_infinite> state;
    std::unique_ptr<class_model_infinite> model;
    std::unique_ptr<class_edges_infinite> edges;
    // This class should have these responsibilities:
    //  - Initialize the tensors
    //  - Move/manage center position
    //  - Rebuild edges
    //  - Activate sites
    //  - Manage measurements cache

    class_tensors_infinite();
    void initialize(ModelType model_type);

    [[nodiscard]] size_t get_length() const;
    [[nodiscard]] size_t get_position() const;
    [[nodiscard]] bool   is_real() const;
    [[nodiscard]] bool   has_nan() const;
    void                 assert_validity() const;
    void                 insert_site_pair();
    void                 do_all_measurements() const;
    void                 clear_measurements() const;
    void                 clear_cache() const;

    mutable tensors_measure_infinite measurements;
};
