#pragma once
#include <config/enums.h>
#include <list>
#include <measure/tensors_measure_finite.h>
#include <memory>
class class_state_finite;
class class_model_finite;
class class_edges_finite;

class class_tensors_finite {
    public:
    std::unique_ptr<class_state_finite> state;
    std::unique_ptr<class_model_finite> model;
    std::unique_ptr<class_edges_finite> edges;
    // This class should have these responsibilities:
    //  - Initialize the tensors
    //  - Move/manage center position
    //  - Rebuild edges
    //  - Activate sites
    //  - Manage measurements cache

    class_tensors_finite();
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
    std::list<size_t> active_sites;
    void              sync_active_sites();
    void              activate_sites(const size_t threshold, const size_t max_sites, const size_t min_sites = 2);
    void              do_all_measurements() const;
    void              move_center_point();
    void              rebuild_edges();
    void              assert_validity() const;

    void clear_measurements() const;
    void clear_cache() const;

    mutable tensors_measure_finite measurements;
};
