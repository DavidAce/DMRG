//
// Created by david on 2019-01-29.
//

#pragma once

#include <Eigen/Core>
#include <complex>
#include <config/enums.h>
#include <list>
#include <measure/state_measure_finite.h>
#include <memory>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>

/**
 * \class class_finite_state
 *
 * \brief Stores the finite mps state and components to operate on it, such as environments and model mpo's.
 *
 * The finite state mps is always partitioned into two lists, corresponding to the two sides left and right of the current position.
 * Each site contains a "Vidal" site, i.e., a bond matrix \f$ \Lambda \f$ and an mps tensor \f$\Gamma \f$ in Vidal's notation.
 * In the left side the sites are left normalized, with the bond matrix on the left \f$ \Lambda \Gamma \ = A\f$.
 * In the right side the sites are right normalized, with the bond matrix on the right \f$ \Gamma \Lambda = B\f$.
 * At the center is the center bond matrix MPS_C, a.k.a. \f$ \Lambda^C \f$. For a state with 10 sites the layout is seen below
 *
 * \code
 * |Finite state> =  ---MPS_L(0)--- ... MPS_L(4)---MPS_L(5)---MPS_C---MPS_R(6)---MPS_R(7)---...---MPS_R(9)
 * \endcode
 *
 *  The numbers in parentheses denote the position in the chain, note that this isn't the same as the position in the corresponding containers.
 */
class class_mps_site;
class class_tensors_finite;

class class_state_finite {
    public:
    using Scalar = std::complex<double>;

    private:
    struct Cache {
        std::optional<Eigen::Tensor<Scalar, 3>> multisite_mps = std::nullopt;
    };

    int                 direction = 1;
    mutable Cache             cache;
    mutable std::vector<bool> tag_normalized_sites;

    public:
    std::list<std::unique_ptr<class_mps_site>> mps_sites;
    std::vector<size_t>                        active_sites;
    mutable state_measure_finite               measurements;

    public:
    class_state_finite();
    ~class_state_finite();                                          // Read comment on implementation
    class_state_finite(class_state_finite &&other);                 // default move ctor
    class_state_finite &operator=(class_state_finite &&other);      // default move assign
    class_state_finite(const class_state_finite &other);            // copy ctor
    class_state_finite &operator=(const class_state_finite &other); // copy assign

    void initialize(ModelType modeltype, size_t model_size, size_t position = 0);

    const Eigen::Tensor<Scalar, 1> &midchain_bond() const;
    const Eigen::Tensor<Scalar, 1> &current_bond() const;

    long find_largest_chi() const;

    void                   set_positions();
    size_t                 get_length() const;
    size_t                 get_position() const;
    void                   flip_direction();
    int                    get_direction() const;
    std::vector<std::string> get_labels() const;
    Eigen::DSizes<long, 3> dimensions_2site() const;
    long                   size_2site() const;

    bool position_is_the_middle() const;
    bool position_is_the_middle_any_direction() const;
    bool position_is_left_edge(size_t nsite = 2) const;
    bool position_is_right_edge(size_t nsite = 2) const;
    bool position_is_any_edge(size_t nsite = 2) const;
    bool position_is_at(size_t pos) const;
    bool is_real() const;
    bool has_nan() const;
    void assert_validity() const;

    const class_mps_site &get_mps_site(size_t pos) const;
    class_mps_site &      get_mps_site(size_t pos);
    const class_mps_site &get_mps_site() const;
    class_mps_site &      get_mps_site();

    // For multisite
    Eigen::DSizes<long, 3> active_dimensions() const;
    long                   active_problem_size() const;

    Eigen::Tensor<Scalar, 3>        get_multisite_mps(const std::vector<size_t> &sites) const;
    const Eigen::Tensor<Scalar, 3> &get_multisite_mps() const;

    public:
    void                set_truncation_error(size_t pos, double error);
    void                set_truncation_error(double error);
    void                set_truncation_error_LC(double error);
    double              get_truncation_error(size_t pos) const;
    double              get_truncation_error() const;
    double              get_truncation_error_LC() const;
    double              get_truncation_error_midchain() const;
    std::vector<double> get_truncation_errors() const;
    std::vector<double> get_truncation_errors_active() const;

    size_t num_sites_truncated(double truncation_threshold) const;
    size_t num_bonds_reached_chi(long chi_level) const;
    bool   is_bond_limited(long chi_limit, double truncation_threshold) const;

    void clear_measurements(LogPolicy logPolicy = LogPolicy::NORMAL) const;
    void clear_cache(LogPolicy logPolicy = LogPolicy::NORMAL) const;
    void do_all_measurements() const;

    void tag_active_sites_normalized(bool tag) const;
    void tag_all_sites_normalized(bool tag) const;
    void tag_site_normalized(size_t pos, bool tag) const;
    bool is_normalized_on_all_sites() const;
    bool is_normalized_on_any_sites() const;
    bool is_normalized_on_active_sites() const;
    bool is_normalized_on_non_active_sites() const;
};
