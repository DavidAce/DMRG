//
// Created by david on 2019-01-29.
//

#pragma once

#include <complex>
#include <config/enums.h>
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

    int                       direction = 1;
    mutable Cache             cache;
    mutable std::vector<bool> tag_normalized_sites;
    std::string               name;
    AlgorithmType             algo = AlgorithmType::ANY;

    public:
    std::vector<std::unique_ptr<class_mps_site>> mps_sites;
    std::vector<size_t>                          active_sites;
    mutable state_measure_finite                 measurements;

    public:
    class_state_finite();
    ~class_state_finite();                                          // Read comment on implementation
    class_state_finite(class_state_finite &&other);                 // default move ctor
    class_state_finite &operator=(class_state_finite &&other);      // default move assign
    class_state_finite(const class_state_finite &other);            // copy ctor
    class_state_finite &operator=(const class_state_finite &other); // copy assign

    void initialize(ModelType modeltype, size_t model_size, size_t position = 0);

    void                      set_name(const std::string &statename);
    [[nodiscard]] std::string get_name() const;

    void                        set_algorithm(const AlgorithmType &algo_type);
    [[nodiscard]] AlgorithmType get_algorithm() const;

    const Eigen::Tensor<Scalar, 1> &midchain_bond() const;
    const Eigen::Tensor<Scalar, 1> &current_bond() const;

    template<typename T = size_t>
    [[nodiscard]] T get_length() const;
    template<typename T = size_t>
    [[nodiscard]] T get_position() const;

    void                                   set_positions();
    void                                   flip_direction();
    [[nodiscard]] int                      get_direction() const;
    [[nodiscard]] std::vector<std::string> get_labels() const;
    [[nodiscard]] std::array<long, 3>      dimensions_1site() const;
    [[nodiscard]] std::array<long, 3>      dimensions_2site() const;
    [[nodiscard]] std::array<long, 3>      dimensions_nsite() const;
    [[nodiscard]] long                     size_1site() const;
    [[nodiscard]] long                     size_2site() const;
    [[nodiscard]] long                     size_nsite() const;
    [[nodiscard]] long                     find_largest_chi() const;
    [[nodiscard]] bool                     position_is_the_middle() const;
    [[nodiscard]] bool                     position_is_the_middle_any_direction() const;
    [[nodiscard]] bool                     position_is_outward_edge_left(size_t nsite = 1) const;
    [[nodiscard]] bool                     position_is_outward_edge_right(size_t nsite = 1) const;
    [[nodiscard]] bool                     position_is_outward_edge(size_t nsite = 1) const;
    [[nodiscard]] bool                     position_is_inward_edge_left(size_t nsite = 1) const;
    [[nodiscard]] bool                     position_is_inward_edge_right(size_t nsite = 1) const;
    [[nodiscard]] bool                     position_is_inward_edge(size_t nsite = 1) const;
    [[nodiscard]] bool                     position_is_at(long pos) const;
    [[nodiscard]] bool                     position_is_at(long pos, int dir) const;
    [[nodiscard]] bool                     position_is_at(long pos, int dir, bool isCenter) const;
    [[nodiscard]] bool                     has_center_point() const;
    [[nodiscard]] bool                     is_real() const;
    [[nodiscard]] bool                     has_nan() const;

    void assert_validity() const;

    template<typename T = size_t>
    const class_mps_site &get_mps_site(T pos) const;
    template<typename T = size_t>
    class_mps_site &            get_mps_site(T pos);
    const class_mps_site &      get_mps_site() const;
    class_mps_site &            get_mps_site();
    std::vector<class_mps_site> get_mps_sites(const std::vector<size_t> &sites) const;
    void                        set_mps_sites(const std::vector<class_mps_site> &new_mps);
    // For multisite
    std::array<long, 3> active_dimensions() const;
    long                active_problem_size() const;
    std::vector<long>   get_spin_dims(const std::vector<size_t> &sites) const;
    std::vector<long>   get_spin_dims() const;

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

    std::vector<size_t> get_active_ids() const;
};
