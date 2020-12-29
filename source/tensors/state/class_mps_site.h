//
// Created by david on 2019-07-06.
//

#pragma once

#include <complex>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>

class class_mps_site {
    public:
    using Scalar = std::complex<double>;

    private:
    std::optional<Eigen::Tensor<Scalar, 3>>         M                   = std::nullopt; /*!< \f$M\f$ A or B tensor (note: not a Gamma tensor!) */
    std::optional<Eigen::Tensor<Scalar, 1>>         L                   = std::nullopt; /*!< \f$\Lambda\f$*/
    std::optional<Eigen::Tensor<Scalar, 1>>         LC                  = std::nullopt; /*!< \f$\Lambda_C\f$ Center lambda, if this is a center matrix*/
    mutable std::optional<Eigen::Tensor<Scalar, 3>> MC                  = std::nullopt;
    std::optional<long>                             position            = std::nullopt;
    double                                          truncation_error    = 0;
    double                                          truncation_error_LC = 0;
    std::string                                     label;
    mutable std::optional<std::size_t>              unique_id;
    mutable std::optional<Eigen::Tensor<Scalar, 3>> U_stash                  = std::nullopt; /*!< \f$U\f$ A "U" matrix from SVD stored temporarily  */
    mutable std::optional<Eigen::Tensor<Scalar, 1>> S_stash                  = std::nullopt; /*!< \f$S\f$ A "S" matrix from SVD stored temporarily  */
    mutable std::optional<Eigen::Tensor<Scalar, 3>> V_stash                  = std::nullopt; /*!< \f$V\f$ A "V" matrix from SVD stored temporarily  */
    mutable std::optional<double>                   truncation_error_S_stash = std::nullopt;

    public:
    ~class_mps_site(); // Read comment on implementation
    class_mps_site(const Eigen::Tensor<Scalar, 3> &M_, const Eigen::Tensor<Scalar, 1> &L_, size_t pos, double error, std::string  label_);
    class_mps_site(const Eigen::Tensor<Scalar, 3> &M_, std::optional<Eigen::Tensor<Scalar, 1>> L_, size_t pos, double error, std::string  label_);
    class_mps_site();                                       // ctor
    class_mps_site(const class_mps_site &other);            // default copy ctor
    class_mps_site(class_mps_site &&other) noexcept ;                 // default move ctor
    class_mps_site &operator=(class_mps_site &&other) noexcept ;      // default move assign
    class_mps_site &operator=(const class_mps_site &other); // default copy assign

    [[nodiscard]] bool                            is_real() const;
    [[nodiscard]] bool                            has_nan() const;
    void                                          assert_validity() const;
    [[nodiscard]] bool                            isCenter() const;
    [[nodiscard]] bool                            has_L() const;
    [[nodiscard]] bool                            has_M() const;
    [[nodiscard]] bool                            has_LC() const;
    [[nodiscard]] Eigen::DSizes<long, 3>          dimensions() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &get_M_bare() const; /*!< Gets the A or B matrix without LC attached */
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &get_M() const;      /*!< Gets A or B matrix. If this is a center (A matrix) it attaches LC */
    [[nodiscard]] const Eigen::Tensor<Scalar, 1> &get_L() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 1> &get_LC() const;
    [[nodiscard]] Eigen::Tensor<Scalar, 3> &      get_M_bare();
    [[nodiscard]] Eigen::Tensor<Scalar, 3> &      get_M();
    [[nodiscard]] Eigen::Tensor<Scalar, 1> &      get_L();
    [[nodiscard]] Eigen::Tensor<Scalar, 1> &      get_LC();
    [[nodiscard]] double                          get_truncation_error() const;
    [[nodiscard]] double                          get_truncation_error_LC() const;
    [[nodiscard]] std::string                     get_label() const;
    [[nodiscard]] std::tuple<long, long, long>    get_dims() const;
    [[nodiscard]] long                            spin_dim() const;
    [[nodiscard]] long                            get_chiL() const;
    [[nodiscard]] long                            get_chiR() const;

    template<typename T = size_t>  [[nodiscard]] T get_position() const;

    void set_M(const Eigen::Tensor<Scalar, 3> &M_);
    void set_L(const Eigen::Tensor<Scalar, 1> &L_, double error = 0);
    void set_L(const std::pair<Eigen::Tensor<Scalar, 1>, double> &L_and_error);
    void set_LC(const Eigen::Tensor<Scalar, 1> &LC_, double error = 0);
    void set_LC(const std::pair<Eigen::Tensor<Scalar, 1>, double> &LC_and_error);
    void set_truncation_error(double error);
    void set_truncation_error_LC(double error);
    void set_label(const std::string &label_);
    void set_position(long position_);
    void set_mps(const Eigen::Tensor<Scalar, 3> &M_, const Eigen::Tensor<Scalar, 1> &L_, double error, const std::string & label_);

    void unset_LC();
    void merge_mps(const class_mps_site &other);
    void apply_mpo(const Eigen::Tensor<Scalar, 4> &mpo);
    void apply_mpo(const Eigen::Tensor<Scalar, 2> &mpo);

    void stash_U(const Eigen::Tensor<Scalar, 3> &U) const;
    void stash_S(const Eigen::Tensor<Scalar, 1> &S, double error) const;
    void stash_S(const std::pair<Eigen::Tensor<Scalar, 1>, double> &S_and_error) const;
    void stash_V(const Eigen::Tensor<Scalar, 3> &V) const;

    bool has_stash_U() const;
    bool has_stash_S() const;
    bool has_stash_V() const;


    Eigen::Tensor<Scalar, 3>                    unstash_U() const;
    std::pair<Eigen::Tensor<Scalar, 1>, double> unstash_S() const;
    Eigen::Tensor<Scalar, 3>                    unstash_V() const;
    void                                        unstash() const;
    void                                        merge_stash(const class_mps_site &other);

    std::size_t get_unique_id() const;
};
