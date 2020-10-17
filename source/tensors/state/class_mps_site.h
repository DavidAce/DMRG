//
// Created by david on 2019-07-06.
//

#pragma once

#include <complex>
#include <optional>
#include <general/nmspc_tensor_extra.h>

class class_mps_site {
    public:
    using Scalar = std::complex<double>;

    private:
    std::optional<Eigen::Tensor<Scalar, 3>>         M                   = std::nullopt; /*!< \f$M\f$ A or B tensor (note: not a Gamma tensor!) */
    std::optional<Eigen::Tensor<Scalar, 1>>         L                   = std::nullopt; /*!< \f$\Lambda\f$*/
    std::optional<Eigen::Tensor<Scalar, 1>>         LC                  = std::nullopt; /*!< \f$\Lambda_C\f$ Center lambda, if this is a center matrix*/
    mutable std::optional<Eigen::Tensor<Scalar, 3>> MC                  = std::nullopt;
    std::optional<size_t>                           position            = std::nullopt;
    double                                          truncation_error    = 0;
    double                                          truncation_error_LC = 0;

    public:
    ~class_mps_site(); // Read comment on implementation
    class_mps_site(const Eigen::Tensor<Scalar, 3> &M_, const Eigen::Tensor<Scalar, 1> &L_, size_t pos, double error = 0);
    class_mps_site(const Eigen::Tensor<Scalar, 3> &M_, std::optional<Eigen::Tensor<Scalar, 1>> L_, size_t pos, double error = 0);
    class_mps_site();                                           // ctor
    class_mps_site(const class_mps_site &other);                // default copy ctor
    class_mps_site(class_mps_site &&other);                     // default move ctor
    class_mps_site &operator=(class_mps_site &&other);          // default move assign
    class_mps_site &operator=(const class_mps_site &other);     // default copy assign

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
    [[nodiscard]] std::tuple<long, long, long>    get_dims() const;
    [[nodiscard]] long                            spin_dim() const;
    [[nodiscard]] long                            get_chiL() const;
    [[nodiscard]] long                            get_chiR() const;
    [[nodiscard]] size_t                          get_position() const;
    void                                          set_position(size_t position_);
    void                                          set_mps(const Eigen::Tensor<Scalar, 3> &M_, const Eigen::Tensor<Scalar, 1> &L_, double error = 0);

    void   set_M(const Eigen::Tensor<Scalar, 3> &M_);
    void   set_L(const Eigen::Tensor<Scalar, 1> &L_, double error = 0);
    void   set_LC(const Eigen::Tensor<Scalar, 1> &LC_, double error = 0);
    void   set_truncation_error(double error);
    void   set_truncation_error_LC(double error);
    double get_truncation_error() const;
    double get_truncation_error_LC() const;
    void   unset_LC() {
        LC.reset();
        MC.reset();
    }
    void merge_mps(const class_mps_site & other);
    void apply_mpo(const Eigen::Tensor<Scalar, 4> &mpo);
    void apply_mpo(const Eigen::Tensor<Scalar, 2> &mpo);
};
