#pragma once

#include "MpsStash.h"
#include <complex>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>

class MpsSite {
    public:
    using Scalar = std::complex<double>;

    private:
    std::optional<Eigen::Tensor<Scalar, 3>>                M                   = std::nullopt; /*!< \f$M\f$ A or B tensor (note: not a Gamma tensor!) */
    std::optional<Eigen::Tensor<Scalar, 1>>                L                   = std::nullopt; /*!< \f$\Lambda\f$*/
    std::optional<Eigen::Tensor<Scalar, 1>>                LC                  = std::nullopt; /*!< \f$\Lambda_C\f$ Center lambda, if this is a center matrix*/
    mutable std::optional<Eigen::Tensor<Scalar, 3>>        MC                  = std::nullopt;
    std::optional<size_t>                                  position            = std::nullopt;
    double                                                 truncation_error    = 0;
    double                                                 truncation_error_LC = 0;
    std::string                                            label;
    mutable std::optional<std::size_t>                     unique_id;
    mutable std::optional<stash<Eigen::Tensor<Scalar, 3>>> U_stash = std::nullopt; /*!< \f$U\f$ A "U" matrix from SVD stored temporarily  */
    mutable std::optional<stash<Eigen::Tensor<Scalar, 1>>> S_stash = std::nullopt; /*!< \f$S\f$ A "S" matrix from SVD stored temporarily  */
    mutable std::optional<stash<Eigen::Tensor<Scalar, 1>>> C_stash = std::nullopt; /*!< \f$S\f$ A "C" matrix from SVD stored temporarily  */
    mutable std::optional<stash<Eigen::Tensor<Scalar, 3>>> V_stash = std::nullopt; /*!< \f$V\f$ A "V" matrix from SVD stored temporarily  */

    public:
    ~MpsSite(); // Read comment on implementation
    MpsSite(const Eigen::Tensor<Scalar, 3> &M_, const Eigen::Tensor<Scalar, 1> &L_, size_t pos, double error, const std::string &label_);
    MpsSite(const Eigen::Tensor<Scalar, 3> &M_, std::optional<Eigen::Tensor<Scalar, 1>> L_, size_t pos, double error, const std::string &label_);
    MpsSite();                                // ctor
    MpsSite(const MpsSite &other);            // default copy ctor
    MpsSite(MpsSite &&other);                 // default move ctor
    MpsSite &operator=(MpsSite &&other);      // default move assign
    MpsSite &operator=(const MpsSite &other); // default copy assign

    [[nodiscard]] bool                            is_real() const;
    [[nodiscard]] bool                            has_nan() const;
    void                                          assert_validity() const;
    void                                          assert_dimensions() const;
    void                                          assert_identity() const;
    [[nodiscard]] bool                            isCenter() const;
    [[nodiscard]] bool                            has_L() const;
    [[nodiscard]] bool                            has_M() const;
    [[nodiscard]] bool                            has_LC() const;
    [[nodiscard]] Eigen::DSizes<long, 3>          dimensions() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &get_M_bare() const; /*!< Gets the A or B matrix without LC attached */
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &get_M() const;      /*!< Gets A or B matrix. If this is a center (A matrix) it attaches LC */
    [[nodiscard]] const Eigen::Tensor<Scalar, 1> &get_L() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 1> &get_LC() const;
    [[nodiscard]] Eigen::Tensor<Scalar, 3>       &get_M_bare();
    [[nodiscard]] Eigen::Tensor<Scalar, 3>       &get_M();
    [[nodiscard]] Eigen::Tensor<Scalar, 1>       &get_L();
    [[nodiscard]] Eigen::Tensor<Scalar, 1>       &get_LC();
    [[nodiscard]] double                          get_truncation_error() const;
    [[nodiscard]] double                          get_truncation_error_LC() const;
    [[nodiscard]] std::string                     get_label() const;
    [[nodiscard]] std::tuple<long, long, long>    get_dims() const;
    [[nodiscard]] long                            spin_dim() const;
    [[nodiscard]] long                            get_chiL() const;
    [[nodiscard]] long                            get_chiR() const;

    template<typename T = size_t>
    [[nodiscard]] T get_position() const;

    void set_M(const Eigen::Tensor<Scalar, 3> &M_);
    void set_L(const Eigen::Tensor<Scalar, 1> &L_, double error = 0);
    void set_L(const std::pair<Eigen::Tensor<Scalar, 1>, double> &L_and_error);
    void set_LC(const Eigen::Tensor<Scalar, 1> &LC_, double error = 0);
    void set_LC(const std::pair<Eigen::Tensor<Scalar, 1>, double> &LC_and_error);
    void set_truncation_error(double error);
    void set_truncation_error_LC(double error);
    void set_label(const std::string &label_);
    void set_position(size_t position_);
    void set_mps(const Eigen::Tensor<Scalar, 3> &M_, const Eigen::Tensor<Scalar, 1> &L_, double error, const std::string &label_);

    void unset_LC();
    void unset_L();
    void fuse_mps(const MpsSite &other);
    void apply_mpo(const Eigen::Tensor<Scalar, 4> &mpo);
    void apply_mpo(const Eigen::Tensor<Scalar, 2> &mpo);

    void stash_U(const Eigen::Tensor<Scalar, 3> &U, size_t dst) const;
    void stash_S(const Eigen::Tensor<Scalar, 1> &S, double error, size_t dst) const;
    void stash_S(const std::pair<Eigen::Tensor<Scalar, 1>, double> &S_and_error, size_t dst) const;
    void stash_C(const Eigen::Tensor<Scalar, 1> &S, double error, size_t dst) const;
    void stash_C(const std::pair<Eigen::Tensor<Scalar, 1>, double> &S_and_error, size_t dst) const;
    void stash_V(const Eigen::Tensor<Scalar, 3> &V, size_t dst) const;
    void drop_stash() const;
    void take_stash(const MpsSite &other);

    std::optional<stash<Eigen::Tensor<Scalar, 3>>> &get_U_stash() const;
    std::optional<stash<Eigen::Tensor<Scalar, 1>>> &get_S_stash() const;
    std::optional<stash<Eigen::Tensor<Scalar, 1>>> &get_C_stash() const;
    std::optional<stash<Eigen::Tensor<Scalar, 3>>> &get_V_stash() const;

    std::size_t get_unique_id() const;
};
