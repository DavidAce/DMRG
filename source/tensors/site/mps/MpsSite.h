#pragma once
#include "math/float.h"
#include "MpsStash.h"
#include <complex>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>


template<typename T>
concept is_valid_tensor3 = std::is_same_v<T, Eigen::Tensor<cx64, 3>> ||                   //
                           std::is_same_v<T, Eigen::Tensor<fp64, 3>> ||                   //
                           std::is_same_v<T, Eigen::TensorMap<Eigen::Tensor<cx64, 3>>> || //
                           std::is_same_v<T, Eigen::TensorMap<Eigen::Tensor<fp64, 3>>>;   //
//                                std::is_same_v<T, Eigen::TensorMap<const Eigen::Tensor<cx64, 3>>> || //
//                                std::is_same_v<T, Eigen::TensorMap<const Eigen::Tensor<fp64, 3>>>;

template<typename T>
concept is_valid_tensor1 = std::is_same_v<T, Eigen::Tensor<cx64, 1>> ||                   //
                           std::is_same_v<T, Eigen::Tensor<fp64, 1>> ||                   //
                           std::is_same_v<T, Eigen::TensorMap<Eigen::Tensor<cx64, 1>>> || //
                           std::is_same_v<T, Eigen::TensorMap<Eigen::Tensor<fp64, 1>>>;   //||       //
//                           std::is_same_v<T, Eigen::TensorMap<const Eigen::Tensor<cx64, 1>>> || //
//                           std::is_same_v<T, Eigen::TensorMap<const Eigen::Tensor<fp64, 1>>>;

class MpsSite {
    public:
    using value_type = cx64;

    private:
    std::optional<Eigen::Tensor<cx64, 3>>                M                   = std::nullopt; /*!< \f$M\f$ A or B tensor (note: not a Gamma tensor!) */
    std::optional<Eigen::Tensor<cx64, 1>>                L                   = std::nullopt; /*!< \f$\Lambda\f$*/
    std::optional<Eigen::Tensor<cx64, 1>>                LC                  = std::nullopt; /*!< \f$\Lambda_C\f$ Center lambda, if this is a center matrix*/
    mutable std::optional<Eigen::Tensor<cx64, 3>>        MC                  = std::nullopt;
    std::optional<size_t>                                position            = std::nullopt;
    double                                               truncation_error    = -2.0;
    double                                               truncation_error_LC = -2.0;
    std::string                                          label;
    mutable std::optional<std::size_t>                   unique_id;
    mutable std::optional<stash<Eigen::Tensor<cx64, 3>>> U_stash = std::nullopt; /*!< \f$U\f$ A "U" matrix from SVD stored temporarily  */
    mutable std::optional<stash<Eigen::Tensor<cx64, 1>>> S_stash = std::nullopt; /*!< \f$S\f$ A "S" matrix from SVD stored temporarily  */
    mutable std::optional<stash<Eigen::Tensor<cx64, 1>>> C_stash = std::nullopt; /*!< \f$S\f$ A "C" matrix from SVD stored temporarily  */
    mutable std::optional<stash<Eigen::Tensor<cx64, 3>>> V_stash = std::nullopt; /*!< \f$V\f$ A "V" matrix from SVD stored temporarily  */

    public:
    ~MpsSite(); // Read comment on implementation
    explicit MpsSite(const Eigen::Tensor<cx64, 3> &M_, const std::optional<Eigen::Tensor<cx64, 1>> &L_, size_t pos, double error, std::string_view label_);
    explicit MpsSite(const Eigen::Tensor<fp64, 3> &M_, const std::optional<Eigen::Tensor<fp64, 1>> &L_, size_t pos, double error, std::string_view label_);
    template<typename T3, typename T1>
    requires is_valid_tensor3<T3> && is_valid_tensor1<T1>
    MpsSite(const T3 &M_, const T1 &L_, size_t pos, double error, std::string_view label_);
    template<typename T3>
    requires is_valid_tensor3<T3>
    MpsSite(const T3 &M_, size_t pos, std::string_view label_);
    MpsSite();                                    // ctor
    MpsSite(const MpsSite &other);                // default copy ctor
    MpsSite(MpsSite &&other) noexcept;            // default move ctor
    MpsSite &operator=(MpsSite &&other) noexcept; // default move assign
    MpsSite &operator=(const MpsSite &other);     // default copy assign

    [[nodiscard]] bool                          is_real() const;
    [[nodiscard]] bool                          has_nan() const;
    [[nodiscard]] bool                          is_normalized(double prec = 1e-10) const;
    void                                        assert_validity() const;
    void                                        assert_dimensions() const;
    void                                        assert_normalized(double prec = 1e-10) const;
    [[nodiscard]] bool                          isCenter() const;
    [[nodiscard]] bool                          has_L() const;
    [[nodiscard]] bool                          has_M() const;
    [[nodiscard]] bool                          has_LC() const;
    [[nodiscard]] Eigen::DSizes<long, 3>        dimensions() const;
    [[nodiscard]] const Eigen::Tensor<cx64, 3> &get_M_bare() const; /*!< Gets the A or B matrix without LC attached */
    [[nodiscard]] const Eigen::Tensor<cx64, 3> &get_M() const;      /*!< Gets A or B matrix. If this is a center (A matrix) it attaches LC */
    [[nodiscard]] const Eigen::Tensor<cx64, 1> &get_L() const;
    [[nodiscard]] const Eigen::Tensor<cx64, 1> &get_LC() const;
    [[nodiscard]] Eigen::Tensor<cx64, 3>       &get_M_bare();
    [[nodiscard]] Eigen::Tensor<cx64, 3>       &get_M();
    [[nodiscard]] Eigen::Tensor<cx64, 1>       &get_L();
    [[nodiscard]] Eigen::Tensor<cx64, 1>       &get_LC();
    [[nodiscard]] double                        get_truncation_error() const;
    [[nodiscard]] double                        get_truncation_error_LC() const;
    [[nodiscard]] std::string_view              get_label() const;
    [[nodiscard]] std::string                   get_tag() const;
    [[nodiscard]] std::tuple<long, long, long>  get_dims() const;
    [[nodiscard]] long                          spin_dim() const;
    [[nodiscard]] long                          get_chiL() const;
    [[nodiscard]] long                          get_chiR() const;

    template<typename T = size_t>
    [[nodiscard]] T get_position() const;

    template<typename T>
    [[nodiscard]] bool is_at_position(T pos) const;

    template<typename T3>
    requires is_valid_tensor3<T3>
    void set_M(const T3 &M_);
    template<typename T1>
    requires is_valid_tensor1<T1>
    void set_L(const T1 &L_, double error = 0);
    template<typename T1>
    requires is_valid_tensor1<T1>
    void set_LC(const T1 &LC_, double error = 0);

    void set_L(const std::pair<Eigen::Tensor<cx64, 1>, double> &L_and_error);
    void set_L(const std::pair<Eigen::Tensor<fp64, 1>, double> &L_and_error);
    void set_LC(const std::pair<Eigen::Tensor<cx64, 1>, double> &LC_and_error);
    void set_LC(const std::pair<Eigen::Tensor<fp64, 1>, double> &LC_and_error);
    void set_truncation_error(double error);
    void set_truncation_error_LC(double error);
    void set_label(std::string_view label_);
    void set_position(size_t position_);
    void set_mps(const Eigen::Tensor<cx64, 3> &M_, const Eigen::Tensor<cx64, 1> &L_, double error, std::string_view label_);
    void set_mps(const Eigen::Tensor<cx64, 3> &M_, const Eigen::Tensor<fp64, 1> &L_, double error, std::string_view label_);

    void unset_LC();
    void unset_L();
    void unset_truncation_error();
    void unset_truncation_error_LC();
    void fuse_mps(const MpsSite &other);
    void apply_mpo(const Eigen::Tensor<cx64, 4> &mpo, bool adjoint = false);
    void apply_mpo(const Eigen::Tensor<cx64, 2> &mpo, bool adjoint = false);

    template<typename T3>
    requires is_valid_tensor3<T3>
    void stash_U(const T3 &U, size_t dst) const;
    template<typename T1>
    requires is_valid_tensor1<T1>
    void stash_S(const T1 &S, double error, size_t dst) const;
    template<typename T1>
    requires is_valid_tensor1<T1>
    void stash_C(const T1 &S, double error, size_t dst) const;
    template<typename T3>
    requires is_valid_tensor3<T3>
    void stash_V(const T3 &V, size_t dst) const;

    void stash_S(const std::pair<Eigen::Tensor<cx64, 1>, double> &S_and_error, size_t dst) const;
    void stash_C(const std::pair<Eigen::Tensor<cx64, 1>, double> &S_and_error, size_t dst) const;
    void drop_stash() const;
    void drop_stashed_errors() const;
    void take_stash(const MpsSite &other);

    std::optional<stash<Eigen::Tensor<cx64, 3>>> &get_U_stash() const;
    std::optional<stash<Eigen::Tensor<cx64, 1>>> &get_S_stash() const;
    std::optional<stash<Eigen::Tensor<cx64, 1>>> &get_C_stash() const;
    std::optional<stash<Eigen::Tensor<cx64, 3>>> &get_V_stash() const;

    void convert_AL_to_A(const Eigen::Tensor<cx64, 1> &LR);
    void convert_LB_to_B(const Eigen::Tensor<cx64, 1> &LL);

    std::size_t get_unique_id() const;
};
