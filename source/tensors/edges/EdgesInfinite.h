#pragma once
#pragma once
#include <complex>
#include <config/enums.h>
#include <math/tenx/fwd_decl.h>
#include <memory>
#include <optional>

class EnvEne;
class EnvVar;

class EdgesInfinite {
    public:
    using Scalar = std::complex<double>;

    private:
    std::unique_ptr<EnvEne> eneL;
    std::unique_ptr<EnvEne> eneR;
    std::unique_ptr<EnvVar> varL;
    std::unique_ptr<EnvVar> varR;

    public:
    EdgesInfinite();
    ~EdgesInfinite();                                     // Read comment on implementation
    EdgesInfinite(EdgesInfinite &&other);                 // default move ctor
    EdgesInfinite &operator=(EdgesInfinite &&other);      // default move assign
    EdgesInfinite(const EdgesInfinite &other);            // copy ctor
    EdgesInfinite &operator=(const EdgesInfinite &other); // copy assign

    void                 initialize();
    void                 eject_edges();
    void                 assert_validity() const;
    [[nodiscard]] size_t get_length() const;
    size_t               get_position() const; // pos of eneL or varL
    [[nodiscard]] bool   is_real() const;
    [[nodiscard]] bool   has_nan() const;

    // This is a reference wrapper for an edge pair
    template<typename env_type>
    struct env_pair {
        env_type &L;
        env_type &R;
        void      assert_validity() const;
        //        env_pair(env_type &L_, env_type &R_) : L(L_), R(R_) {}
    };

    [[nodiscard]] env_pair<const EnvEne> get_ene() const;
    [[nodiscard]] env_pair<const EnvVar> get_var() const;
    [[nodiscard]] env_pair<EnvEne>       get_ene();
    [[nodiscard]] env_pair<EnvVar>       get_var();

    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 3>> get_ene_blk() const;
    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 3>> get_var_blk() const;
    [[nodiscard]] env_pair<Eigen::Tensor<Scalar, 3>>       get_ene_blk();
    [[nodiscard]] env_pair<Eigen::Tensor<Scalar, 3>>       get_var_blk();
};
