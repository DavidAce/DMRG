#pragma once
#pragma once
#include <complex>
#include <config/enums.h>
#include <list>
#include <memory>
#include <optional>
#include <general/eigen_tensor_fwd_decl.h>

class class_env_ene;
class class_env_var;

class class_edges_infinite {
    public:
    using Scalar = std::complex<double>;

    private:
    std::unique_ptr<class_env_ene> eneL;
    std::unique_ptr<class_env_ene> eneR;
    std::unique_ptr<class_env_var> varL;
    std::unique_ptr<class_env_var> varR;

    public:
    class_edges_infinite();
    ~class_edges_infinite();                                                // Read comment on implementation
    class_edges_infinite(class_edges_infinite &&other);                     // default move ctor
    class_edges_infinite &operator=(class_edges_infinite &&other);          // default move assign
    class_edges_infinite(const class_edges_infinite &other);                // copy ctor
    class_edges_infinite &operator=(const class_edges_infinite &other);     // copy assign

    void initialize();
    void eject_edges();
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

    [[nodiscard]] env_pair<const class_env_ene> get_ene() const;
    [[nodiscard]] env_pair<const class_env_var> get_var() const;
    [[nodiscard]] env_pair<class_env_ene>       get_ene();
    [[nodiscard]] env_pair<class_env_var>       get_var();

    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 3>> get_ene_blk() const;
    [[nodiscard]] env_pair<const Eigen::Tensor<Scalar, 3>> get_var_blk() const;
    [[nodiscard]] env_pair<Eigen::Tensor<Scalar, 3>>       get_ene_blk();
    [[nodiscard]] env_pair<Eigen::Tensor<Scalar, 3>>       get_var_blk();
};
