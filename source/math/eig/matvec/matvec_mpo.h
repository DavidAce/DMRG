#pragma once
#include "../enums.h"
#include "math/float.h"
#include <array>
#include <complex>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

namespace tid {
    class ur;
}

struct primme_params;

template<typename T>
class MatVecMPO {
    static_assert(std::is_same_v<T, real> or std::is_same_v<T, cplx>);

    public:
    using Scalar     = T;
    using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorType = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    constexpr static bool         can_shift_invert = true;
    constexpr static bool         can_shift        = true;
    constexpr static bool         can_compress     = true;
    constexpr static eig::Storage storage          = eig::Storage::MPS;
    eig::Factorization            factorization    = eig::Factorization::NONE;

    private:
    Eigen::Tensor<T, 3> envL;
    Eigen::Tensor<T, 3> envR;
    Eigen::Tensor<T, 4> mpo;
    std::array<long, 3> shape_mps;
    long                size_mps;
    eig::Form           form = eig::Form::SYMM;
    eig::Side           side = eig::Side::R;

    // Shift and shift-invert mode stuff
    std::complex<double> sigma         = cplx(0.0, 0.0); // The shift
    bool                 readyShift    = false;          // Flag to make sure the shift has occurred
    bool                 readyFactorOp = false;          // Flag to make sure LU factorization has occurred
    bool                 readyCompress = false;          // Flag to check if compression has occurred

    Eigen::LDLT<MatrixType>         ldlt; // Stores the ldlt matrix factorization on shift-invert
    Eigen::LLT<MatrixType>          llt;  // Stores the llt matrix factorization on shift-invert
    Eigen::PartialPivLU<MatrixType> lu;   // Stores the lu matrix factorization on shift-invert

    public:
    MatVecMPO() = default;
    template<typename S>
    MatVecMPO(const Eigen::Tensor<S, 3> &envL_, /*!< The left block tensor.  */
              const Eigen::Tensor<S, 3> &envR_, /*!< The right block tensor.  */
              const Eigen::Tensor<S, 4> &mpo_   /*!< The Hamiltonian MPO's  */
    );
    // Functions used in Arpack++ solver
    [[nodiscard]] int rows() const; /*!< Linear size\f$d^2 \times \chi_L \times \chi_R \f$  */
    [[nodiscard]] int cols() const; /*!< Linear size\f$d^2 \times \chi_L \times \chi_R \f$  */

    void FactorOP();                      //  Would normally factor (A-sigma*I) into PLU --> here it does nothing
    void MultOPv(T *mps_in_, T *mps_out); //  Computes the matrix-vector product x_out <- inv(A-sigma*I)*x_in.
    void MultOPv(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *err);
    void MultAx(T *mps_in_, T *mps_out_); //  Computes the matrix-vector multiplication x_out <- A*x_in.
    void MultAx(T *mps_in, T *mps_out, T *mpo_ptr, T *envL_ptr, T *envR_ptr, std::array<long, 3> shape_mps_,
                std::array<long, 4> shape_mpo_); //  Computes the matrix-vector multiplication x_out <- A*x_in.
    void MultAx(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *err);

    void compress(); //  Compresses the MPO virtual bond so that x_out <-- A*x_in can happen faster

    // Various utility functions
    long num_mv = 0;
    long num_op = 0;
    void print() const;
    void reset();
    void set_shift(std::complex<double> shift);
    void set_mode(eig::Form form_);
    void set_side(eig::Side side_);
    void set_readyCompress(bool compressed);

    [[nodiscard]] T                               get_shift() const;
    [[nodiscard]] eig::Form                       get_form() const;
    [[nodiscard]] eig::Side                       get_side() const;
    [[nodiscard]] eig::Type                       get_type() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 4> &get_mpo() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &get_envL() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &get_envR() const;
    [[nodiscard]] long                            get_size() const;
    [[nodiscard]] std::array<long, 3>             get_shape_mps() const;
    [[nodiscard]] std::array<long, 4>             get_shape_mpo() const;
    [[nodiscard]] std::array<long, 3>             get_shape_envL() const;
    [[nodiscard]] std::array<long, 3>             get_shape_envR() const;
    [[nodiscard]] Eigen::Tensor<Scalar, 6>        get_tensor() const;
    [[nodiscard]] MatrixType                      get_matrix() const;

    [[nodiscard]] bool isReadyFactorOp() const;
    [[nodiscard]] bool isReadyShift() const;
    [[nodiscard]] bool isReadyCompress() const;

    // Timers
    std::unique_ptr<tid::ur> t_factorOP;
    std::unique_ptr<tid::ur> t_genMat;
    std::unique_ptr<tid::ur> t_multOPv;
    std::unique_ptr<tid::ur> t_multAx;
};
