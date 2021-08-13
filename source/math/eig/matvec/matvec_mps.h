#pragma once
#include "../enums.h"
#include <array>
#include <complex>
#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>
namespace tid {
    class ur;
}

struct primme_params;

template<class Scalar_>
class MatVecMps {
    public:
    using Scalar                                   = Scalar_;
    constexpr static bool         can_shift_invert = true;
    constexpr static bool         can_shift        = true;
    constexpr static bool         can_compress     = true;
    constexpr static eig::Storage storage          = eig::Storage::MPS;

    private:
    Eigen::Tensor<Scalar, 3> envL;
    Eigen::Tensor<Scalar, 3> envR;
    Eigen::Tensor<Scalar, 4> mpo;
    std::array<long, 3>      shape_mps;

    long      mps_size;
    eig::Form form = eig::Form::SYMM;
    eig::Side side = eig::Side::R;

    // Shift and shift-invert mode stuff
    std::complex<double> sigma         = 0.0;   // The shift
    bool                 readyShift    = false; // Flag to make sure the shift has occurred
    bool                 readyFactorOp = false; // Flag to make sure LU factorization has occurred
    bool                 readyCompress = false; // Flag to check if compression has occurred

    public:
    MatVecMps() = default;
    template<typename T>
    MatVecMps(const Eigen::Tensor<T, 3> &envL_, /*!< The left block tensor.  */
              const Eigen::Tensor<T, 3> &envR_, /*!< The right block tensor.  */
              const Eigen::Tensor<T, 4> &mpo_   /*!< The Hamiltonian MPO's  */
    );

    // Functions used in in Arpack++ solver
    [[nodiscard]] int rows() const { return static_cast<int>(mps_size); }; /*!< Linear size\f$d^2 \times \chi_L \times \chi_R \f$  */
    [[nodiscard]] int cols() const { return static_cast<int>(mps_size); }; /*!< Linear size\f$d^2 \times \chi_L \times \chi_R \f$  */

    void FactorOP();                                  //  Would normally factor (A-sigma*I) into PLU --> here it does nothing
    void MultOPv(Scalar_ *mps_in_, Scalar_ *mps_out); //  Computes the matrix-vector product x_out <- inv(A-sigma*I)*x_in.
    void MultAx(Scalar_ *mps_in_, Scalar_ *mps_out_); //  Computes the matrix-vector multiplication x_out <- A*x_in.
    void MultAx(Scalar_ *mps_in, Scalar_ *mps_out, Scalar_ *mpo_ptr, Scalar_ *envL_ptr, Scalar_ *envR_ptr, std::array<long, 3> shape_mps_,
                std::array<long, 4> shape_mpo_); //  Computes the matrix-vector multiplication x_out <- A*x_in.
    void MultAx(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *err);

    void compress(); //  Compresses the MPO virtual bond so that x_out <-- A*x_in can happen faster

    // Various utility functions
    int  counter = 0;
    void print() const;
    void reset();
    void set_shift(std::complex<double> sigma_);
    void set_mode(eig::Form form_);
    void set_side(eig::Side side_);

    template<typename T = eig::cplx>
    [[nodiscard]] T         get_shift() const;
    [[nodiscard]] eig::Form get_form() const;
    [[nodiscard]] eig::Side get_side() const;
    [[nodiscard]] eig::Type get_type() const;

    [[nodiscard]] const Eigen::Tensor<Scalar, 4> &get_mpo() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &get_envL() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3> &get_envR() const;
    [[nodiscard]] std::array<long, 3>             get_shape_mps() const;
    [[nodiscard]] std::array<long, 4>             get_shape_mpo() const;
    [[nodiscard]] std::array<long, 3>             get_shape_envL() const;
    [[nodiscard]] std::array<long, 3>             get_shape_envR() const;

    [[nodiscard]] bool isReadyFactorOp() const;
    [[nodiscard]] bool isReadyShift() const;
    [[nodiscard]] bool isReadyCompress() const;

    // Profiling
    std::unique_ptr<tid::ur> t_factorOP;
    std::unique_ptr<tid::ur> t_multOPv;
    std::unique_ptr<tid::ur> t_multAx;
};
