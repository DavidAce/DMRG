#pragma once
#include "../enums.h"
#include <array>
#include <complex>
#include <memory>
#include <vector>
namespace tid {
    class ur;
}

class primme_params;

template<class Scalar_>
class MatVecMPO {
    public:
    using Scalar                                   = Scalar_;
    constexpr static bool         can_shift_invert = true;
    constexpr static bool         can_shift        = true;
    constexpr static bool         can_compress     = true;
    constexpr static eig::Storage storage          = eig::Storage::TENSOR;

    private:
    const Scalar_      *envL;
    const Scalar_      *envR;
    const Scalar_      *mpo;
    std::array<long, 3> shape_mps;
    std::array<long, 4> shape_mpo;
    std::array<long, 3> shape_envL;
    std::array<long, 3> shape_envR;
    std::vector<Scalar> mpo_internal;
    std::vector<Scalar> envL_internal;
    std::vector<Scalar> envR_internal;

    long      mps_size;
    eig::Form form = eig::Form::SYMM;
    eig::Side side = eig::Side::R;

    // Shift and shift-invert mode stuff
    std::complex<double> sigma         = 0.0;   // The real part of the shift
    bool                 readyShift    = false; // Flag to make sure the shift has occurred
    bool                 readyFactorOp = false; // Flag to make sure LU factorization has occurred
    bool                 readyCompress = false; // Flag to check if compression has occurred

    public:
    MatVecMPO(const Scalar_      *envL_,      /*!< The left block tensor.  */
              const Scalar_      *envR_,      /*!< The right block tensor.  */
              const Scalar_      *mpo_,       /*!< The Hamiltonian MPO's  */
              std::array<long, 3> shape_mps_, /*!< An array containing the dimensions of the multisite mps  */
              std::array<long, 4> shape_mpo_  /*!< An array containing the dimensions of the multisite mpo  */
    );

    // Functions used in in Arpack++ solver
    [[nodiscard]] int rows() const { return static_cast<int>(mps_size); }; /*!< Linear size\f$d^2 \times \chi_L \times \chi_R \f$  */
    [[nodiscard]] int cols() const { return static_cast<int>(mps_size); }; /*!< Linear size\f$d^2 \times \chi_L \times \chi_R \f$  */

    void FactorOP();                                  //  Would normally factor (A-sigma*I) into PLU --> here it does nothing
    void MultOPv(Scalar_ *mps_in_, Scalar_ *mps_out); //  Computes the matrix-vector product x_out <- inv(A-sigma*I)*x_in.
    void MultAx(Scalar_ *mps_in_, Scalar_ *mps_out_); //  Computes the matrix-vector multiplication x_out <- A*x_in.
    void MultAx(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *err);
    void compress(); //  Compresses the MPO virtual bond so that x_out <-- A*x_in can happen faster

    // Various utility functions
    int  counter = 0;
    void print() const;
    void set_shift(std::complex<double> sigma_);
    void set_mode(eig::Form form_);
    void set_side(eig::Side side_);

    [[nodiscard]] const eig::Form &get_form() const;
    [[nodiscard]] const eig::Side &get_side() const;

    [[nodiscard]] const Scalar       *get_mpo() const { return mpo; }
    [[nodiscard]] const Scalar       *get_envL() const { return envL; }
    [[nodiscard]] const Scalar       *get_envR() const { return envR; }
    [[nodiscard]] std::array<long, 4> get_shape_mpo() const { return shape_mpo; }
    [[nodiscard]] std::array<long, 3> get_shape_envL() const { return shape_envL; }
    [[nodiscard]] std::array<long, 3> get_shape_envR() const { return shape_envR; }

    [[nodiscard]] bool isReadyFactorOp() const { return readyFactorOp; }
    [[nodiscard]] bool isReadyShift() const { return readyShift; }
    [[nodiscard]] bool isReadyCompress() const { return readyCompress; }

    // Profiling
    std::unique_ptr<tid::ur> t_factorOP;
    std::unique_ptr<tid::ur> t_multOPv;
    std::unique_ptr<tid::ur> t_multAx;
};
