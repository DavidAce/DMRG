#pragma once
#include "../enums.h"
#include "math/float.h"
#include <array>
#include <complex>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

namespace tid {
    class ur;
}
class MpoSite;

template<typename T>
struct env_pair;

struct primme_params;

template<typename T>
class MatVecMPOS {
    static_assert(std::is_same_v<T, real> or std::is_same_v<T, cplx>);

    public:
    using Scalar         = T;
    using MatrixType     = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorType     = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    using SparseType     = Eigen::SparseMatrix<T>;
    using SparseTypeRowM = Eigen::SparseMatrix<T, Eigen::RowMajor>;

    constexpr static bool         can_shift_invert = false;
    constexpr static bool         can_shift        = false;
    constexpr static bool         can_precondition = true;
    constexpr static eig::Storage storage          = eig::Storage::MPS;
    eig::Factorization            factorization    = eig::Factorization::NONE;
    eig::Preconditioner           preconditioner   = eig::Preconditioner::NONE;

    private:
    std::vector<Eigen::Tensor<T, 4>>                mpos, mpos_shf;
    Eigen::Tensor<T, 3>                             mpoTL, mpoLR; // Top left and lower right after svd.
    Eigen::Tensor<T, 3>                             mpoTR, mpoLL; // Top right and lower left after svd.
    Eigen::Tensor<T, 3>                             envL;
    Eigen::Tensor<T, 3>                             envR;
    std::array<long, 3>                             shape_mps;
    long                                            size_mps;
    std::vector<long>                               spindims;
    eig::Form                                       form = eig::Form::SYMM;
    eig::Side                                       side = eig::Side::R;
    VectorType                                      diagonal;   // The diagonal elements of the matrix, used in the diagonal and tridiagonal preconditioners
    VectorType                                      diaglower;  // The sub-diagonal elements of the matrix, used in the tridiagonal preconditioner
    VectorType                                      diagupper;  // The super-diagonal elements of the matrix, used in the tridiagonal preconditioner
    VectorType                                      diagtemp;   // Scratch memory for the tridiagonal solver
    SparseType                                      diagband;   // The diagonal band stored as a sparse matrix
    Eigen::SimplicialLLT<SparseType, Eigen::Lower>  lltSolver;  // The solver for the diagonal band preconditioner
    // Eigen::SimplicialLDLT<SparseType, Eigen::Lower> ldltSolver; // The solver for the diagonal band preconditioner
    // Eigen::SimplicialLDLT<SparseType, Eigen::Lower> bandSolver; // The solver for the diagonal band preconditioner
    // Eigen::SparseLU<SparseType, Eigen::COLAMDOrdering<int>> qrSolver; // The solver for the diagonal band preconditioner
    // Eigen::ConjugateGradient<SparseType, Eigen::Lower | Eigen::Upper> cgSolver; // The solver for the diagonal band preconditioner
    // Eigen::BiCGSTAB<SparseTypeRowM> bandSolver; // The solver for the diagonal band preconditioner
    SparseType sparseMatrix;
    VectorType solverGuess;
    T          get_matrix_element(long I, long J) const;
    VectorType get_diagonal_new(long offset) const;
    VectorType get_diagonal_old(long offset) const;
    VectorType get_diagonal(long offset) const;
    void       thomas(long rows, T *x, T *const dl, T *const dm, T *const du);
    void       thomas2(long rows, T *x, T *const dl, T *const dm, T *const du);
    // void                             thomas(const long rows, const VectorType &x, const VectorType &dl, const VectorType &dm, const VectorType &du);

    // Shift stuff
    std::complex<double>  sigma         = cplx(0.0, 0.0); // The shift
    bool                  readyShift    = false;          // Flag to make sure the shift has occurred
    constexpr static bool readyFactorOp = false;          // Flag to check if factorization has occurred
    bool                  readyCalcPc   = false;
    long                  pcBandwidth   = 4l;

    public:
    MatVecMPOS() = default;
    template<typename EnvType>
    MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos, /*!< The Hamiltonian MPO's  */
               const env_pair<const EnvType &>                          &envs  /*!< The left and right environments.  */
    );
    // Functions used in Arpack++ solver
    [[nodiscard]] int rows() const; /*!< Linear size\f$d^2 \times \chi_L \times \chi_R \f$  */
    [[nodiscard]] int cols() const; /*!< Linear size\f$d^2 \times \chi_L \times \chi_R \f$  */

    void FactorOP();                                    //  Factorizes (A-sigma*I) (or finds its diagonal elements)
    void MultOPv(T *mps_in_, T *mps_out, T eval = 0.0); //  Applies the preconditioner as the matrix-vector product x_out <- inv(A-sigma*I)*x_in.
    void MultOPv(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *err); //  Applies the preconditioner
    void MultAx(T *mps_in_, T *mps_out_); //  Computes the matrix-vector multiplication x_out <- A*x_in.
    void MultAx(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *err);

    void CalcPc(T shift = 0.0);                         //  Calculates the diagonal or tridiagonal part of A
    void MultPc(T *mps_in_, T *mps_out, T shift = 0.0); //  Applies the preconditioner as the matrix-vector product x_out <- inv(A-sigma*I)*x_in.
    void MultPc(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *err); //  Applies the preconditioner

    // Various utility functions
    long num_mv = 0;
    long num_op = 0;
    long num_pc = 0;
    void print() const;
    void reset();
    void set_shift(std::complex<double> shift);
    void set_mode(eig::Form form_);
    void set_side(eig::Side side_);
    void set_lltBandwidth(long bandwidth); // the llt preconditioner bandwidth (default 8) (tridiagonal has bandwidth == 1)

    [[nodiscard]] T                                       get_shift() const;
    [[nodiscard]] eig::Form                               get_form() const;
    [[nodiscard]] eig::Side                               get_side() const;
    [[nodiscard]] eig::Type                               get_type() const;
    [[nodiscard]] const std::vector<Eigen::Tensor<T, 4>> &get_mpos() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3>         &get_envL() const;
    [[nodiscard]] const Eigen::Tensor<Scalar, 3>         &get_envR() const;
    [[nodiscard]] long                                    get_size() const;
    [[nodiscard]] std::array<long, 3>                     get_shape_mps() const;
    [[nodiscard]] std::vector<std::array<long, 4>>        get_shape_mpo() const;
    [[nodiscard]] std::array<long, 3>                     get_shape_envL() const;
    [[nodiscard]] std::array<long, 3>                     get_shape_envR() const;
    [[nodiscard]] Eigen::Tensor<Scalar, 6>                get_tensor() const;
    [[nodiscard]] MatrixType                              get_matrix() const;
    [[nodiscard]] SparseType                              get_sparse_matrix() const;
    [[nodiscard]] double                                  get_sparsity() const;
    [[nodiscard]] long                                    get_non_zeros() const;

    [[nodiscard]] bool isReadyShift() const;
    [[nodiscard]] bool isReadyFactorOp() const;

    // Timers
    std::unique_ptr<tid::ur> t_factorOP; // Factorization time
    std::unique_ptr<tid::ur> t_genMat;
    std::unique_ptr<tid::ur> t_multOPv;
    std::unique_ptr<tid::ur> t_multPc; // Preconditioner time
    std::unique_ptr<tid::ur> t_multAx; // Matvec time
};
