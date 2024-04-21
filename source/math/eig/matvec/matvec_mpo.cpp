#include "matvec_mpo.h"
#include "../log.h"
#include "math/eig/solver.h"
#include "math/svd.h"
#include "math/tenx.h"
#include "tid/tid.h"
#include "tools/common/contraction.h"
#include <Eigen/Cholesky>
#include <primme/primme.h>

namespace eig {

#ifdef NDEBUG
    inline constexpr bool debug = false;
#else
    inline constexpr bool debug = true;
#endif
}

template<typename T>
template<typename S>
MatVecMPO<T>::MatVecMPO(const Eigen::Tensor<S, 3> &envL_, /*!< The left block tensor.  */
                        const Eigen::Tensor<S, 3> &envR_, /*!< The right block tensor.  */
                        const Eigen::Tensor<S, 4> &mpo_   /*!< The Hamiltonian MPO's  */
) {
    static_assert(std::is_same_v<S, real> or std::is_same_v<S, cplx>);

    if constexpr(std::is_same_v<T, S>) {
        mpo  = mpo_;
        envL = envL_;
        envR = envR_;
    } else if constexpr(std::is_same_v<T, real> and std::is_same_v<S, cplx>) {
        // This should only be done if we know for a fact that there is no imaginary component.
        if constexpr(eig::debug) {
            if(not tenx::isReal(mpo_)) throw except::runtime_error("mpo is not real");
            if(not tenx::isReal(envL_)) throw except::runtime_error("envL is not real");
            if(not tenx::isReal(envR_)) throw except::runtime_error("envR is not real");
        }
        mpo  = mpo_.real();
        envL = envL_.real();
        envR = envR_.real();
    } else if constexpr(std::is_same_v<T, cplx> and std::is_same_v<S, real>) {
        mpo  = mpo_.template cast<cplx>();
        envL = envL_.template cast<cplx>();
        envR = envR_.template cast<cplx>();
    }

    shape_mps  = {mpo.dimension(2), envL.dimension(0), envR.dimension(0)};
    size_mps   = shape_mps[0] * shape_mps[1] * shape_mps[2];
    t_factorOP = std::make_unique<tid::ur>("Time FactorOp");
    t_genMat   = std::make_unique<tid::ur>("Time genMat");
    t_multOPv  = std::make_unique<tid::ur>("Time MultOpv");
    t_multAx   = std::make_unique<tid::ur>("Time MultAx");
}

template MatVecMPO<cplx>::MatVecMPO(const Eigen::Tensor<real, 3> &envL_, const Eigen::Tensor<real, 3> &envR_, const Eigen::Tensor<real, 4> &mpo_);
template MatVecMPO<real>::MatVecMPO(const Eigen::Tensor<real, 3> &envL_, const Eigen::Tensor<real, 3> &envR_, const Eigen::Tensor<real, 4> &mpo_);
template MatVecMPO<cplx>::MatVecMPO(const Eigen::Tensor<cplx, 3> &envL_, const Eigen::Tensor<cplx, 3> &envR_, const Eigen::Tensor<cplx, 4> &mpo_);
template MatVecMPO<real>::MatVecMPO(const Eigen::Tensor<cplx, 3> &envL_, const Eigen::Tensor<cplx, 3> &envR_, const Eigen::Tensor<cplx, 4> &mpo_);

template<typename T>
int MatVecMPO<T>::rows() const {
    return safe_cast<int>(size_mps);
}

template<typename T>
int MatVecMPO<T>::cols() const {
    return safe_cast<int>(size_mps);
}

template<typename T>
void MatVecMPO<T>::FactorOP() {
    auto t_token = t_factorOP->tic_token();
    if(readyFactorOp) { return; }
    if(factorization == eig::Factorization::NONE) {
        readyFactorOp = true;
        return;
    }
    MatrixType A_matrix = get_matrix();
    if(not readyShift and std::abs(get_shift()) != 0.0) { A_matrix.diagonal() -= VectorType::Constant(rows(), get_shift()); }

    if(factorization == eig::Factorization::LDLT) {
        eig::log->debug("LDLT Factorization");
        ldlt.compute(A_matrix);
    } else if(factorization == eig::Factorization::LLT) {
        eig::log->debug("LLT Factorization");
        llt.compute(A_matrix);
    } else if(factorization == eig::Factorization::LU) {
        eig::log->debug("LU Factorization");
        lu.compute(A_matrix);
    } else if(factorization == eig::Factorization::NONE) {
        /* We don't actually invert a matrix here: we let an iterative matrix-free solver apply OP^-1 x */
        if(not readyShift) throw std::runtime_error("Cannot FactorOP with Factorization::NONE: Shift value sigma has not been set on the MPO.");
    }
    eig::log->debug("Finished factorization");
    readyFactorOp = true;
}

template<typename T>
void MatVecMPO<T>::MultAx(T *mps_in_, T *mps_out_) {
    auto                                  token = t_multAx->tic_token();
    Eigen::TensorMap<Eigen::Tensor<T, 3>> mps_in(mps_in_, shape_mps);
    Eigen::TensorMap<Eigen::Tensor<T, 3>> mps_out(mps_out_, shape_mps);
    tools::common::contraction::matrix_vector_product(mps_out, mps_in, mpo, envL, envR);
    num_mv++;
}

template<typename T>
void MatVecMPO<T>::MultAx(T *mps_in, T *mps_out, T *mpo_ptr, T *envL_ptr, T *envR_ptr, std::array<long, 3> shape_mps_, std::array<long, 4> shape_mpo_) {
    auto                token       = t_multAx->tic_token();
    std::array<long, 3> shape_envL_ = {shape_mps_[1], shape_mps_[1], shape_mpo_[0]};
    std::array<long, 3> shape_envR_ = {shape_mps_[2], shape_mps_[2], shape_mpo_[1]};
    tools::common::contraction::matrix_vector_product(mps_out, mps_in, shape_mps, mpo_ptr, shape_mpo_, envL_ptr, shape_envL_, envR_ptr, shape_envR_);
    num_mv++;
}

template<typename T>
void MatVecMPO<T>::MultAx(void *x, int *ldx, void *y, int *ldy, int *blockSize, [[maybe_unused]] primme_params *primme, [[maybe_unused]] int *err) {
    auto token = t_multAx->tic_token();
    for(int i = 0; i < *blockSize; i++) {
        T                                    *mps_in_ptr  = static_cast<T *>(x) + *ldx * i;
        T                                    *mps_out_ptr = static_cast<T *>(y) + *ldy * i;
        Eigen::TensorMap<Eigen::Tensor<T, 3>> mps_in(mps_in_ptr, shape_mps);
        Eigen::TensorMap<Eigen::Tensor<T, 3>> mps_out(mps_out_ptr, shape_mps);
        tools::common::contraction::matrix_vector_product(mps_out, mps_in, mpo, envL, envR);
    }

    num_mv += *blockSize;
    *err = 0;
}

template<typename T>
void MatVecMPO<T>::MultOPv(T *mps_in_, T *mps_out_) {
    auto                                  token = t_multOPv->tic_token();
    Eigen::TensorMap<Eigen::Tensor<T, 3>> mps_in(mps_in_, shape_mps);
    Eigen::TensorMap<Eigen::Tensor<T, 3>> mps_out(mps_out_, shape_mps);
    switch(side) {
        case eig::Side::R: {
            tools::common::contraction::matrix_inverse_vector_product(mps_out, mps_in, mpo, envL, envR);
            break;
        }
        case eig::Side::L: {
            throw std::runtime_error("Left sided matrix-free MultOPv has not been implemented");
            break;
        }
        case eig::Side::LR: {
            throw std::runtime_error("eigs cannot handle sides L and R simultaneously");
        }
    }
    num_op++;
}

template<typename T>
void MatVecMPO<T>::MultOPv(void *x, int *ldx, void *y, int *ldy, int *blockSize, [[maybe_unused]] primme_params *primme, int *err) {
    auto token = t_multOPv->tic_token();
    switch(side) {
        case eig::Side::R: {
            for(int i = 0; i < *blockSize; i++) {
                T *x_ptr = static_cast<T *>(x) + *ldx * i;
                T *y_ptr = static_cast<T *>(y) + *ldy * i;
                if(factorization == eig::Factorization::NONE) {
                    Eigen::TensorMap<Eigen::Tensor<T, 3>> x_map(x_ptr, shape_mps);
                    Eigen::TensorMap<Eigen::Tensor<T, 3>> y_map(y_ptr, shape_mps);
                    tools::common::contraction::matrix_inverse_vector_product(y_map, x_map, mpo, envL, envR);
                } else if(factorization == eig::Factorization::LDLT) {
                    Eigen::Map<VectorType> x_map(x_ptr, *ldx);
                    Eigen::Map<VectorType> y_map(y_ptr, *ldy);
                    y_map.noalias() = ldlt.solve(x_map);
                } else if(factorization == eig::Factorization::LLT) {
                    Eigen::Map<VectorType> x_map(x_ptr, *ldx);
                    Eigen::Map<VectorType> y_map(y_ptr, *ldy);
                    y_map.noalias() = llt.solve(x_map);
                } else if(factorization == eig::Factorization::LU) {
                    Eigen::Map<VectorType> x_map(x_ptr, *ldx);
                    Eigen::Map<VectorType> y_map(y_ptr, *ldy);
                    y_map.noalias() = lu.solve(x_map);
                }
                num_op++;
            }
            break;
        }
        case eig::Side::L: {
            throw std::runtime_error("Left sided matrix-free MultOPv has not been implemented");
            break;
        }
        case eig::Side::LR: {
            throw std::runtime_error("eigs cannot handle sides L and R simultaneously");
        }
    }
    *err = 0;
}

template<typename T>
void MatVecMPO<T>::print() const {}

template<typename T>
void MatVecMPO<T>::reset() {
    if(t_factorOP) t_factorOP->reset();
    if(t_multOPv) t_multOPv->reset();
    if(t_genMat) t_genMat->reset();
    if(t_multAx) t_multAx->reset();
    num_mv = 0;
    num_op = 0;
}

template<typename T>
void MatVecMPO<T>::set_shift(std::complex<double> shift) {
    // Here we set an energy shift directly on the MPO.
    // This only works if the MPO is not compressed already.
    if(readyShift) return;
    if(sigma == shift) return;

            eig::log->trace("Setting shift = {:.16f} + i{:.16f}", std::real(shift), std::imag(shift));
            sigma = shift; // We can shift the diagonal of the full matrix instead
            return;



    // The MPO is a rank4 tensor ijkl where the first 2 ij indices draw a simple
    // rank2 matrix, where each element is also a matrix with the size
    // determined by the last 2 indices kl.
    // When we shift an MPO, all we do is subtract a diagonal matrix from
    // the botton left corner of the ij-matrix.
    auto dims = mpo.dimensions();
    if(dims[2] != dims[3]) throw except::logic_error("MPO has different spin dimensions up and down: {}", dims);
    auto spindim = dims[2];
    long offset1 = dims[0] - 1;

    // Setup extents and handy objects
    std::array<long, 4> offset4{offset1, 0, 0, 0};
    std::array<long, 4> extent4{1, 1, spindim, spindim};
    std::array<long, 2> extent2{spindim, spindim};
    auto                id = tenx::TensorIdentity<T>(spindim);
    // We undo the previous sigma and then subtract the new one. We are aiming for [A - I*shift]
    if constexpr(std::is_same_v<T, real>)
        mpo.slice(offset4, extent4).reshape(extent2) += id * std::real(sigma - shift);
    else
        mpo.slice(offset4, extent4).reshape(extent2) += id * (sigma - shift);
    sigma = shift;
    eig::log->debug("Shifted MPO dimensions {}", mpo.dimensions());

    readyShift = true;
}

template<typename T>
void MatVecMPO<T>::set_mode(const eig::Form form_) {
    form = form_;
}
template<typename T>
void MatVecMPO<T>::set_side(const eig::Side side_) {
    side = side_;
}

template<typename T>
T MatVecMPO<T>::get_shift() const {
    if constexpr(std::is_same_v<T, real>)
        return std::real(sigma);
    else
        return sigma;
}

template<typename T>
eig::Form MatVecMPO<T>::get_form() const {
    return form;
}
template<typename T>
eig::Side MatVecMPO<T>::get_side() const {
    return side;
}
template<typename T>
eig::Type MatVecMPO<T>::get_type() const {
    if constexpr(std::is_same_v<T, real>)
        return eig::Type::REAL;
    else if constexpr(std::is_same_v<T, cplx>)
        return eig::Type::CPLX;
    else
        throw std::runtime_error("Unsupported type");
}

template<typename T>
const Eigen::Tensor<T, 4> &MatVecMPO<T>::get_mpo() const {
    return mpo;
}
template<typename T>
const Eigen::Tensor<T, 3> &MatVecMPO<T>::get_envL() const {
    return envL;
}
template<typename T>
const Eigen::Tensor<T, 3> &MatVecMPO<T>::get_envR() const {
    return envR;
}

template<typename T>
long MatVecMPO<T>::get_size() const {
    return size_mps;
}

template<typename T>
std::array<long, 3> MatVecMPO<T>::get_shape_mps() const {
    return shape_mps;
}
template<typename T>
std::array<long, 4> MatVecMPO<T>::get_shape_mpo() const {
    return mpo.dimensions();
}
template<typename T>
std::array<long, 3> MatVecMPO<T>::get_shape_envL() const {
    return envL.dimensions();
}
template<typename T>
std::array<long, 3> MatVecMPO<T>::get_shape_envR() const {
    return envR.dimensions();
}
template<typename T>
Eigen::Tensor<T, 6> MatVecMPO<T>::get_tensor() const {
    auto t_token = t_genMat->tic_token();
    eig::log->debug("Generating tensor");

    auto d0      = shape_mps[0];
    auto d1      = shape_mps[1];
    auto d2      = shape_mps[2];
    auto &threads = tenx::threads::get();
    Eigen::Tensor<T, 6> tensor;
    tensor.resize(tenx::array6{d0, d1, d2, d0, d1, d2});
    tensor.device(*threads->dev) =
        get_envL().contract(get_mpo(), tenx::idx({2}, {0})).contract(get_envR(), tenx::idx({2}, {2})).shuffle(tenx::array6{2, 0, 4, 3, 1, 5});
    return tensor;
}

template<typename T>
typename MatVecMPO<T>::MatrixType MatVecMPO<T>::get_matrix() const {
    return tenx::MatrixCast(get_tensor(), rows(), cols());
}

template<typename T>
bool MatVecMPO<T>::isReadyFactorOp() const {
    return readyFactorOp;
}
template<typename T>
bool MatVecMPO<T>::isReadyShift() const {
    return readyShift;
}

// Explicit instantiations
template class MatVecMPO<double>;
template class MatVecMPO<std::complex<double>>;
