#include "matvec_mpo.h"
#include "../log.h"
#include "math/eig/solver.h"
#include "math/svd.h"
#include "math/tenx.h"
#include "tid/tid.h"
#include "tools/common/contraction.h"
#include <primme/primme.h>

namespace eig {

#ifdef NDEBUG
    inline constexpr bool debug = false;
#else
    inline constexpr bool debug = true;
#endif
}

// template<typename Scalar_>
// MatVecMPO<Scalar_>::~MatVecMPO() noexcept = default;

template<typename Scalar_>
template<typename T>
MatVecMPO<Scalar_>::MatVecMPO(const Eigen::Tensor<T, 3> &envL_, /*!< The left block tensor.  */
                              const Eigen::Tensor<T, 3> &envR_, /*!< The right block tensor.  */
                              const Eigen::Tensor<T, 4> &mpo_   /*!< The Hamiltonian MPO's  */
) {
    if constexpr(std::is_same_v<Scalar_, T>) {
        mpo  = mpo_;
        envL = envL_;
        envR = envR_;
    } else if constexpr(std::is_same_v<Scalar_, eig::real> and std::is_same_v<T, eig::cplx>) {
        // This should only be done if we know for a fact that there is no imaginary component.
        if constexpr(eig::debug) {
            if(not tenx::isReal(mpo_)) throw std::runtime_error("mpo is not real");
            if(not tenx::isReal(envL_)) throw std::runtime_error("envL is not real");
            if(not tenx::isReal(envR_)) throw std::runtime_error("envR is not real");
        }
        mpo  = mpo_.real();
        envL = envL_.real();
        envR = envR_.real();
    } else if constexpr(std::is_same_v<Scalar_, eig::cplx> and std::is_same_v<T, eig::real>) {
        mpo  = mpo_.template cast<eig::cplx>();
        envL = envL_.template cast<eig::cplx>();
        envR = envR_.template cast<eig::cplx>();
    }

    shape_mps  = {mpo.dimension(2), envL.dimension(0), envR.dimension(0)};
    mps_size   = shape_mps[0] * shape_mps[1] * shape_mps[2];
    t_factorOP = std::make_unique<tid::ur>("Time FactorOp");
    t_multOPv  = std::make_unique<tid::ur>("Time MultOpv");
    t_multAx   = std::make_unique<tid::ur>("Time MultAx");
}

template MatVecMPO<eig::cplx>::MatVecMPO(const Eigen::Tensor<eig::real, 3> &envL_, const Eigen::Tensor<eig::real, 3> &envR_,
                                         const Eigen::Tensor<eig::real, 4> &mpo_);
template MatVecMPO<eig::real>::MatVecMPO(const Eigen::Tensor<eig::real, 3> &envL_, const Eigen::Tensor<eig::real, 3> &envR_,
                                         const Eigen::Tensor<eig::real, 4> &mpo_);
template MatVecMPO<eig::cplx>::MatVecMPO(const Eigen::Tensor<eig::cplx, 3> &envL_, const Eigen::Tensor<eig::cplx, 3> &envR_,
                                         const Eigen::Tensor<eig::cplx, 4> &mpo_);
template MatVecMPO<eig::real>::MatVecMPO(const Eigen::Tensor<eig::cplx, 3> &envL_, const Eigen::Tensor<eig::cplx, 3> &envR_,
                                         const Eigen::Tensor<eig::cplx, 4> &mpo_);

template<typename Scalar>
void MatVecMPO<Scalar>::FactorOP()
/* We don't actually invert a matrix here: we let an iterative matrix-free solver apply OP^-1 x */
{
    if(readyFactorOp) return; // happens only once
    if(not readyShift) throw std::runtime_error("Cannot FactorOP: Shift value sigma has not been set.");
    t_factorOP->tic();
    readyFactorOp = true;
    t_factorOP->toc();
}

template<typename T>
void MatVecMPO<T>::MultAx(T *mps_in_, T *mps_out_) {
    auto                                       token = t_multAx->tic_token();
    Eigen::TensorMap<Eigen::Tensor<Scalar, 3>> mps_in(mps_in_, shape_mps);
    Eigen::TensorMap<Eigen::Tensor<Scalar, 3>> mps_out(mps_out_, shape_mps);
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
        T                                         *mps_in_ptr  = static_cast<T *>(x) + *ldx * i;
        T                                         *mps_out_ptr = static_cast<T *>(y) + *ldy * i;
        Eigen::TensorMap<Eigen::Tensor<Scalar, 3>> mps_in(mps_in_ptr, shape_mps);
        Eigen::TensorMap<Eigen::Tensor<Scalar, 3>> mps_out(mps_out_ptr, shape_mps);
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
                T                                    *mps_in_  = static_cast<T *>(x) + *ldx * i;
                T                                    *mps_out_ = static_cast<T *>(y) + *ldy * i;
                Eigen::TensorMap<Eigen::Tensor<T, 3>> mps_in(mps_in_, shape_mps);
                Eigen::TensorMap<Eigen::Tensor<T, 3>> mps_out(mps_out_, shape_mps);
                tools::common::contraction::matrix_inverse_vector_product(mps_out, mps_in, mpo, envL, envR);
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

template<typename Scalar>
void MatVecMPO<Scalar>::print() const {}

template<typename Scalar>
void MatVecMPO<Scalar>::reset() {
    if(t_factorOP) t_factorOP->reset();
    if(t_multOPv) t_multOPv->reset();
    if(t_multAx) t_multAx->reset();
    num_mv = 0;
    num_op = 0;
}

template<typename T>
void MatVecMPO<T>::set_shift(std::complex<double> shift) {
    // Here we set an energy shift directly on the MPO.
    // This only works if the MPO is not compressed already.
    readyShift = sigma == shift;

    if(readyShift) return; // This only happens once!!
    if(readyCompress) throw std::runtime_error("Cannot shift the matrix: it is already compressed!");

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

    auto id = tenx::TensorIdentity<Scalar>(spindim);
    // We undo the previous sigma and then subtract the new one. We are aiming for [A - I*shift]
    if constexpr(std::is_same_v<T, eig::real>)
        mpo.slice(offset4, extent4).reshape(extent2) += id * std::real(sigma - shift);
    else
        mpo.slice(offset4, extent4).reshape(extent2) += id * (sigma - shift);

    sigma = shift;
    eig::log->debug("Shifted MPO dimensions {}", mpo.dimensions());

    readyShift = true;
}

template<typename T>
void MatVecMPO<T>::compress() {
    if(readyCompress) return;
    svd::settings svd_settings;
    svd_settings.svd_lib        = SVDLib::lapacke;
    svd_settings.use_bdc        = false;
    svd_settings.threshold      = 1e-12;
    svd_settings.threshold_tr   = 1e-4;
    svd_settings.switchsize_bdc = 4096;
    svd::solver svd(svd_settings);

    Eigen::Tensor<T, 4> mpo_tmp = mpo;
    Eigen::Tensor<T, 4> mpo_l2r, mpo_r2l;
    {
        // Compress left to right
        auto [U_l2r, S, V_l2r] = svd.split_mpo_l2r(mpo_tmp);
        mpo_l2r                = U_l2r.contract(tenx::asDiagonal(S), tenx::idx({1}, {0})).shuffle(std::array<long, 4>{0, 3, 1, 2});

        // Contract V_l2r into the right environment
        Eigen::Tensor<Scalar, 3> envR_tmp = envR.contract(V_l2r, tenx::idx({2}, {1}));
        envR                              = envR_tmp;
    }
    {
        // Compress right to left
        auto [U_r2l, S, V_r2l] = svd.split_mpo_r2l(mpo_l2r);
        mpo_r2l                = tenx::asDiagonal(S).contract(V_r2l, tenx::idx({1}, {0}));

        // Contract U_r2l into the left environment
        Eigen::Tensor<Scalar, 3> envL_tmp = envL.contract(U_r2l, tenx::idx({2}, {0}));
        envL                              = envL_tmp;
    }
    {
        // The new mpo is what remains
        mpo = mpo_r2l;
    }
    readyCompress = true;
    eig::log->debug("Compressed MPO dimensions {} -> {}", mpo.dimensions(), mpo_tmp.dimensions());
}

template<typename Scalar>
void MatVecMPO<Scalar>::set_mode(const eig::Form form_) {
    form = form_;
}
template<typename Scalar>
void MatVecMPO<Scalar>::set_side(const eig::Side side_) {
    side = side_;
}

template<typename Scalar>
template<typename T>
T MatVecMPO<Scalar>::get_shift() const {
    if constexpr(std::is_same_v<T, eig::cplx>)
        return sigma;
    else if constexpr(std::is_floating_point_v<T>)
        return std::real(sigma);
    else
        return static_cast<T>(sigma);
}

template eig::real MatVecMPO<eig::real>::get_shift<eig::real>() const;
template eig::real MatVecMPO<eig::cplx>::get_shift<eig::real>() const;
template eig::cplx MatVecMPO<eig::real>::get_shift<eig::cplx>() const;
template eig::cplx MatVecMPO<eig::cplx>::get_shift<eig::cplx>() const;

template<typename Scalar>
eig::Form MatVecMPO<Scalar>::get_form() const {
    return form;
}
template<typename Scalar>
eig::Side MatVecMPO<Scalar>::get_side() const {
    return side;
}
template<typename Scalar>
eig::Type MatVecMPO<Scalar>::get_type() const {
    if constexpr(std::is_same_v<Scalar, eig::real>)
        return eig::Type::REAL;
    else if constexpr(std::is_same_v<Scalar, eig::cplx>)
        return eig::Type::CPLX;
    else
        throw std::runtime_error("Unsupported type");
}

template<typename Scalar>
const Eigen::Tensor<Scalar, 4> &MatVecMPO<Scalar>::get_mpo() const {
    return mpo;
}
template<typename Scalar>
const Eigen::Tensor<Scalar, 3> &MatVecMPO<Scalar>::get_envL() const {
    return envL;
}
template<typename Scalar>
const Eigen::Tensor<Scalar, 3> &MatVecMPO<Scalar>::get_envR() const {
    return envR;
}
template<typename Scalar>
std::array<long, 3> MatVecMPO<Scalar>::get_shape_mps() const {
    return shape_mps;
}
template<typename Scalar>
std::array<long, 4> MatVecMPO<Scalar>::get_shape_mpo() const {
    return mpo.dimensions();
}
template<typename Scalar>
std::array<long, 3> MatVecMPO<Scalar>::get_shape_envL() const {
    return envL.dimensions();
}
template<typename Scalar>
std::array<long, 3> MatVecMPO<Scalar>::get_shape_envR() const {
    return envR.dimensions();
}
template<typename Scalar>
Eigen::Tensor<Scalar, 6> MatVecMPO<Scalar>::get_tensor() const {
    auto                     d0 = shape_mps[0];
    auto                     d1 = shape_mps[1];
    auto                     d2 = shape_mps[2];
    Eigen::Tensor<Scalar, 6> tensor;
    tensor.resize(tenx::array6{d0, d1, d2, d0, d1, d2});
    tensor.device(tenx::omp::getDevice()) =
        get_envL().contract(get_mpo(), tenx::idx({2}, {0})).contract(get_envR(), tenx::idx({2}, {2})).shuffle(tenx::array6{2, 0, 4, 3, 1, 5});
    return tensor;
}

template<typename Scalar>
typename MatVecMPO<Scalar>::MatrixType MatVecMPO<Scalar>::get_matrix() const {
    return tenx::MatrixCast(get_tensor(), rows(), cols());
    //    return get_tensor().reshape(tenx::array2{rows(), cols()});
}

template<typename Scalar>
bool MatVecMPO<Scalar>::isReadyFactorOp() const {
    return readyFactorOp;
}
template<typename Scalar>
bool MatVecMPO<Scalar>::isReadyShift() const {
    return readyShift;
}
template<typename Scalar>
bool MatVecMPO<Scalar>::isReadyCompress() const {
    return readyCompress;
}

// Explicit instantiations
template class MatVecMPO<double>;
template class MatVecMPO<std::complex<double>>;
