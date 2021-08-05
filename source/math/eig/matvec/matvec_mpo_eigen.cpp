#include "matvec_mpo_eigen.h"
#include "../log.h"
#include <math/svd.h>
#include <tid/tid.h>
#include <tools/common/contraction.h>
#include <math/tenx.h>

template<typename Scalar_>
template<typename T>
MatVecMPOEigen<Scalar_>::MatVecMPOEigen(const Eigen::Tensor<T,3> & envL_,      /*!< The left block tensor.  */
                                        const Eigen::Tensor<T,3> & envR_,      /*!< The right block tensor.  */
                                        const Eigen::Tensor<T,4> & mpo_       /*!< The Hamiltonian MPO's  */
                              ) {
    if constexpr(std::is_same_v<Scalar_, T>){
        mpo = mpo_;
        envL = envL_;
        envR = envR_;
    }else if constexpr(std::is_same_v<Scalar_, eig::real> and std::is_same_v<T, eig::cplx>){
        // This should only be done if we know for a fact that there is no imaginary component.
        mpo = mpo_.real();
        envL = envL_.real();
        envR = envR_.real();
    }else if constexpr(std::is_same_v<Scalar_, eig::cplx> and std::is_same_v<T, eig::real>){
        mpo = mpo_.template cast<eig::cplx>();
        envL = envL_.template cast<eig::cplx>();
        envR = envR_.template cast<eig::cplx>();
    }


    shape_mps  = {mpo.dimension(2), envL.dimension(0), envR.dimension(0)};
    mps_size   = shape_mps[0] * shape_mps[1] * shape_mps[2];
    t_factorOP = std::make_unique<tid::ur>("Time FactorOp");
    t_multOPv  = std::make_unique<tid::ur>("Time MultOpv");
    t_multAx   = std::make_unique<tid::ur>("Time MultAx");
}

template MatVecMPOEigen<eig::cplx>::MatVecMPOEigen(const Eigen::Tensor<eig::real,3> & envL_,const Eigen::Tensor<eig::real,3> & envR_,const Eigen::Tensor<eig::real,4> &  mpo_);
template MatVecMPOEigen<eig::real>::MatVecMPOEigen(const Eigen::Tensor<eig::real,3> & envL_,const Eigen::Tensor<eig::real,3> & envR_,const Eigen::Tensor<eig::real,4> &  mpo_);
template MatVecMPOEigen<eig::cplx>::MatVecMPOEigen(const Eigen::Tensor<eig::cplx,3> & envL_,const Eigen::Tensor<eig::cplx,3> & envR_,const Eigen::Tensor<eig::cplx,4> &  mpo_);
template MatVecMPOEigen<eig::real>::MatVecMPOEigen(const Eigen::Tensor<eig::cplx,3> & envL_,const Eigen::Tensor<eig::cplx,3> & envR_,const Eigen::Tensor<eig::cplx,4> &  mpo_);


template<typename Scalar>
void MatVecMPOEigen<Scalar>::FactorOP()
/* We don't actually invert a matrix here: we let an iterative matrix-free solver apply OP^-1 x */
{
    if(readyFactorOp) return; // happens only once
    if(not readyShift) throw std::runtime_error("Cannot FactorOP: Shift value sigma has not been set.");
    t_factorOP->tic();
    readyFactorOp = true;
    t_factorOP->toc();
}

template<typename T>
void MatVecMPOEigen<T>::MultAx(T *mps_in_, T *mps_out_) {
    auto token = t_multAx->tic_token();
    Eigen::TensorMap<Eigen::Tensor<Scalar,3>> mps_in(mps_in_, shape_mps);
    Eigen::TensorMap<Eigen::Tensor<Scalar,3>> mps_out(mps_out_, shape_mps);
    tools::common::contraction::matrix_vector_product(mps_out, mps_in, mpo, envL, envR);
    counter++;
}

template<typename T>
void MatVecMPOEigen<T>::MultAx(void *x, int *ldx, void *y, int *ldy, int *blockSize, [[maybe_unused]] primme_params *primme, [[maybe_unused]] int *err) {
    for(int i = 0; i < *blockSize; i++) {
        T *mps_in  = static_cast<T *>(x) + *ldx * i;
        T *mps_out = static_cast<T *>(y) + *ldy * i;
        MultAx(mps_in, mps_out);
        counter++;
    }
    *err = 0;
}

template<typename T>
void MatVecMPOEigen<T>::MultOPv(T *mps_in_, T *mps_out_) {
    t_multOPv->tic();
    Eigen::TensorMap<Eigen::Tensor<Scalar,3>> mps_in(mps_in_, shape_mps);
    Eigen::TensorMap<Eigen::Tensor<Scalar,3>> mps_out(mps_out_, shape_mps);
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
    t_multOPv->toc();
    counter++;
}

template<typename Scalar>
void MatVecMPOEigen<Scalar>::print() const {}

template<typename T>
void MatVecMPOEigen<T>::set_shift(std::complex<double> sigma_) {
    if(readyShift) return; // This only happens once!!
    if(readyCompress) throw std::runtime_error("Cannot shift the matrix: it is already compressed!");
    sigma = sigma_;

    // The MPO is a rank4 tensor ijkl where the first 2 ij indices draw a simple
    // rank2 matrix, where each element is also a matrix with the size
    // determined by the last 2 indices kl.
    // When we shift an MPO, all we do is subtract a diagonal matrix from
    // the botton left corner of the ij-matrix.
    auto dims = mpo.dimensions();
    if(dims[2] != dims[3])throw std::logic_error(fmt::format("MPO has different spin dimensions up and down: {}", dims));
    auto spindim = dims[2];

    // Setup extents and handy objects
    std::array<long, 4>                              offset4{dims[0] - 1, 0, 0, 0};
    std::array<long, 4>                              extent4{1, 1, spindim, spindim};
    std::array<long, 2>                              extent2{spindim, spindim};

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> sigma_Id;
    auto sigma_ID = tenx::TensorIdentity<Scalar>(spindim);

    if constexpr(std::is_same_v<T, eig::real>)
        mpo.slice(offset4, extent4).reshape(extent2) -=  tenx::TensorIdentity<Scalar>(spindim) * std::real(sigma);
    else
        mpo.slice(offset4, extent4).reshape(extent2) -=  tenx::TensorIdentity<Scalar>(spindim) * sigma;
    readyShift = true;
}


template<typename T>
void MatVecMPOEigen<T>::compress() {
    if(readyCompress) return;

    svd::settings svd_settings;
    svd_settings.use_lapacke = true;
    svd_settings.use_bdc     = false;
    svd_settings.threshold   = 1e-12;
    svd_settings.switchsize  = 4096;
    svd::solver svd(svd_settings);

    Eigen::Tensor<T, 4> mpo_tmp = mpo;
    Eigen::Tensor<T, 4> mpo_l2r, mpo_r2l;
    {
        // Compress left to right
        auto [U_l2r, S, V_l2r] = svd.split_mpo_l2r(mpo_tmp);
        mpo_l2r                = U_l2r.contract(tenx::asDiagonal(S), tenx::idx({1}, {0})).shuffle(std::array<long, 4>{0, 3, 1, 2});

        // Contract V_l2r into the right environment
        Eigen::Tensor<Scalar,3> envR_tmp = envR.contract(V_l2r, tenx::idx({2}, {1}));
        envR = envR_tmp;
    }
    {
        // Compress right to left
        auto [U_r2l, S, V_r2l] = svd.split_mpo_r2l(mpo_l2r);
        mpo_r2l                = tenx::asDiagonal(S).contract(V_r2l, tenx::idx({1}, {0}));

        // Contract U_r2l into the left environment
        Eigen::Tensor<Scalar,3> envL_tmp = envL.contract(U_r2l, tenx::idx({2}, {0}));
        envL = envL_tmp;
    }
    {
        // The new mpo is what remains
        mpo = mpo_r2l;
    }
    readyCompress = true;
    eig::log->trace("Compressed MPO² dimensions {}", mpo.dimensions());
}

template<typename Scalar>
void MatVecMPOEigen<Scalar>::set_mode(const eig::Form form_) {
    form = form_;
}
template<typename Scalar>
void MatVecMPOEigen<Scalar>::set_side(const eig::Side side_) {
    side = side_;
}

template<typename Scalar>
eig::Form MatVecMPOEigen<Scalar>::get_form() const {
    return form;
}
template<typename Scalar>
eig::Side MatVecMPOEigen<Scalar>::get_side() const {
    return side;
}
template<typename Scalar>
eig::Type MatVecMPOEigen<Scalar>::get_type() const {
    if constexpr(std::is_same_v<Scalar, eig::real>)
    return eig::Type::REAL;
    else if constexpr(std::is_same_v<Scalar, eig::cplx>)
    return eig::Type::CPLX;
    else
        throw std::runtime_error("Unsupported type");
}

template<typename Scalar> const Eigen::Tensor<Scalar,4> & MatVecMPOEigen<Scalar>::get_mpo() const { return mpo; }
template<typename Scalar> const Eigen::Tensor<Scalar,3> & MatVecMPOEigen<Scalar>::get_envL() const { return envL; }
template<typename Scalar> const Eigen::Tensor<Scalar,3> & MatVecMPOEigen<Scalar>::get_envR() const { return envR; }
template<typename Scalar> std::array<long, 4> MatVecMPOEigen<Scalar>::get_shape_mpo() const { return mpo.dimensions(); }
template<typename Scalar> std::array<long, 3> MatVecMPOEigen<Scalar>::get_shape_envL() const { return envL.dimensions(); }
template<typename Scalar> std::array<long, 3> MatVecMPOEigen<Scalar>::get_shape_envR() const { return envR.dimensions(); }

template<typename Scalar> bool MatVecMPOEigen<Scalar>::isReadyFactorOp() const { return readyFactorOp; }
template<typename Scalar> bool MatVecMPOEigen<Scalar>::isReadyShift() const { return readyShift; }
template<typename Scalar> bool MatVecMPOEigen<Scalar>::isReadyCompress() const { return readyCompress; }


// Explicit instantiations
template class MatVecMPOEigen<double>;
    template class MatVecMPOEigen<std::complex<double>>;
