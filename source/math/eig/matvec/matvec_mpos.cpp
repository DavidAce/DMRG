#include "matvec_mpos.h"
#include "../log.h"
#include "general/sfinae.h"
#include "math/eig/solver.h"
#include "math/svd.h"
#include "math/tenx.h"
#include "tensors/edges/EdgesFinite.h"
#include "tensors/site/env/EnvEne.h"
#include "tensors/site/env/EnvVar.h"
#include "tensors/site/mpo/MpoSite.h"
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
template<typename EnvType>
MatVecMPOS<T>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, /*!< The Hamiltonian MPO's  */
                          const env_pair<const EnvType &>                          &envs_  /*!< The left and right environments.  */
) {
    static_assert(sfinae::is_any_v<EnvType, EnvEne, EnvVar>);
    mpos.reserve(mpos_.size());
    if constexpr(std::is_same_v<T, cplx>) {
        for(const auto &mpo_ : mpos_) {
            if constexpr(std::is_same_v<EnvType, EnvEne>)
                mpos.emplace_back(mpo_.get().MPO());
            else
                mpos.emplace_back(mpo_.get().MPO2());
        }
        envL = envs_.L.get_block();
        envR = envs_.R.get_block();
    } else if constexpr(std::is_same_v<T, real>) {
        for(const auto &mpo_ : mpos_) {
            if(not tenx::isReal(mpo_.get().MPO())) throw except::runtime_error("mpo is not real");
            if constexpr(std::is_same_v<EnvType, EnvEne>)
                mpos.emplace_back(mpo_.get().MPO().real());
            else
                mpos.emplace_back(mpo_.get().MPO2().real());
        }
        if constexpr(eig::debug) {
            if(not tenx::isReal(envs_.L.get_block())) throw except::runtime_error("envL is not real");
            if(not tenx::isReal(envs_.R.get_block())) throw except::runtime_error("envR is not real");
        }
        envL = envs_.L.get_block().real();
        envR = envs_.R.get_block().real();
    }
    long spin_dim = 1;
    for(const auto &mpo : mpos) spin_dim *= mpo.dimension(2);


    shape_mps  = {spin_dim, envL.dimension(0), envR.dimension(0)};
    size_mps   = spin_dim * envL.dimension(0) * envR.dimension(0);
    t_factorOP = std::make_unique<tid::ur>("Time FactorOp");
    t_genMat   = std::make_unique<tid::ur>("Time genMat");
    t_multOPv  = std::make_unique<tid::ur>("Time MultOpv");
    t_multAx   = std::make_unique<tid::ur>("Time MultAx");
}

template MatVecMPOS<cplx>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvEne &> &envs_);
template MatVecMPOS<real>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvEne &> &envs_);
template MatVecMPOS<cplx>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvVar &> &envs_);
template MatVecMPOS<real>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvVar &> &envs_);

template<typename T>
int MatVecMPOS<T>::rows() const {
    return safe_cast<int>(size_mps);
}

template<typename T>
int MatVecMPOS<T>::cols() const {
    return safe_cast<int>(size_mps);
}

template<typename T>
void MatVecMPOS<T>::FactorOP() {
    throw std::runtime_error("template<typename T> void MatVecMPOS<T>::FactorOP(): Not implemented");
}

template<typename T>
void MatVecMPOS<T>::MultAx(T *mps_in_, T *mps_out_) {
    auto token   = t_multAx->tic_token();
    auto mps_in  = Eigen::TensorMap<Eigen::Tensor<T, 3>>(mps_in_, shape_mps);
    auto mps_out = Eigen::TensorMap<Eigen::Tensor<T, 3>>(mps_out_, shape_mps);
    if(mpos.size() == 1) {
        tools::common::contraction::matrix_vector_product(mps_out, mps_in, mpos.front(), envL, envR);
    } else {
        auto &threads = tenx::threads::get();

        constexpr auto shf6            = std::array<long, 6>{0, 5, 1, 2, 3, 4};
        auto           shp4            = std::array<long, 4>{shape_mps[0], shape_mps[2], envL.dimension(1), envL.dimension(2)};
        auto           mps_tmp1        = Eigen::Tensor<T, 4>(shp4);
        auto           mps_tmp2        = Eigen::Tensor<T, 4>();
        mps_tmp1.device(*threads->dev) = mps_in.contract(envL, tenx::idx({1}, {0}));
        auto mpodimprod                = [&](size_t fr, size_t to) -> long {
            long prod = 1;
            for(size_t idx = fr; idx < to; ++idx) {
                if(idx >= mpos.size()) break;
                prod *= mpos[idx].dimension(3);
            }
            return prod;
        };
        for(size_t idx = 0; idx < mpos.size(); ++idx) {
            const auto &mpo = mpos[idx];
            auto        d0  = mpodimprod(0, idx);
            auto        d1  = mpo.dimension(3);
            auto        d2  = mpodimprod(idx + 1, mpos.size());
            assert(d0 * d1 * d2 == mps_tmp1.dimension(0));
            auto old_shp6 = std::array<long, 6>{d0, d1, d2, mps_tmp1.dimension(1), mps_tmp1.dimension(2), mps_tmp1.dimension(3)};
            auto new_shp4 = std::array<long, 4>{mps_tmp1.dimension(0), mps_tmp1.dimension(1), mps_tmp1.dimension(2), mpo.dimension(1)};
            mps_tmp2.resize(new_shp4);
            mps_tmp2.device(*threads->dev) = mps_tmp1.reshape(old_shp6).contract(mpo, tenx::idx({5, 1}, {0, 2})).shuffle(shf6).reshape(new_shp4);
            mps_tmp1                       = std::move(mps_tmp2);
        }
        mps_out.device(*threads->dev) = mps_tmp1.contract(envR, tenx::idx({1, 3}, {0, 2}));
    }
    num_mv++;
}

template<typename T>
void MatVecMPOS<T>::MultAx(void *x, int *ldx, void *y, int *ldy, int *blockSize, [[maybe_unused]] primme_params *primme, [[maybe_unused]] int *err) {
    // #pragma omp parallel for for schedule(static, 8)
    for(int i = 0; i < *blockSize; i++) {
        T *mps_in_ptr  = static_cast<T *>(x) + *ldx * i;
        T *mps_out_ptr = static_cast<T *>(y) + *ldy * i;
        MultAx(mps_in_ptr, mps_out_ptr);
    }
    *err = 0;
}

template<typename T>
void MatVecMPOS<T>::MultOPv([[maybe_unused]] T *mps_in_, [[maybe_unused]] T *mps_out_) {
    throw std::runtime_error("template<typename T> void MatVecMPOS<T>::MultOPv(T *mps_in_, T *mps_out_): Not implemented");
}

template<typename T>
void MatVecMPOS<T>::MultOPv([[maybe_unused]] void *x, [[maybe_unused]] int *ldx, [[maybe_unused]] void *y, [[maybe_unused]] int *ldy,
                            [[maybe_unused]] int *blockSize, [[maybe_unused]] primme_params *primme, [[maybe_unused]] int *err) {
    throw std::runtime_error("template<typename T> void MatVecMPOS<T>::MultOPv(void *x, int *ldx, void *y, int *ldy, int *blockSize, [[maybe_unused]] "
                             "primme_params *primme, int *err): Not implemented");
}

template<typename T>
void MatVecMPOS<T>::print() const {}

template<typename T>
void MatVecMPOS<T>::reset() {
    if(t_factorOP) t_factorOP->reset();
    if(t_multOPv) t_multOPv->reset();
    if(t_genMat) t_genMat->reset();
    if(t_multAx) t_multAx->reset();
    num_mv = 0;
    num_op = 0;
}

template<typename T>
void MatVecMPOS<T>::set_shift(std::complex<double> shift) {
    // Here we set an energy shift directly on the MPO.
    // This only works if the MPO is not compressed already.
    if(readyShift) return;
    if(sigma == shift) return;
    auto shift_per_mpo = shift / static_cast<double>(mpos.size());
    auto sigma_per_mpo = sigma / static_cast<double>(mpos.size());
    for(size_t idx = 0; idx < mpos.size(); ++idx) {
        // The MPO is a rank4 tensor ijkl where the first 2 ij indices draw a simple
        // rank2 matrix, where each element is also a matrix with the size
        // determined by the last 2 indices kl.
        // When we shift an MPO, all we do is subtract a diagonal matrix from
        // the botton left corner of the ij-matrix.
        auto &mpo  = mpos[idx];
        auto  dims = mpo.dimensions();
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
            mpo.slice(offset4, extent4).reshape(extent2) += id * std::real(sigma_per_mpo - shift_per_mpo);
        else
            mpo.slice(offset4, extent4).reshape(extent2) += id * (sigma_per_mpo - shift_per_mpo);
        eig::log->debug("Shifted MPO {} energy by {:.16f}", idx, shift_per_mpo);
    }
    sigma = shift;

    readyShift = true;
}

template<typename T>
void MatVecMPOS<T>::set_mode(const eig::Form form_) {
    form = form_;
}
template<typename T>
void MatVecMPOS<T>::set_side(const eig::Side side_) {
    side = side_;
}

template<typename T>
T MatVecMPOS<T>::get_shift() const {
    if constexpr(std::is_same_v<T, real>)
        return std::real(sigma);
    else
        return sigma;
}

template<typename T>
eig::Form MatVecMPOS<T>::get_form() const {
    return form;
}
template<typename T>
eig::Side MatVecMPOS<T>::get_side() const {
    return side;
}
template<typename T>
eig::Type MatVecMPOS<T>::get_type() const {
    if constexpr(std::is_same_v<T, real>)
        return eig::Type::REAL;
    else if constexpr(std::is_same_v<T, cplx>)
        return eig::Type::CPLX;
    else
        throw std::runtime_error("Unsupported type");
}

template<typename T>
const std::vector<Eigen::Tensor<T, 4>> &MatVecMPOS<T>::get_mpos() const {
    return mpos;
}
template<typename T>
const Eigen::Tensor<T, 3> &MatVecMPOS<T>::get_envL() const {
    return envL;
}
template<typename T>
const Eigen::Tensor<T, 3> &MatVecMPOS<T>::get_envR() const {
    return envR;
}

template<typename T>
long MatVecMPOS<T>::get_size() const {
    return size_mps;
}

template<typename T>
std::array<long, 3> MatVecMPOS<T>::get_shape_mps() const {
    return shape_mps;
}
template<typename T>
std::vector<std::array<long, 4>> MatVecMPOS<T>::get_shape_mpo() const {
    std::vector<std::array<long, 4>> shapes;
    for(const auto &mpo : mpos) shapes.emplace_back(mpo.dimensions());
    return shapes;
}

template<typename T>
std::array<long, 3> MatVecMPOS<T>::get_shape_envL() const {
    return envL.dimensions();
}
template<typename T>
std::array<long, 3> MatVecMPOS<T>::get_shape_envR() const {
    return envR.dimensions();
}
template<typename T>
Eigen::Tensor<T, 6> MatVecMPOS<T>::get_tensor() const {
    throw std::runtime_error("template<typename T> void MatVecMPOS<T>::get_tensor(): Not implemented");
}

template<typename T>
typename MatVecMPOS<T>::MatrixType MatVecMPOS<T>::get_matrix() const {
    return tenx::MatrixCast(get_tensor(), rows(), cols());
}

template<typename T>
bool MatVecMPOS<T>::isReadyFactorOp() const {
    return readyFactorOp;
}
template<typename T>
bool MatVecMPOS<T>::isReadyShift() const {
    return readyShift;
}

// Explicit instantiations
template class MatVecMPOS<double>;
template class MatVecMPOS<std::complex<double>>;
