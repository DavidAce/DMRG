#define DMRG_ENABLE_TBLIS
#include "matvec_mpos.h"
#include "../log.h"
#include "general/sfinae.h"
#include "math/eig/solver.h"
#include "math/linalg/tensor.h"
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
#include <unsupported/Eigen/CXX11/Tensor>

#if defined(DMRG_ENABLE_TBLIS)
    #include <tblis/tblis.h>
    #include <tblis/util/thread.h>
    #include <tci/tci_config.h>
#endif
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

    shape_mps = {spin_dim, envL.dimension(0), envR.dimension(0)};
    size_mps  = spin_dim * envL.dimension(0) * envR.dimension(0);

    // If we have 5 or fewer mpos, it is faster to just merge them once and apply them in one contraction.
    if(mpos.size() <= 5) {
        constexpr auto contract_idx    = tenx::idx({1}, {0});
        constexpr auto shuffle_idx     = tenx::array6{0, 3, 1, 4, 2, 5};
        auto          &threads         = tenx::threads::get();
        auto           contracted_mpos = mpos.front();
        for(size_t idx = 0; idx + 1 < mpos.size(); ++idx) {
            const auto &mpoL = idx == 0 ? mpos[idx] : contracted_mpos;
            const auto &mpoR = mpos[idx + 1];
            auto new_dims    = std::array{mpoL.dimension(0), mpoR.dimension(1), mpoL.dimension(2) * mpoR.dimension(2), mpoL.dimension(3) * mpoR.dimension(3)};
            auto temp        = Eigen::Tensor<T, 4>(new_dims);
            temp.device(*threads->dev) = mpoL.contract(mpoR, contract_idx).shuffle(shuffle_idx).reshape(new_dims);
            contracted_mpos            = std::move(temp);
        }
        mpos = {contracted_mpos}; // Replace by a single pre-contracted mpo
                                  //         {                         // A quick experiment
                                  // #pragma message "Remove a quick experiment with matvec_mpos"
                                  //             auto                svd           = svd::solver();
                                  //             auto                dim4          = contracted_mpos.dimensions();
                                  //             auto                dim2          = tenx::array2{dim4[2] * dim4[0], dim4[3] * dim4[1]};
                                  //             Eigen::Tensor<T, 4> shuffled_mpos = contracted_mpos.shuffle(tenx::array4{2, 0, 3, 1});
                                  //             auto [U, S, V]                    = svd.schmidt(shuffled_mpos, svd::config(dim2[0], 1e-16));
        //             tools::log->info("MPO decomposition dim4 {} -> {} {} | trnc {:.3e}", dim4, U.dimensions(), V.dimensions(), svd.get_truncation_error());
        //             tools::log->info("S: \n{}\n", linalg::tensor::to_string(S, 16));
        //         }
        //         { // A quick experiment
        // #pragma message "Remove a quick experiment with matvec_mpos"
        //             auto                svd           = svd::solver();
        //             auto                dim4          = contracted_mpos.dimensions();
        //             auto                dim2          = tenx::array2{dim4[2],  dim4[0]* dim4[3] * dim4[1]};
        //             Eigen::Tensor<T, 2> shuffled_mpos = contracted_mpos.shuffle(tenx::array4{2, 3, 0, 1}).reshape(dim2);
        //             auto [U, S, V]                    = svd.decompose(shuffled_mpos, svd::config(dim2[0], 1e-16));
        //             tools::log->info("MPO decomposition dim4 {} -> {} {} | trnc {:.3e}", dim4, U.dimensions(), V.dimensions(), svd.get_truncation_error());
        //             tools::log->info("S: \n{}\n", linalg::tensor::to_string(S, 16));
        //         }

    } else {
        // We pre-shuffle each mpo to speed up the sequential contraction
        for(const auto &mpo : mpos) mpos_shf.emplace_back(mpo.shuffle(tenx::array4{2, 3, 0, 1}));
    }

    t_factorOP = std::make_unique<tid::ur>("Time FactorOp");
    t_genMat   = std::make_unique<tid::ur>("Time genMat");
    t_multOPv  = std::make_unique<tid::ur>("Time MultOpv");
    t_multAx   = std::make_unique<tid::ur>("Time MultAx");
    t_multPc   = std::make_unique<tid::ur>("Time MultPc");
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

std::vector<long> indices(long x, long rank, long size) {
    std::vector<long> indices(rank);
    for(int i = 0; i < rank; i++) {
        indices[i] = x % size;
        x /= size;
    }
    return indices;
}

template<size_t rank>
constexpr std::array<long, rank> indices(long flatindex, const std::array<long, rank> &dimensions) {
    std::array<long, rank> indices;
    for(size_t i = 0; i < rank; i++) {
        indices[i] = flatindex % dimensions[i];
        flatindex /= dimensions[i];
    }
    return indices;
}

std::vector<long> indices(long flatindex, size_t rank, const std::vector<long> &dimensions) {
    std::vector<long> indices(rank);
    for(size_t i = 0; i < rank; i++) {
        indices[i] = flatindex % dimensions[i];
        flatindex /= dimensions[i];
    }
    return indices;
}

template<typename T>
void MatVecMPOS<T>::thomas(long size, T *x, T *const dl, T *const dm, T *const du) {
    /*
     solves Ax = d, where A is a tridiagonal matrix consisting of vectors a, b, c
     X = number of equations
     x[] = initially contains the input v, and returns x. indexed from [0, ..., X - 1]
     dl[] = lower diagonal, indexed from [1, ..., X - 1]
     dm[] = main diagonal, indexed from [0, ..., X - 1]
     du[] = upper diagonal, indexed from [0, ..., X - 2]
     not performed in this example: manual expensive common subexpression elimination
     */
    diagtemp.resize(size);
    auto &t = diagtemp;

    diagtemp[0] = du[0] / dm[0];
    x[0]        = x[0] / dm[0];

    /* loop from 1 to X - 1 inclusive */
    for(long ix = 1; ix < size - 1; ix++) {
        // if(ix < size - 1) {  }
        auto denom = (dm[ix] - dl[ix] * t[ix - 1]);
        t[ix]      = du[ix] / denom;
        x[ix]      = (x[ix] - dl[ix] * x[ix - 1]) / denom;
    }

    /* loop from X - 2 to 0 inclusive */
    for(long ix = size - 2; ix >= 0; ix--) x[ix] -= t[ix] * x[ix + 1];
}

template<typename T>
void MatVecMPOS<T>::thomas2(long size, T *x, T *const a, T *const b, T *const c) {
    /*
    // size is the number of unknowns

    |b0 c0 0 ||x0| |d0|
    |a1 b1 c1||x1|=|d1|
    |0  a2 b2||x2| |d2|

    */
    size--;
    c[0] /= b[0];
    x[0] /= b[0];
    for(long i = 1; i < size; i++) {
        c[i] /= b[i] - a[i] * c[i - 1];
        x[i] = (x[i] - a[i] * x[i - 1]) / (b[i] - a[i] * c[i - 1]);
    }

    x[size] = (x[size] - a[size] * x[size - 1]) / (b[size] - a[size] * c[size - 1]);

    for(long i = size; i-- > 0;) { x[i] -= c[i] * x[i + 1]; }
}

template<typename T>
typename MatVecMPOS<T>::VectorType MatVecMPOS<T>::get_diagonal(long offset) {
    assert(std::abs(offset) <= 1);
    auto  diagonal = VectorType(size_mps);
    auto &threads  = tenx::threads::get();
    auto  ext_j    = std::array<long, 3>{1, 1, envL.dimension(2)};
    auto  ext_k    = std::array<long, 3>{1, 1, envR.dimension(2)};
    auto  spindims = std::vector<long>();
    spindims.reserve(mpos.size());
    for(const auto &mpo : mpos) spindims.emplace_back(mpo.dimension(2));
    auto rowindices = std::array<long, 3>{};
    auto colindices = std::array<long, 3>{};
    for(long diagidx = 0; diagidx < size_mps; ++diagidx) {
        if(diagidx + offset < 0 or diagidx + offset >= size_mps) { diagonal(diagidx) = 0.0; }
        if(offset == 0) {
            rowindices = indices(diagidx, shape_mps);
            colindices = rowindices;
        } else {
            rowindices = indices(diagidx, shape_mps);
            colindices = indices(diagidx + offset, shape_mps);
        }
        auto ir = rowindices[0];
        auto jr = rowindices[1];
        auto kr = rowindices[2];
        auto ic = colindices[0];
        auto jc = colindices[1];
        auto kc = colindices[2];

        // Index i is special, since it is composed of all the mpo physical indices.
        // auto idxs  = indices(i, mpos.size(), mpos.front().dimension(2)); // Maps i to tensor indices to select the physical indices in the
        // mpos
        auto irxs  = indices(ir, mpos.size(), spindims); // Maps ir to tensor indices to select the upper phsical indices in the mpos
        auto icxs  = indices(ic, mpos.size(), spindims); // Maps ic to tensor indices to select the lower phsical indices in the mpos
        auto mpo_i = Eigen::Tensor<Scalar, 4>();
        for(size_t mdx = 0; mdx < mpos.size(); ++mdx) {
            if(mdx == 0) {
                auto &mpoL = mpos[mdx];
                auto  offL = std::array<long, 4>{0, 0, irxs[mdx], icxs[mdx]};
                auto  extL = std::array<long, 4>{mpoL.dimension(0), mpoL.dimension(1), 1, 1};
                mpo_i      = mpoL.slice(offL, extL);
            }
            if(mdx + 1 == mpos.size()) break;
            auto dimR = mpos[mdx + 1].dimensions();
            auto offR = std::array<long, 4>{0, 0, irxs[mdx + 1], icxs[mdx + 1]};
            auto extR = std::array<long, 4>{dimR[0], dimR[1], 1, 1};
            auto shp4 = std::array<long, 4>{mpo_i.dimension(0), dimR[1], 1, 1};
            auto temp = Eigen::Tensor<Scalar, 4>(shp4);
            temp.device(*threads->dev) =
                mpo_i.contract(mpos[mdx + 1].slice(offR, extR), tenx::idx({1}, {0})).shuffle(tenx::array6{0, 3, 1, 4, 2, 5}).reshape(shp4);
            mpo_i = temp;
        }

        auto off_j  = std::array<long, 3>{jr, jc, 0};
        auto off_k  = std::array<long, 3>{kr, kc, 0};
        auto envL_j = envL.slice(off_j, ext_j);
        auto envR_k = envR.slice(off_k, ext_k);

        Eigen::Tensor<Scalar, 6> elem(1, 1, 1, 1, 1, 1);
        elem.device(*threads->dev) = envL_j.contract(mpo_i, tenx::idx({2}, {0})).contract(envR_k, tenx::idx({2}, {2})).shuffle(tenx::array6{2, 0, 4, 3, 1, 5});
        diagonal(diagidx)          = elem.coeff(0);
    }
    return diagonal;
}

template<typename T>
void MatVecMPOS<T>::FactorOP() {
    throw except::runtime_error("MatVecMPOS<T>::FactorOP(): not implemented");
}

template<typename T>
void MatVecMPOS<T>::MultOPv([[maybe_unused]] T *mps_in_, [[maybe_unused]] T *mps_out_,[[maybe_unused]] T eval) {
    throw except::runtime_error("void MatVecMPOS<T>::MultOPv(...) not implemented");
}

template<typename T>
void MatVecMPOS<T>::MultOPv([[maybe_unused]] void *x, [[maybe_unused]] int *ldx, [[maybe_unused]] void *y, [[maybe_unused]] int *ldy,
                            [[maybe_unused]] int *blockSize, [[maybe_unused]] primme_params *primme, [[maybe_unused]] int *err) {
    throw except::runtime_error("MatVecMPOS<T>::MultOPv(...): not implemented");
}

template<typename T>
void MatVecMPOS<T>::MultAx(T *mps_in_, T *mps_out_) {
    auto token   = t_multAx->tic_token();
    auto mps_in  = Eigen::TensorMap<Eigen::Tensor<T, 3>>(mps_in_, shape_mps);
    auto mps_out = Eigen::TensorMap<Eigen::Tensor<T, 3>>(mps_out_, shape_mps);
    if(mpos.size() == 1) {
        tools::common::contraction::matrix_vector_product(mps_out, mps_in, mpos.front(), envL, envR);
    } else {
        tools::common::contraction::matrix_vector_product(mps_out, mps_in, mpos_shf, envL, envR);
    }
    num_mv++;
}

template<typename T>
void MatVecMPOS<T>::MultAx(void *x, int *ldx, void *y, int *ldy, int *blockSize, [[maybe_unused]] primme_params *primme, [[maybe_unused]] int *err) {
    for(int i = 0; i < *blockSize; i++) {
        T *mps_in_ptr  = static_cast<T *>(x) + *ldx * i;
        T *mps_out_ptr = static_cast<T *>(y) + *ldy * i;
        MultAx(mps_in_ptr, mps_out_ptr);
    }
    *err = 0;
}

template<typename T>
void MatVecMPOS<T>::CalcPc() {
    // auto token = t_multPc->tic_token();
    if(preconditioner == eig::Preconditioner::DIAG) {
        if(diagonal.size() == 0) {
            eig::log->debug("MatVecMPOS<T>::CalcPc(): calculating the diagonal preconditioner ... ");
            diagonal = get_diagonal(0);
            eig::log->trace("MatVecMPOS<T>::CalcPc(): calculating the diagonal preconditioner ... done");
        }
    } else if(preconditioner == eig::Preconditioner::TRIDIAG) {
        if(diagonal.size() == 0 or diaglower.size() == 0 or diagupper.size() == 0) {
            eig::log->debug("MatVecMPOS<T>::CalcPc(): calculating the tridiagonal preconditioner ... ");
            if(diagonal.size() == 0) diagonal = get_diagonal(0);
            if(diaglower.size() == 0) diaglower = get_diagonal(-1);
            if(diagupper.size() == 0) diagupper = get_diagonal(1);
            eig::log->trace("MatVecMPOS<T>::CalcPc(): calculating the tridiagonal preconditioner ... done");
        }
    }
}

template<typename T>
void MatVecMPOS<T>::MultPc([[maybe_unused]] void *x, [[maybe_unused]] int *ldx, [[maybe_unused]] void *y, [[maybe_unused]] int *ldy,
                           [[maybe_unused]] int *blockSize, [[maybe_unused]] primme_params *primme, [[maybe_unused]] int *err) {
    for(int i = 0; i < *blockSize; i++) {
        T *mps_in_ptr  = static_cast<T *>(x) + *ldx * i;
        T *mps_out_ptr = static_cast<T *>(y) + *ldy * i;
        MultPc(mps_in_ptr, mps_out_ptr, primme->ShiftsForPreconditioner[i]);
    }
    *err = 0;
}

template<typename T>
void MatVecMPOS<T>::MultPc([[maybe_unused]] T *mps_in_, [[maybe_unused]] T *mps_out_, T shift) {
    if(preconditioner == eig::Preconditioner::NONE) return;
    if(diagonal.size() == 0) { CalcPc(); }
    auto token   = t_multPc->tic_token();
    auto mps_in  = Eigen::Map<VectorType>(mps_in_, size_mps);
    auto mps_out = Eigen::Map<VectorType>(mps_out_, size_mps);
    using ArrayType = Eigen::Array<T,Eigen::Dynamic,1>;
    ArrayType diagshf = (diagonal.array() - shift).cwiseAbs();
    // Use a small cutoff for small differences following Eq 2.11 in https://sdm.lbl.gov/~kewu/ps/thesis.pdf
    static constexpr auto eps = std::numeric_limits<double>::epsilon();
    for(auto &d : diagshf) {
        if(d <= eps) d = eps;
        // if(d >= 0 and d <= +eps) d = +eps;
        // if(d <= 0 and d >= -eps) d = -eps;
    }
    if(preconditioner == eig::Preconditioner::DIAG) {
        mps_out = diagshf.array().cwiseInverse().cwiseProduct(mps_in.array());
    } else if(preconditioner == eig::Preconditioner::TRIDIAG) {
        mps_out = mps_in; // copy the values
        thomas(size_mps, mps_out_, diaglower.data(), diagshf.data(), diagupper.data());
    }
    num_pc++;
}

template<typename T>
void MatVecMPOS<T>::print() const {}

template<typename T>
void MatVecMPOS<T>::reset() {
    if(t_factorOP) t_factorOP->reset();
    if(t_multOPv) t_multOPv->reset();
    if(t_genMat) t_genMat->reset();
    if(t_multAx) t_multAx->reset();
    if(t_multPc) t_multPc->reset();
    num_mv = 0;
    num_op = 0;
    num_pc = 0;
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
    if(not mpos_shf.empty()) {
        mpos_shf.clear();
        for(const auto &mpo : mpos) mpos_shf.emplace_back(mpo.shuffle(tenx::array4{2, 3, 0, 1}));
    }

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
