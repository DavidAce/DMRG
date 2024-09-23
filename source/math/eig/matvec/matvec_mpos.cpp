#define DMRG_ENABLE_TBLIS
#include "matvec_mpos.h"
#include "../log.h"
#include "config/settings.h"
#include "general/sfinae.h"
#include "math/eig/solver.h"
#include "math/linalg/matrix.h"
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
#include <Eigen/Eigenvalues>
#include <h5pp/h5pp.h>
#include <primme/primme.h>
#include <queue>
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
    mpos_A.reserve(mpos_.size());
    fullsystem = envs_.L.get_sites() == 0 and envs_.R.get_sites() == 0; //  mpos.size() == settings::model::model_size;

    if constexpr(std::is_same_v<EnvType, EnvEne>) {
        if constexpr(std::is_same_v<T, cplx>) {
            for(const auto &mpo_ : mpos_) mpos_A.emplace_back(mpo_.get().MPO());
            envL_A = envs_.L.get_block();
            envR_A = envs_.R.get_block();
        } else {
            if constexpr(eig::debug) {
                if(not tenx::isReal(envs_.L.get_block())) throw except::runtime_error("envL is not real");
                if(not tenx::isReal(envs_.R.get_block())) throw except::runtime_error("envR is not real");
            }
            for(const auto &mpo_ : mpos_) {
                if(not tenx::isReal(mpo_.get().MPO())) throw except::runtime_error("mpo is not real");
                mpos_A.emplace_back(mpo_.get().MPO().real());
            }
            envL_A = envs_.L.get_block().real();
            envR_A = envs_.R.get_block().real();
        }
    }
    if constexpr(std::is_same_v<EnvType, EnvVar>) {
        if constexpr(std::is_same_v<T, cplx>) {
            for(const auto &mpo_ : mpos_) mpos_A.emplace_back(mpo_.get().MPO2());
            envL_A = envs_.L.get_block();
            envR_A = envs_.R.get_block();
        } else {
            if constexpr(eig::debug) {
                if(not tenx::isReal(envs_.L.get_block())) throw except::runtime_error("envL is not real");
                if(not tenx::isReal(envs_.R.get_block())) throw except::runtime_error("envR is not real");
            }
            for(const auto &mpo_ : mpos_) {
                if(not tenx::isReal(mpo_.get().MPO2())) throw except::runtime_error("mpo is not real");
                mpos_A.emplace_back(mpo_.get().MPO2().real());
            }
            envL_A = envs_.L.get_block().real();
            envR_A = envs_.R.get_block().real();
        }
    }

    long spin_dim = 1;
    for(const auto &mpo : mpos_A) spin_dim *= mpo.dimension(2);
    spindims.reserve(mpos_A.size());
    for(const auto &mpo : mpos_A) spindims.emplace_back(mpo.dimension(2));

    shape_mps = {spin_dim, envL_A.dimension(0), envR_A.dimension(0)};
    size_mps  = spin_dim * envL_A.dimension(0) * envR_A.dimension(0);

    // if(mpos.size() == settings::model::model_size) {
    //     auto t_spm = tid::ur("t_spm");
    //     t_spm.tic();
    //     eig::log->info("making sparse matrix ... ", t_spm.get_last_interval());
    //     sparseMatrix = get_sparse_matrix();
    //     t_spm.toc();
    //     eig::log->info("making sparse matrix ... {:.3e} s | nnz {} / {} = {:.16f}", t_spm.get_last_interval(), sparseMatrix.nonZeros(), sparseMatrix.size(),
    //                    static_cast<double>(sparseMatrix.nonZeros()) / static_cast<double>(sparseMatrix.size()));
    // }

    // If we have 5 or fewer mpos, it is faster to just merge them once and apply them in one contraction.
    if(mpos_A.size() <= 5) {
        constexpr auto contract_idx    = tenx::idx({1}, {0});
        constexpr auto shuffle_idx     = tenx::array6{0, 3, 1, 4, 2, 5};
        auto          &threads         = tenx::threads::get();
        auto           contracted_mpos = mpos_A.front();
        for(size_t idx = 0; idx + 1 < mpos_A.size(); ++idx) {
            const auto &mpoL = idx == 0 ? mpos_A[idx] : contracted_mpos;
            const auto &mpoR = mpos_A[idx + 1];
            auto new_dims    = std::array{mpoL.dimension(0), mpoR.dimension(1), mpoL.dimension(2) * mpoR.dimension(2), mpoL.dimension(3) * mpoR.dimension(3)};
            auto temp        = Eigen::Tensor<T, 4>(new_dims);
            temp.device(*threads->dev) = mpoL.contract(mpoR, contract_idx).shuffle(shuffle_idx).reshape(new_dims);
            contracted_mpos            = std::move(temp);
        }
        mpos_A   = {contracted_mpos}; // Replace by a single pre-contracted mpo
        spindims = {mpos_A.front().dimension(2)};
    } else {
        // We pre-shuffle each mpo to speed up the sequential contraction
        for(const auto &mpo : mpos_A) mpos_A_shf.emplace_back(mpo.shuffle(tenx::array4{2, 3, 0, 1}));
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
template<typename EnvTypeA, typename EnvTypeB>
MatVecMPOS<T>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, /*!< The Hamiltonian MPO's  */
                          const env_pair<const EnvTypeA &>                         &enva_, /*!< The left and right environments.  */
                          const env_pair<const EnvTypeB &>                         &envb_)
    : MatVecMPOS(mpos_, enva_) {
    // static_assert(sfinae::is_any_v<EnvTypeA, EnvVar>);
    // static_assert(sfinae::is_any_v<EnvTypeB, EnvEne>);
    if constexpr(std::is_same_v<EnvTypeB, EnvEne>) {
        if constexpr(std::is_same_v<T, cplx>) {
            for(const auto &mpo_ : mpos_) mpos_B.emplace_back(mpo_.get().MPO());
            envL_B = envb_.L.get_block();
            envR_B = envb_.R.get_block();
        } else {
            for(const auto &mpo_ : mpos_) mpos_B.emplace_back(mpo_.get().MPO().real());
            envL_B = envb_.L.get_block().real();
            envR_B = envb_.R.get_block().real();
        }
    }
    if constexpr(std::is_same_v<EnvTypeB, EnvVar>) {
        if constexpr(std::is_same_v<T, cplx>) {
            for(const auto &mpo_ : mpos_) mpos_B.emplace_back(mpo_.get().MPO2());
            envL_B = envb_.L.get_block();
            envR_B = envb_.R.get_block();
        } else {
            for(const auto &mpo_ : mpos_) mpos_B.emplace_back(mpo_.get().MPO2().real());
            envL_B = envb_.L.get_block().real();
            envR_B = envb_.R.get_block().real();
        }
    }

    if(mpos_B.size() <= 5) {
        constexpr auto contract_idx    = tenx::idx({1}, {0});
        constexpr auto shuffle_idx     = tenx::array6{0, 3, 1, 4, 2, 5};
        auto          &threads         = tenx::threads::get();
        auto           contracted_mpos = mpos_B.front();
        for(size_t idx = 0; idx + 1 < mpos_B.size(); ++idx) {
            const auto &mpoL = idx == 0 ? mpos_B[idx] : contracted_mpos;
            const auto &mpoR = mpos_B[idx + 1];
            auto new_dims    = std::array{mpoL.dimension(0), mpoR.dimension(1), mpoL.dimension(2) * mpoR.dimension(2), mpoL.dimension(3) * mpoR.dimension(3)};
            auto temp        = Eigen::Tensor<T, 4>(new_dims);
            temp.device(*threads->dev) = mpoL.contract(mpoR, contract_idx).shuffle(shuffle_idx).reshape(new_dims);
            contracted_mpos            = std::move(temp);
        }
        mpos_B = {contracted_mpos}; // Replace by a single pre-contracted mpo
    } else {
        // We pre-shuffle each mpo to speed up the sequential contraction
        for(const auto &mpo : mpos_B) mpos_B_shf.emplace_back(mpo.shuffle(tenx::array4{2, 3, 0, 1}));
    }
}

template MatVecMPOS<real>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvVar &> &enva_,
                                      const env_pair<const EnvEne &> &envb);
template MatVecMPOS<cplx>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvVar &> &enva_,
                                      const env_pair<const EnvEne &> &envb);

template MatVecMPOS<real>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvEne &> &enva_,
                                      const env_pair<const EnvVar &> &envb);
template MatVecMPOS<cplx>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvEne &> &enva_,
                                      const env_pair<const EnvVar &> &envb);

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
//
// template<typename T>
// void MatVecMPOS<T>::thomas(long size, T *x, T *const dl, T *const dm, T *const du) {
//     /*
//      solves Ax = d, where A is a tridiagonal matrix consisting of vectors a, b, c
//      X = number of equations
//      x[] = initially contains the input v, and returns x. indexed from [0, ..., X - 1]
//      dl[] = lower diagonal, indexed from [1, ..., X - 1]
//      dm[] = main diagonal, indexed from [0, ..., X - 1]
//      du[] = upper diagonal, indexed from [0, ..., X - 2]
//      not performed in this example: manual expensive common subexpression elimination
//      */
//     diagtemp.resize(size);
//     auto &t = diagtemp;
//
//     diagtemp[0] = du[0] / dm[0];
//     x[0]        = x[0] / dm[0];
//
//     /* loop from 1 to X - 1 inclusive */
//     for(long ix = 1; ix < size - 1; ix++) {
//         // if(ix < size - 1) {  }
//         auto denom = (dm[ix] - dl[ix] * t[ix - 1]);
//         t[ix]      = du[ix] / denom;
//         x[ix]      = (x[ix] - dl[ix] * x[ix - 1]) / denom;
//     }
//
//     /* loop from X - 2 to 0 inclusive */
//     for(long ix = size - 2; ix >= 0; ix--) x[ix] -= t[ix] * x[ix + 1];
// }
//
// template<typename T>
// void MatVecMPOS<T>::thomas2(long size, T *x, T *const a, T *const b, T *const c) {
//     /*
//     // size is the number of unknowns
//
//     |b0 c0 0 ||x0| |d0|
//     |a1 b1 c1||x1|=|d1|
//     |0  a2 b2||x2| |d2|
//
//     */
//     size--;
//     c[0] /= b[0];
//     x[0] /= b[0];
//     for(long i = 1; i < size; i++) {
//         c[i] /= b[i] - a[i] * c[i - 1];
//         x[i] = (x[i] - a[i] * x[i - 1]) / (b[i] - a[i] * c[i - 1]);
//     }
//
//     x[size] = (x[size] - a[size] * x[size - 1]) / (b[size] - a[size] * c[size - 1]);
//
//     for(long i = size; i-- > 0;) { x[i] -= c[i] * x[i + 1]; }
// }
template<typename TI>
TI HammingDist(TI x, TI y) {
    TI dist = 0;
    TI val  = x ^ y; // calculate differ bit
    while(val)       // this dist variable calculates set bits in a loop
    {
        ++dist;
        if(dist > 4) return dist;
        val &= val - 1;
    }
    return dist;
}

template<typename T>
T MatVecMPOS<T>::get_matrix_element(long I, long J, const std::vector<Eigen::Tensor<T, 4>> &MPOS, const Eigen::Tensor<T, 3> &ENVL,
                                    const Eigen::Tensor<T, 3> &ENVR) const {
    if(I < 0 or I >= size_mps) return 0;
    if(J < 0 or J >= size_mps) return 0;

    if(MPOS.empty()) {
        return I == J ? T(1.0) : T(0.0); // Assume an identity matrix
    }

    auto rowindices = std::array<long, 3>{};
    auto colindices = std::array<long, 3>{};

    if(I == J) {
        rowindices = indices(I, shape_mps);
        colindices = rowindices;
    } else {
        rowindices = indices(I, shape_mps);
        colindices = indices(J, shape_mps);
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

    bool shouldBeZero = false;
    // if(fullsystem) {
    //     auto pr = std::popcount(static_cast<unsigned int>(ir));
    //     auto pc = std::popcount(static_cast<unsigned int>(ic));
    //     if(((pr ^ pc) & 1) != 0) { return 0.0; } // different parity
    //     auto hd = std::popcount(static_cast<unsigned int>(ir) ^ static_cast<unsigned int>(ic));
    //     if(hd > 4) { return 0.0; }
    //     auto pd = std::abs(pr - pc);
    //     if(pd > 4) { return 0.0; } // popcount difference is larger than 4
    //     if(pd & 1) { return 0.0; } // popcount difference is odd
    //     // if(((pr ^ pc) & 1) != 0) { shouldBeZero = true; } // different parity
    //     // if(hd > 4) { shouldBeZero = true; }
    //     // if(pd > 4) { shouldBeZero = true; } // popcount difference is larger than 4
    //     // if(pd & 1) { shouldBeZero = true; } // popcount difference is odd
    // }

    auto irxs  = indices(ir, MPOS.size(), spindims); // Maps ir to tensor indices to select the upper physical indices in the mpos
    auto icxs  = indices(ic, MPOS.size(), spindims); // Maps ic to tensor indices to select the lower physical indices in the mpos
    auto mpo_i = Eigen::Tensor<Scalar, 4>();
    auto temp  = Eigen::Tensor<Scalar, 4>();

    constexpr auto shf = tenx::array6{0, 3, 1, 4, 2, 5};
    for(size_t mdx = 0; mdx < MPOS.size(); ++mdx) {
        const auto &mpo = MPOS[mdx];
        auto        dim = mpo.dimensions();
        auto        off = std::array<long, 4>{0, 0, irxs[mdx], icxs[mdx]};
        auto        ext = std::array<long, 4>{dim[0], dim[1], 1, 1};
        if(mdx == 0) {
            mpo_i = mpo.slice(off, ext);
            if(tenx::isZero(mpo_i, std::numeric_limits<double>::epsilon())) { return 0.0; }
            continue;
        }
        auto shp = std::array<long, 4>{mpo_i.dimension(0), dim[1], 1, 1};
        temp.resize(shp);
        temp  = mpo_i.contract(mpo.slice(off, ext), tenx::idx({1}, {0})).shuffle(shf).reshape(shp);
        mpo_i = std::move(temp);
        if(tenx::isZero(mpo_i, std::numeric_limits<double>::epsilon())) {
            // eig::log->info("({}, {}) = < {} | {} > = 0 (mdx {}, pr {}, pc {}, pd {}, hd {})", I, J, irxs, icxs, mdx, pr, pc, pd,
            // HammingDist(static_cast<size_t>(ir), static_cast<size_t>(ic)));
            return 0.0;
        }
    }

    if(fullsystem and mpo_i.size() == 1) {
        if(mpo_i.coeff(0) != 0.0 and shouldBeZero) { eig::log->info("({}, {}) = < {} | {} > = {:.16f}", I, J, irxs, icxs, mpo_i.coeff(0)); }
        // eig::log->info("({}, {}) = < {} | {} > = {:.16f} (pr {}, pc {} diff {} | hd {})", I, J, irxs, icxs, mpo_i.coeff(0), pr, pc, pd, hd);
        return mpo_i.coeff(0);
    }

    auto ext_j = std::array<long, 3>{1, 1, ENVL.dimension(2)};
    auto ext_k = std::array<long, 3>{1, 1, ENVR.dimension(2)};
    auto off_j = std::array<long, 3>{jr, jc, 0};
    auto off_k = std::array<long, 3>{kr, kc, 0};
    // auto envL_j = envL.slice(off_j, ext_j);
    // auto envR_k = envR.slice(off_k, ext_k);
    auto envL_j     = Eigen::Tensor<Scalar, 3>(ENVL.slice(off_j, ext_j));
    auto envR_k     = Eigen::Tensor<Scalar, 3>(ENVR.slice(off_k, ext_k));
    auto envL_j_map = Eigen::Map<const VectorType>(envL_j.data(), envL_j.size());
    auto envR_k_map = Eigen::Map<const VectorType>(envR_k.data(), envR_k.size());
    auto mpo_i_map  = Eigen::Map<const MatrixType>(mpo_i.data(), mpo_i.dimension(0), mpo_i.dimension(1));
    return envL_j_map.transpose() * mpo_i_map * envR_k_map;
    // Eigen::Tensor<Scalar, 6> elem(1, 1, 1, 1, 1, 1);
    // elem = envL_j.contract(mpo_i, tenx::idx({2}, {0})).contract(envR_k, tenx::idx({2}, {2})).shuffle(tenx::array6{2, 0, 4, 3, 1, 5});
    // return elem.coeff(0);
}

template<typename T>
typename MatVecMPOS<T>::VectorType MatVecMPOS<T>::get_diagonal_new(long offset, const std::vector<Eigen::Tensor<T, 4>> &MPOS, const Eigen::Tensor<T, 3> &ENVL,
                                                                   const Eigen::Tensor<T, 3> &ENVR) const {
    if(MPOS.empty()) return VectorType::Ones(size_mps); // Assume an identity matrix
    auto res = VectorType(size_mps);
#pragma omp parallel for
    for(long I = 0; I < size_mps; ++I) { res[I] = get_matrix_element(I, I + offset, MPOS, ENVL, ENVR); }
    return res;
}

template<typename T>
typename MatVecMPOS<T>::MatrixType MatVecMPOS<T>::get_diagonal_block(long offset, long extent, const std::vector<Eigen::Tensor<T, 4>> &MPOS,
                                                                     const Eigen::Tensor<T, 3> &ENVL, const Eigen::Tensor<T, 3> &ENVR) const {
    if(MPOS.empty()) return MatrixType::Ones(extent, extent);
    auto res = MatrixType(extent, extent);
    if(offset >= size_mps) { return res; }
    extent = std::min(extent, size_mps);
#pragma omp parallel for collapse(2)
    for(long J = 0; J < extent; J++) {
        for(long I = J; I < extent; I++) {
            if(I + offset >= size_mps) continue;
            if(J + offset >= size_mps) continue;
            // res.template selfadjointView<Eigen::Lower>()(I, J) = get_matrix_element(I + offset, J + offset); // Lower part is sufficient
            auto elem = get_matrix_element(I + offset, J + offset, MPOS, ENVL, ENVR);
            res(I, J) = elem;
            if constexpr(std::is_same_v<T, cplx>)
                res(J, I) = std::conj(elem);
            else
                res(J, I) = elem;
        }
    }
    return res;
}

template<typename T>
typename MatVecMPOS<T>::VectorType MatVecMPOS<T>::get_row(long row_idx, const std::vector<Eigen::Tensor<T, 4>> &MPOS, const Eigen::Tensor<T, 3> &ENVL,
                                                          const Eigen::Tensor<T, 3> &ENVR) const {
    auto res = VectorType(size_mps);
#pragma omp parallel for
    for(long J = 0; J < size_mps; ++J) { res[J] = get_matrix_element(row_idx, J, MPOS, ENVL, ENVR); }
    return res;
}
template<typename T>
typename MatVecMPOS<T>::VectorType MatVecMPOS<T>::get_col(long col_idx, const std::vector<Eigen::Tensor<T, 4>> &MPOS, const Eigen::Tensor<T, 3> &ENVL,
                                                          const Eigen::Tensor<T, 3> &ENVR) const {
    auto res = VectorType(size_mps);
#pragma omp parallel for
    for(long I = 0; I < size_mps; ++I) { res[I] = get_matrix_element(I, col_idx, MPOS, ENVL, ENVR); }
    return res;
}

// template<typename T>
// typename MatVecMPOS<T>::VectorType MatVecMPOS<T>::get_diagonal(long offset, const std::vector<Eigen::Tensor<T, 4>> &MPOS, const Eigen::Tensor<T, 3> &ENVL,
// const Eigen::Tensor<T, 3> &ENVR) const {
//     // auto &threads = tenx::threads::get();
//     auto res   = VectorType(size_mps);
//     auto ext_j = std::array<long, 3>{1, 1, envL_A.dimension(2)};
//     auto ext_k = std::array<long, 3>{1, 1, envR_A.dimension(2)};
//
//     auto rowindices = std::array<long, 3>{};
//     auto colindices = std::array<long, 3>{};
//     for(long diagidx = 0; diagidx < size_mps; ++diagidx) {
//         if(diagidx + offset < 0 or diagidx + offset >= size_mps) {
//             res[diagidx] = 0.0;
//             continue;
//         }
//         if(offset == 0) {
//             rowindices = indices(diagidx, shape_mps);
//             colindices = rowindices;
//         } else {
//             rowindices = indices(diagidx, shape_mps);
//             colindices = indices(diagidx + offset, shape_mps);
//         }
//         auto ir = rowindices[0];
//         auto jr = rowindices[1];
//         auto kr = rowindices[2];
//         auto ic = colindices[0];
//         auto jc = colindices[1];
//         auto kc = colindices[2];
//
//         // Index i is special, since it is composed of all the mpo physical indices.
//         // auto idxs  = indices(i, mpos.size(), mpos.front().dimension(2)); // Maps i to tensor indices to select the physical indices in the
//         // mpos
//         auto irxs   = indices(ir, mpos_A.size(), spindims); // Maps ir to tensor indices to select the upper phsical indices in the mpos
//         auto icxs   = indices(ic, mpos_A.size(), spindims); // Maps ic to tensor indices to select the lower phsical indices in the mpos
//         auto mpo_i  = Eigen::Tensor<Scalar, 4>();
//         auto temp   = Eigen::Tensor<Scalar, 4>();
//         bool isZero = false;
//         for(size_t mdx = 0; mdx < mpos_A.size(); ++mdx) {
//             if(mdx == 0) {
//                 auto &mpoL = mpos_A[mdx];
//                 auto  offL = std::array<long, 4>{0, 0, irxs[mdx], icxs[mdx]};
//                 auto  extL = std::array<long, 4>{mpoL.dimension(0), mpoL.dimension(1), 1, 1};
//                 mpo_i      = mpoL.slice(offL, extL);
//             }
//             if(mdx + 1 == mpos_A.size()) break;
//             auto dimR  = mpos_A[mdx + 1].dimensions();
//             auto offR  = std::array<long, 4>{0, 0, irxs[mdx + 1], icxs[mdx + 1]};
//             auto extR  = std::array<long, 4>{dimR[0], dimR[1], 1, 1};
//             auto shp4  = std::array<long, 4>{mpo_i.dimension(0), dimR[1], 1, 1};
//             auto mpo_r = Eigen::Tensor<Scalar, 4>(mpos_A[mdx + 1].slice(offR, extR));
//             temp.resize(shp4);
//             auto mpo_i_map     = Eigen::Map<const MatrixType>(mpo_i.data(), mpo_i.dimension(0), mpo_i.dimension(1));
//             auto mpo_R_map     = Eigen::Map<const MatrixType>(mpo_r.data(), extR[0], extR[1]);
//             auto temp_map      = Eigen::Map<MatrixType>(temp.data(), shp4[0], shp4[1]);
//             temp_map.noalias() = mpo_i_map * mpo_R_map;
//
//             mpo_i = std::move(temp);
//             if(tenx::isZero(mpo_i, std::numeric_limits<double>::epsilon())) {
//                 isZero = true;
//                 break;
//             }
//             // temp  = mpo_i.contract(mpo_r, tenx::idx({1}, {0})).shuffle(tenx::array6{0, 3, 1, 4, 2, 5}).reshape(shp4);
//         }
//         if(isZero) {
//             res[diagidx] = 0.0;
//             continue;
//         }
//
//         auto off_j = std::array<long, 3>{jr, jc, 0};
//         auto off_k = std::array<long, 3>{kr, kc, 0};
//
//         auto envL_j     = Eigen::Tensor<Scalar, 3>(envL_A.slice(off_j, ext_j));
//         auto envR_k     = Eigen::Tensor<Scalar, 3>(envR_A.slice(off_k, ext_k));
//         auto envL_j_map = Eigen::Map<const VectorType>(envL_j.data(), envL_j.size());
//         auto envR_k_map = Eigen::Map<const VectorType>(envR_k.data(), envR_k.size());
//         auto mpo_i_map  = Eigen::Map<const MatrixType>(mpo_i.data(), mpo_i.dimension(0), mpo_i.dimension(1));
//         res[diagidx]    = envL_j_map.transpose() * mpo_i_map * envR_k_map;
//
//         // Eigen::Tensor<Scalar, 6> elem(1, 1, 1, 1, 1, 1);
//         // elem         = envL_j.contract(mpo_i, tenx::idx({2}, {0})).contract(envR_k, tenx::idx({2}, {2})).shuffle(tenx::array6{2, 0, 4, 3, 1, 5});
//         // res[diagidx] = elem.coeff(0);
//     }
//     // auto   res2    = get_diagonal_new(offset);
//     // long   idxd    = -1;
//     // double maxdiff = (res - res2).cwiseAbs().maxCoeff(&idxd);
//     // if(maxdiff > 1e-12) {
//     //     eig::log->error("diagonal mismatch offset {}: maxdiff {:.3e} (idx {}) | res {:20.16f} res2 {:20.16f}", offset, maxdiff, idxd, res[idxd],
//     res2[idxd]);
//     // }
//     return res;
// }

template<typename T>
void MatVecMPOS<T>::FactorOP() {
    throw except::runtime_error("MatVecMPOS<T>::FactorOP(): not implemented");
}

template<typename T>
void MatVecMPOS<T>::MultOPv([[maybe_unused]] T *mps_in_, [[maybe_unused]] T *mps_out_, [[maybe_unused]] T eval) {
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
    if(mpos_A.size() == 1) {
        tools::common::contraction::matrix_vector_product(mps_out, mps_in, mpos_A.front(), envL_A, envR_A);
    } else {
        tools::common::contraction::matrix_vector_product(mps_out, mps_in, mpos_A_shf, envL_A, envR_A);
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
void MatVecMPOS<T>::MultBx(T *mps_in_, T *mps_out_) {
    if(mpos_B.empty() and mpos_B_shf.empty()) return;
    auto token = t_multAx->tic_token();

    auto mps_in  = Eigen::TensorMap<Eigen::Tensor<T, 3>>(mps_in_, shape_mps);
    auto mps_out = Eigen::TensorMap<Eigen::Tensor<T, 3>>(mps_out_, shape_mps);
    auto tmp_mps = Eigen::Tensor<T, 3>(shape_mps);
    if(mpos_B.size() == 1) {
        // tools::common::contraction::matrix_vector_product(tmp_mps, mps_in, mpos_B.front(), envL_B, envR_B);
        tools::common::contraction::matrix_vector_product(mps_out, mps_in, mpos_B.front(), envL_B, envR_B);
    } else if(!mpos_B_shf.empty()) {
        // tools::common::contraction::matrix_vector_product(tmp_mps, mps_in, mpos_B_shf, envL_B, envR_B);
        tools::common::contraction::matrix_vector_product(mps_out, mps_in, mpos_B_shf, envL_B, envR_B);
    }
    // if(mpos_B.size() == 1) {
    //     tools::common::contraction::matrix_vector_product(tmp_mps, mps_in, mpos_B.front(), envL_B, envR_B);
    //     tools::common::contraction::matrix_vector_product(mps_out, tmp_mps, mpos_B.front(), envL_B, envR_B);
    // } else if(!mpos_B_shf().empty()) {
    //     tools::common::contraction::matrix_vector_product(tmp_mps, mps_in, mpos_B_shf, envL_B, envR_B);
    //     tools::common::contraction::matrix_vector_product(mps_out, tmp_mps, mpos_B_shf, envL_B, envR_B);
    // }
    num_mv++;
}

template<typename T>
void MatVecMPOS<T>::MultBx(void *x, int *ldx, void *y, int *ldy, int *blockSize, [[maybe_unused]] primme_params *primme, [[maybe_unused]] int *err) {
    for(int i = 0; i < *blockSize; i++) {
        T *mps_in_ptr  = static_cast<T *>(x) + *ldx * i;
        T *mps_out_ptr = static_cast<T *>(y) + *ldy * i;
        MultBx(mps_in_ptr, mps_out_ptr);
    }
    *err = 0;
}

// template<typename T>
// std::vector<std::pair<T, long>> MatVecMPOS<T>::get_k_smallest(const VectorType &vec, size_t k) const {
//     using idx_pair = std::pair<Scalar, long>;
//     std::priority_queue<T> pq;
//     for(auto d : vec) {
//         if(pq.size() >= k && pq.top() > d) {
//             pq.push(d);
//             pq.pop();
//         } else if(pq.size() < k) {
//             pq.push(d);
//         }
//     }
//     Scalar                kth_element = pq.top();
//     std::vector<idx_pair> result;
//     for(long i = 0; i < vec.size(); i++)
//         if(vec[i] <= kth_element) { result.emplace_back(idx_pair{vec[i], i}); }
//     return result;
// }
template<typename T>
std::vector<long> MatVecMPOS<T>::get_k_smallest(const VectorType &vec, size_t k) const {
    std::priority_queue<T> pq;
    for(auto d : vec) {
        if(pq.size() >= k && pq.top() > d) {
            pq.push(d);
            pq.pop();
        } else if(pq.size() < k) {
            pq.push(d);
        }
    }
    Scalar            kth_element = pq.top();
    std::vector<long> result;
    for(long i = 0; i < vec.size(); i++)
        if(vec[i] <= kth_element) { result.emplace_back(i); }
    return result;
}

template<typename T>
std::vector<long> MatVecMPOS<T>::get_k_largest(const VectorType &vec, size_t k) const {
    using idx_pair = std::pair<Scalar, long>;
    std::priority_queue<idx_pair, std::vector<idx_pair>, std::greater<idx_pair>> q;
    for(long i = 0; i < vec.size(); ++i) {
        if(q.size() < k)
            q.emplace(vec[i], i);
        else if(q.top().first < vec[i]) {
            q.pop();
            q.emplace(vec[i], i);
        }
    }
    k = q.size();
    std::vector<long> res(k);
    for(size_t i = 0; i < k; ++i) {
        res[k - i - 1] = q.top().second;
        q.pop();
    }
    return res;
}
template<typename Derived>
Eigen::Matrix<double, Eigen::Dynamic, 1> cond(const Eigen::MatrixBase<Derived> &m) {
    auto solver = Eigen::BDCSVD(m.eval());
    auto rank   = solver.nonzeroSingularValues();
    return solver.singularValues().head(rank);
}

template<typename T>
void MatVecMPOS<T>::CalcPc(T shift) {
    if(readyCalcPc) return;
    if(preconditioner == eig::Preconditioner::NONE) {
        eig::log->info("MatVecMPOS<T>::CalcPc(): no preconditioner chosen");
        return;
    }
    // auto nnz = get_non_zeros();
    // tools::log->info("nonzeros: {} / {} = {:.16f}", nnz, size_mps * size_mps, static_cast<double>(nnz) / static_cast<double>(size_mps * size_mps));

    // long jcbBlockSize = std::min(jcbMaxBlockSize, shape_mps[0] * shape_mps[1]);
    long jcbBlockSize = jcbMaxBlockSize;
    if(jcbBlockSize == 1) {
        eig::log->trace("MatVecMPOS<T>::CalcPc(): calculating the jacobi preconditioner ... ");
        jcbDiagA = get_diagonal_new(0, mpos_A, envL_A, envR_A);
        jcbDiagB = get_diagonal_new(0, mpos_B, envL_B, envR_B);
        eig::log->debug("MatVecMPOS<T>::CalcPc(): calculating the jacobi preconditioner ... done");
    } else if(jcbBlockSize > 1) {
        // if(jcbBlockSize < shape_mps[0]) jcbBlockSize = std::min(shape_mps[0], 2048l); // Make sure we cover index 0, otherwise the blocks may not be psd
        long nblocks = 1 + ((size_mps - 1) / jcbBlockSize); // ceil: note that the last block may be smaller than blocksize!
        eig::log->debug("MatVecMPOS<T>::CalcPc(): calculating the block jacobi preconditioner | diagonal blocksize {} | nblocks {} ...", jcbBlockSize, nblocks);
        std::vector<double> sparsity;
#pragma omp parallel for ordered schedule(dynamic, 1)
        for(long blkidx = 0; blkidx < nblocks; ++blkidx) {
            long offset = blkidx * jcbBlockSize;
            long extent = std::min((blkidx + 1) * jcbBlockSize - offset, size_mps - offset);

            MatrixType blockA = get_diagonal_block(offset, extent, mpos_A, envL_A, envR_A);
            MatrixType blockB = get_diagonal_block(offset, extent, mpos_B, envL_B, envR_B);
            MatrixType blockI = blockA - blockB * shift; // This is what we invert

            bool lltSuccess  = false;
            bool ldltSuccess = false;
            bool luSuccess   = false;
            bool qrSuccess   = false;
            switch(factorization) {
                case eig::Factorization::NONE: {
                    eig::log->warn("MatvecMPOS::CalcPc(): No factorization has been set for the preconditioner");
                    break;
                }
                case eig::Factorization::LLT: {
                    if(auto llt = Eigen::LLT<MatrixType, Eigen::Lower>(blockI); llt.info() == Eigen::Success) {
                        blockI     = llt.solve(MatrixType::Identity(extent, extent));
                        lltSuccess = llt.info() == Eigen::Success;
                        if(lltSuccess) break;
                        blockI = blockA - blockB * shift; // Reset
                        eig::log->info("llt solve failed on block {}/{}", blkidx, nblocks);
                        [[fallthrough]];
                    }
                }
                case eig::Factorization::LDLT: {
                    if(auto ldlt = Eigen::LDLT<MatrixType, Eigen::Lower>(blockI); ldlt.info() == Eigen::Success) {
                        blockI      = ldlt.solve(MatrixType::Identity(extent, extent));
                        ldltSuccess = ldlt.info() == Eigen::Success;
                        if(ldltSuccess) break;
                        blockI = blockA - blockB * shift; // Reset
                        eig::log->info("ldlt solve failed on block {}/{}", blkidx, nblocks);
                        [[fallthrough]];
                    }
                }
                case eig::Factorization::LU: {
                    auto lu   = Eigen::PartialPivLU<MatrixType>(blockI);
                    blockI    = lu.solve(MatrixType::Identity(extent, extent));
                    luSuccess = true;
                    break;
                }
                case eig::Factorization::QR: {
                    auto qr   = Eigen::HouseholderQR<MatrixType>(blockI);
                    blockI    = qr.solve(MatrixType::Identity(extent, extent));
                    qrSuccess = true;
                    break;
                }
            }

            if(!lltSuccess and !ldltSuccess and !luSuccess and !qrSuccess) {
                eig::log->warn("factorization {} (and others) failed on block {}/{} ... resorting to Eigen::ColPivHouseholderQR",
                               eig::FactorizationToString(factorization), blkidx, nblocks);
                blockI = Eigen::ColPivHouseholderQR<MatrixType>(blockI).inverse(); // Should work on any matrix
            }
#pragma omp ordered
            {
                denseJcbBlocks.emplace_back(offset, blockI); //
                double sp = static_cast<double>(blockI.cwiseAbs().count()) / static_cast<double>(blockI.size());
                sparsity.emplace_back(sp);
            }
        }
        auto spavg = std::accumulate(sparsity.begin(), sparsity.end(), 0.0) / static_cast<double>(sparsity.size());
        eig::log->debug(
            "MatVecMPOS<T>::CalcPc(): calculating the block jacobi preconditioner | diagonal blocksize {} | nblocks {} ... done (avg sparsity {:.3e})",
            jcbBlockSize, nblocks, spavg);
    }

    readyCalcPc = true;
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
    auto mps_in  = Eigen::Map<VectorType>(mps_in_, size_mps);
    auto mps_out = Eigen::Map<VectorType>(mps_out_, size_mps);

    if(jcbMaxBlockSize == 1) {
        if(not readyCalcPc) { CalcPc(0.0); }
        auto token = t_multPc->tic_token();

        // Diagonal jacobi preconditioner
        denseJcbDiagonal = ((jcbDiagA.array()).cwiseInverse().cwiseProduct(jcbDiagB.array())).matrix();
        // Use a small cutoff for small differences following Eq 2.11 in https://sdm.lbl.gov/~kewu/ps/thesis.pdf
        //         diagshf                   = (diagonal.array() - shift).cwiseAbs();
        static constexpr auto eps = std::numeric_limits<double>::epsilon();
        for(auto &d : denseJcbDiagonal) {
            if(d >= 0 and d <= +eps) d = +eps;
            if(d <= 0 and d >= -eps) d = -eps;
        }
        mps_out = denseJcbDiagonal.array().cwiseProduct(mps_in.array());
        num_pc++;
    } else if(jcbMaxBlockSize > 1) {
        if(not readyCalcPc) { CalcPc(shift); }
        auto token = t_multPc->tic_token();

#pragma omp parallel for
        for(size_t idx = 0; idx < denseJcbBlocks.size(); ++idx) {
            const auto &[offset, block]               = denseJcbBlocks[idx];
            long extent                               = block.rows();
            mps_out.segment(offset, extent).noalias() = block.template selfadjointView<Eigen::Lower>() * mps_in.segment(offset, extent);
        }
#pragma omp parallel for
        for(size_t idx = 0; idx < sparseJcbBlocks.size(); ++idx) {
            const auto &[offset, block]               = sparseJcbBlocks[idx];
            long extent                               = block.rows();
            mps_out.segment(offset, extent).noalias() = block.template selfadjointView<Eigen::Lower>() * mps_in.segment(offset, extent);
        }

        num_pc++;
    }
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
    auto shift_per_mpo = shift / static_cast<double>(mpos_A.size());
    auto sigma_per_mpo = sigma / static_cast<double>(mpos_A.size());
    for(size_t idx = 0; idx < mpos_A.size(); ++idx) {
        // The MPO is a rank4 tensor ijkl where the first 2 ij indices draw a simple
        // rank2 matrix, where each element is also a matrix with the size
        // determined by the last 2 indices kl.
        // When we shift an MPO, all we do is subtract a diagonal matrix from
        // the botton left corner of the ij-matrix.
        auto &mpo  = mpos_A[idx];
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
    if(not mpos_A_shf.empty()) {
        mpos_A_shf.clear();
        for(const auto &mpo : mpos_A) mpos_A_shf.emplace_back(mpo.shuffle(tenx::array4{2, 3, 0, 1}));
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
void MatVecMPOS<T>::set_jcbMaxBlockSize(std::optional<long> size) {
    if(size.has_value()) jcbMaxBlockSize = std::clamp(size.value(), 1l, size_mps);
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
    return mpos_A;
}
template<typename T>
const Eigen::Tensor<T, 3> &MatVecMPOS<T>::get_envL() const {
    return envL_A;
}
template<typename T>
const Eigen::Tensor<T, 3> &MatVecMPOS<T>::get_envR() const {
    return envR_A;
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
    for(const auto &mpo : mpos_A) shapes.emplace_back(mpo.dimensions());
    return shapes;
}

template<typename T>
std::array<long, 3> MatVecMPOS<T>::get_shape_envL() const {
    return envL_A.dimensions();
}
template<typename T>
std::array<long, 3> MatVecMPOS<T>::get_shape_envR() const {
    return envR_A.dimensions();
}

template<typename T>
Eigen::Tensor<T, 6> MatVecMPOS<T>::get_tensor() const {
    if(mpos_A.size() == 1) {
        auto t_token = t_genMat->tic_token();
        eig::log->debug("Generating tensor");

        auto                d0      = shape_mps[0];
        auto                d1      = shape_mps[1];
        auto                d2      = shape_mps[2];
        auto               &threads = tenx::threads::get();
        Eigen::Tensor<T, 6> tensor;
        tensor.resize(tenx::array6{d0, d1, d2, d0, d1, d2});
        tensor.device(*threads->dev) =
            envL_A.contract(mpos_A.front(), tenx::idx({2}, {0})).contract(envR_A, tenx::idx({2}, {2})).shuffle(tenx::array6{2, 0, 4, 3, 1, 5});

        return tensor;
    }
    throw std::runtime_error("MatVecMPOS<T>::get_tensor(): Not implemented for mpos.size() > 1");
}
template<typename T>
Eigen::Tensor<T, 6> MatVecMPOS<T>::get_tensor_shf() const {
    Eigen::Tensor<T, 6> tensor = get_tensor();
    return tensor.shuffle(tenx::array6{1, 2, 0, 4, 5, 3});
}

template<typename T>
Eigen::Tensor<T, 6> MatVecMPOS<T>::get_tensor_ene() const {
    if(mpos_B.size() == 1) {
        auto t_token = t_genMat->tic_token();
        eig::log->debug("Generating tensor");

        auto                d0      = shape_mps[0];
        auto                d1      = shape_mps[1];
        auto                d2      = shape_mps[2];
        auto               &threads = tenx::threads::get();
        Eigen::Tensor<T, 6> tensor;
        tensor.resize(tenx::array6{d0, d1, d2, d0, d1, d2});
        tensor.device(*threads->dev) =
            envL_B.contract(mpos_B.front(), tenx::idx({2}, {0})).contract(envR_B, tenx::idx({2}, {2})).shuffle(tenx::array6{2, 0, 4, 3, 1, 5});
        return tensor;
    }
    throw std::runtime_error("MatVecMPOS<T>::get_tensor_ene(): Not implemented for mpos.size() > 1");
}

template<typename T>
typename MatVecMPOS<T>::MatrixType MatVecMPOS<T>::get_matrix() const {
    return tenx::MatrixCast(get_tensor(), rows(), cols());
}
template<typename T>
typename MatVecMPOS<T>::MatrixType MatVecMPOS<T>::get_matrix_shf() const {
    return tenx::MatrixCast(get_tensor_shf(), rows(), cols());
}
template<typename T>
typename MatVecMPOS<T>::MatrixType MatVecMPOS<T>::get_matrix_ene() const {
    return tenx::MatrixCast(get_tensor_ene(), rows(), cols());
}

template<typename T>
typename MatVecMPOS<T>::SparseType MatVecMPOS<T>::get_sparse_matrix() const {
    // Fill lower
    std::vector<Eigen::Triplet<T, long>> trip;
    trip.reserve(size_mps);
    // #pragma omp parallel for collapse(2)
    for(long J = 0; J < size_mps; J++) {
#pragma omp parallel for
        for(long I = J; I < size_mps; I++) {
            auto elem = get_matrix_element(I, J, mpos_A, envL_A, envR_A);
            if(std::abs(elem) > std::numeric_limits<double>::epsilon()) {
#pragma omp critical
                { trip.emplace_back(Eigen::Triplet<T, long>{I, J, elem}); }
            }
        }
    }
    SparseType sparseMatrix(size_mps, size_mps);
    sparseMatrix.setFromTriplets(trip.begin(), trip.end());
    return sparseMatrix;
}

template<typename T>
double MatVecMPOS<T>::get_sparsity() const {
    auto sp = get_sparse_matrix();
    auto n  = static_cast<double>(size_mps);
    return (static_cast<double>(sp.nonZeros()) * 2.0 - n) / (n * n);
}
template<typename T>
long MatVecMPOS<T>::get_non_zeros() const {
    auto sp = get_sparse_matrix();
    return sp.nonZeros();
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
