// #define DMRG_ENABLE_TBLIS
#include "matvec_mpos.h"
#include "../log.h"
#include "config/settings.h"
#include "debug/info.h"
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
#include <unsupported/Eigen/src/IterativeSolvers/IncompleteLU.h>

// #if defined(DMRG_ENABLE_TBLIS)
//     #include <tblis/tblis.h>
//     #include <tblis/util/thread.h>
//     #include <tci/tci_config.h>
// #endif
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
        if constexpr(std::is_same_v<T, cx64>) {
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
        if constexpr(std::is_same_v<T, cx64>) {
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

template MatVecMPOS<cx64>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvEne &> &envs_);
template MatVecMPOS<fp64>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvEne &> &envs_);
template MatVecMPOS<cx64>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvVar &> &envs_);
template MatVecMPOS<fp64>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvVar &> &envs_);

template<typename T>
template<typename EnvTypeA, typename EnvTypeB>
MatVecMPOS<T>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, /*!< The Hamiltonian MPO's  */
                          const env_pair<const EnvTypeA &>                         &enva_, /*!< The left and right environments.  */
                          const env_pair<const EnvTypeB &>                         &envb_)
    : MatVecMPOS(mpos_, enva_) {
    // static_assert(sfinae::is_any_v<EnvTypeA, EnvVar>);
    // static_assert(sfinae::is_any_v<EnvTypeB, EnvEne>);
    if constexpr(std::is_same_v<EnvTypeB, EnvEne>) {
        if constexpr(std::is_same_v<T, cx64>) {
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
        if constexpr(std::is_same_v<T, cx64>) {
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

template MatVecMPOS<fp64>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvVar &> &enva_,
                                      const env_pair<const EnvEne &> &envb);
template MatVecMPOS<cx64>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvVar &> &enva_,
                                      const env_pair<const EnvEne &> &envb);

template MatVecMPOS<fp64>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvEne &> &enva_,
                                      const env_pair<const EnvVar &> &envb);
template MatVecMPOS<cx64>::MatVecMPOS(const std::vector<std::reference_wrapper<const MpoSite>> &mpos_, const env_pair<const EnvEne &> &enva_,
                                      const env_pair<const EnvVar &> &envb);

template<typename T>
int MatVecMPOS<T>::rows() const {
    return safe_cast<int>(size_mps);
}

template<typename T>
int MatVecMPOS<T>::cols() const {
    return safe_cast<int>(size_mps);
}

std::vector<long> get_offset(long x, long rank, long size) {
    std::vector<long> indices(rank);
    for(int i = 0; i < rank; i++) {
        indices[i] = x % size;
        x /= size;
    }
    return indices;
}

template<size_t rank>
constexpr std::array<long, rank> get_offset(long flatindex, const std::array<long, rank> &dimensions) {
    std::array<long, rank> indices;
    for(size_t i = 0; i < rank; i++) {
        indices[i] = flatindex % (dimensions[i]);
        flatindex /= (dimensions[i]);
    }
    return indices;
}

std::vector<long> get_offset(long flatindex, size_t rank, const std::vector<long> &dimensions) {
    std::vector<long> indices(rank);
    for(size_t i = 0; i < rank; i++) {
        indices[i] = flatindex % (dimensions[i]);
        flatindex /= (dimensions[i]);
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
        rowindices = get_offset(I, shape_mps);
        colindices = rowindices;
    } else {
        rowindices = get_offset(I, shape_mps);
        colindices = get_offset(J, shape_mps);
    }
    auto ir = rowindices[0];
    auto jr = rowindices[1];
    auto kr = rowindices[2];
    auto ic = colindices[0];
    auto jc = colindices[1];
    auto kc = colindices[2];
    // Index i is special, since it is composed of all the mpo physical indices.
    // auto idxs  = get_offset(i, mpos.size(), mpos.front().dimension(2)); // Maps i to tensor indices to select the physical indices in the
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

    auto irxs  = get_offset(ir, MPOS.size(), spindims); // Maps ir to tensor indices to select the upper physical indices in the mpos
    auto icxs  = get_offset(ic, MPOS.size(), spindims); // Maps ic to tensor indices to select the lower physical indices in the mpos
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
    // eig::log->info("({:3}, {:3}) = [{:3} {:3} {:3}]  [{:3} {:3} {:3}] = {:.16f}", I, J, ir, jr, kr, ic, jc, kc, result);

    // return result;
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

template<auto rank>
constexpr long ravel_multi_index(const std::array<long, rank> &multi_index, const std::array<long, rank> &dimensions, char order = 'F') noexcept {
    assert(order == 'F' or order == 'C');
    if(order == 'F') {
        long index   = 0;
        long dimprod = 1;
        for(size_t i = 0; i < rank; ++i) {
            index += multi_index[i] * dimprod;
            dimprod *= dimensions[i];
        }
        return index;
    }
    if(order == 'C') {
        long index = 0;
        for(size_t i = 0; i < rank; ++i) index = index * dimensions[i] + multi_index[i];
        return index;
    }

    return -1;
}

template<auto rank>
constexpr std::array<long, rank> get_extent(long N, const std::array<long, rank> &dimensions, const std::array<long, rank> &offsets = {0}) {
    // Finds the largest subindex extents of a linear index, guaranteed to be i
    long offset   = ravel_multi_index(offsets, dimensions);
    long maxcount = std::reduce(dimensions.begin(), dimensions.end(), 1l, std::multiplies<long>());
    assert(N + offset <= maxcount);
    // if(N == maxcount) { return dimensions; }
    if(N + offset > maxcount) throw except::logic_error("N ({}) is outside of bounds for dimensions {}", N, dimensions);
    std::array<long, rank> extents;
    extents.fill(1);
    for(size_t i = 0; i < rank; i++) {
        long count  = std::reduce(extents.begin(), extents.end(), 1l, std::multiplies<long>());
        long newdim = std::min(static_cast<long>(std::ceil(static_cast<double>(N) / static_cast<double>(count))), dimensions[i]);
        assert(newdim >= 0);
        long newcount = count * newdim;
        if(newcount + offset <= maxcount) extents[i] = newdim;
        // extents[i] = std::min(static_cast<long>(std::ceil(static_cast<double>(N) / static_cast<double>(count))), (dimensions[i] - offsets[i]));
    }
    return extents;
}

template<auto rank>
constexpr std::array<long, rank> get_extent(const std::array<long, rank> &I0, const std::array<long, rank> &IN, const std::array<long, rank> &dimensions) {
    // Finds the largest subindex extents of a linear index, guaranteed to be i
    auto extent = dimensions;
    if(I0 == IN) {
        extent.fill(1);
        return extent;
    }
    auto INN = get_offset(1l + ravel_multi_index(IN, dimensions), dimensions);
    if(INN == std::array<long, rank>{0}) INN = dimensions;
    for(size_t idx = rank - 1; idx < rank; --idx) {
        extent[idx] = std::clamp(INN[idx] - I0[idx] + 1, 1l, dimensions[idx] - I0[idx]);
        if(INN[idx] > I0[idx]) break;
    }
    return extent;
}

long round_dn(long num, long multiple) { return (num / multiple) * multiple; }

long round_up(long num, long multiple) {
    if(multiple == 0) return num;
    long remainder = num % multiple;
    if(remainder == 0) return num;
    return num + multiple - remainder;
}

template<typename T>
typename MatVecMPOS<T>::MatrixType MatVecMPOS<T>::get_diagonal_block_old(long offset, long extent, const std::vector<Eigen::Tensor<T, 4>> &MPOS,
                                                                         const Eigen::Tensor<T, 3> &ENVL, const Eigen::Tensor<T, 3> &ENVR) const {
    if(MPOS.empty()) return MatrixType::Identity(extent, extent);
    extent = std::min(extent, size_mps - offset);
    if(offset >= size_mps) { return {}; }
    auto res = MatrixType(extent, extent);
#pragma omp parallel for collapse(2)
    for(long J = 0; J < extent; J++) {
        for(long I = J; I < extent; I++) {
            if(I + offset >= size_mps) continue;
            if(J + offset >= size_mps) continue;
            auto elem = get_matrix_element(I + offset, J + offset, MPOS, ENVL, ENVR);
            res(I, J) = elem;
            if constexpr(std::is_same_v<T, cx64>)
                res(J, I) = std::conj(elem);
            else
                res(J, I) = elem;
        }
    }
    return res;
}

template<typename T>
typename MatVecMPOS<T>::MatrixType MatVecMPOS<T>::get_diagonal_block(long offset, long extent, const std::vector<Eigen::Tensor<T, 4>> &MPOS,
                                                                     const Eigen::Tensor<T, 3> &ENVL, const Eigen::Tensor<T, 3> &ENVR) const {
    if(MPOS.empty()) return MatrixType::Identity(extent, extent);
    extent = std::min(extent, size_mps - offset);
    if(offset >= size_mps) { return {}; }
    auto t_old = tid::ur("old");
    auto t_new = tid::ur("new");

    if(MPOS.size() > 1) {
        auto res = MatrixType(extent, extent);
        t_old.tic();
#pragma omp parallel for collapse(2)
        for(long J = 0; J < extent; J++) {
            for(long I = J; I < extent; I++) {
                if(I + offset >= size_mps) continue;
                if(J + offset >= size_mps) continue;
                auto elem = get_matrix_element(I + offset, J + offset, MPOS, ENVL, ENVR);
                res(I, J) = elem;
                if constexpr(std::is_same_v<T, cx64>)
                    res(J, I) = std::conj(elem);
                else
                    res(J, I) = elem;
            }
        }
        t_old.toc();
        return res;
    } else {
        // t_new.tic();
        MatrixType resC(extent, extent);
        resC.setZero();
        long J0 = offset;
        long JN = offset + extent - 1;
        long JY = J0;
        while(JY <= JN) {
            long I0 = JY;
            long IN = offset + extent - 1;
            long R0 = std::clamp(round_dn(I0, shape_mps[0] * shape_mps[1]), 0l, size_mps - 1);
            long RN = std::clamp(round_up(IN, shape_mps[0] * shape_mps[1]), 0l, size_mps - 1);
            long C0 = JY;
            long CN = JY;
            if(R0 == C0) {
                // We are at a tensor dimension boundary, so we can actually calculate more columns when we start from here
                // There is no point in taking too many columns though, since we discard the top right triangle of the sub block
                CN = std::clamp(round_up(C0 + 1, shape_mps[0] * shape_mps[1]) - 1, C0, JN);
            }
            auto R0_ijk = get_offset(R0, shape_mps);
            auto RN_ijk = get_offset(RN, shape_mps);
            auto C0_ijk = get_offset(C0, shape_mps);
            auto CN_ijk = get_offset(CN, shape_mps);
            auto R_ext  = get_extent(R0_ijk, RN_ijk, shape_mps);
            auto C_ext  = get_extent(C0_ijk, CN_ijk, shape_mps);

            std::array<long, 2> ext_blk2 = {R_ext[0] * R_ext[1] * R_ext[2], C_ext[0] * C_ext[1] * C_ext[2]};
            std::array<long, 3> off_envl = {R0_ijk[1], C0_ijk[1], 0};
            std::array<long, 3> ext_envl = {R_ext[1], C_ext[1], MPOS.front().dimension(0)};
            std::array<long, 3> off_envr = {R0_ijk[2], C0_ijk[2], 0};
            std::array<long, 3> ext_envr = {R_ext[2], C_ext[2], MPOS.front().dimension(1)};
            std::array<long, 4> off_mpos = {off_envl[2], off_envr[2], R0_ijk[0], C0_ijk[0]};
            std::array<long, 4> ext_mpos = {ext_envl[2], ext_envr[2], R_ext[0], C_ext[0]};

            std::array<long, 2> off_resC = {I0 - offset, JY - offset};
            std::array<long, 2> ext_resC = {IN - I0 + 1, CN - C0 + 1};

            auto  blockC  = Eigen::Tensor<T, 2>(ext_blk2);
            auto &threads = tenx::threads::get();

            blockC.device(*threads->dev) = ENVL.slice(off_envl, ext_envl)
                                               .contract(MPOS.front().slice(off_mpos, ext_mpos), tenx::idx({2}, {0}))
                                               .contract(ENVR.slice(off_envr, ext_envr), tenx::idx({2}, {2}))
                                               .shuffle(tenx::array6{2, 0, 4, 3, 1, 5})
                                               .reshape(ext_blk2);

            resC.block(off_resC[0], off_resC[1], ext_resC[0], ext_resC[1]) =
                Eigen::Map<MatrixType>(blockC.data(), blockC.dimension(0), blockC.dimension(1)).block(I0 - R0, JY - C0, ext_resC[0], ext_resC[1]);
            // MatrixType matC = resC.block(off_resC[0], off_resC[1], ext_resC[0], ext_resC[1]);
            // MatrixType matF = res.block(off_resC[0], off_resC[1], ext_resC[0], ext_resC[1]);
            // if(!matC.isApprox(matF)) {
            //     eig::log->info("I0       {}", I0);
            //     eig::log->info("J0       {}", J0);
            //     eig::log->info("IN       {}", IN);
            //     eig::log->info("JN       {}", JN);
            //     eig::log->info("R0       {}", R0);
            //     eig::log->info("RN       {}", RN);
            //     eig::log->info("C0       {}", C0);
            //     eig::log->info("CN       {}", CN);
            //     eig::log->info("R0_ijk   {}", R0_ijk);
            //     eig::log->info("RN_ijk   {}", RN_ijk);
            //     eig::log->info("C0_ijk   {}", C0_ijk);
            //     eig::log->info("CN_ijk   {}", CN_ijk);
            //     eig::log->info("R_ext    {}", R_ext);
            //     eig::log->info("C_ext    {}", C_ext);
            //     eig::log->info("ext_blk6 {}", ext_blk6);
            //     eig::log->info("ext_blk2 {}", ext_blk2);
            //     eig::log->info("off_envl {}", off_envl);
            //     eig::log->info("ext_envl {}", ext_envl);
            //     eig::log->info("off_envr {}", off_envr);
            //     eig::log->info("ext_envr {}", ext_envr);
            //     eig::log->info("off_mpos {}", off_mpos);
            //     eig::log->info("ext_mpos {}", ext_mpos);
            //     eig::log->info("off_resC {}", off_resC);
            //     eig::log->info("ext_resC {}", ext_resC);
            //     for(long r = 0; r < matC.rows(); ++r) {
            //         VectorType vecC = matC.row(r);
            //         VectorType vecF = matF.row(r);
            //         eig::log->info("({},{}:{}):  {} | {} ", r, JY, JY + CN - C0 + 1, vecC, vecF);
            //     }
            //     throw except::logic_error("matC and matF mismatch");
            // }
            JY += (CN - C0 + 1);
        }
        return resC.template selfadjointView<Eigen::Lower>();
    }
}

template<typename T>
typename MatVecMPOS<T>::MatrixType MatVecMPOS<T>::get_diagonal_block(long offset, long extent, T shift, const std::vector<Eigen::Tensor<T, 4>> &MPOS_A,
                                                                     const Eigen::Tensor<T, 3> &ENVL_A, const Eigen::Tensor<T, 3> &ENVR_A,
                                                                     const std::vector<Eigen::Tensor<T, 4>> &MPOS_B, const Eigen::Tensor<T, 3> &ENVL_B,
                                                                     const Eigen::Tensor<T, 3> &ENVR_B) const {
    if(MPOS_A.empty()) return MatrixType::Identity(extent, extent);
    extent = std::min(extent, size_mps - offset);
    if(offset >= size_mps) { return {}; }
    auto t_old = tid::ur("old");
    auto t_new = tid::ur("new");
    auto res   = MatrixType(extent, extent);

    if(MPOS_A.size() > 1) {
        t_old.tic();
#pragma omp parallel for collapse(2)
        for(long J = 0; J < extent; J++) {
            for(long I = J; I < extent; I++) {
                if(I + offset >= size_mps) continue;
                if(J + offset >= size_mps) continue;
                // res.template selfadjointView<Eigen::Lower>()(I, J) = get_matrix_element(I + offset, J + offset); // Lower part is sufficient
                auto elemA = get_matrix_element(I + offset, J + offset, MPOS_A, ENVL_A, ENVR_A);
                auto elemB = get_matrix_element(I + offset, J + offset, MPOS_B, ENVL_B, ENVR_B);
                auto elem  = elemA - shift * elemB;
                res(I, J)  = elem;
                if constexpr(std::is_same_v<T, cx64>)
                    res(J, I) = std::conj(elem);
                else
                    res(J, I) = elem;
            }
        }
        return res;
        t_old.toc();
    }

    {
        // auto dbg = MatrixType(extent, extent);
        // #pragma omp parallel for collapse(2)
        // for(long J = 0; J < extent; J++) {
        //     for(long I = J; I < extent; I++) {
        //         if(I + offset >= size_mps) continue;
        //         if(J + offset >= size_mps) continue;
        //         // res.template selfadjointView<Eigen::Lower>()(I, J) = get_matrix_element(I + offset, J + offset); // Lower part is sufficient
        //         auto elemA = get_matrix_element(I + offset, J + offset, MPOS_A, ENVL_A, ENVR_A);
        //         auto elemB = get_matrix_element(I + offset, J + offset, MPOS_B, ENVL_B, ENVR_B);
        //         auto elem  = elemA - shift * elemB;
        //         dbg(I, J)  = elem;
        //         if constexpr(std::is_same_v<T, cx64>)
        //             dbg(J, I) = std::conj(elem);
        //         else
        //             dbg(J, I) = elem;
        //     }
        // }

        // t_new.tic();
        res.setZero();
        long J0 = offset;
        long JN = offset + extent - 1;
        long JY = J0;
        while(JY <= JN) {
            long I0 = JY;
            long IN = offset + extent - 1;
            long R0 = std::clamp(round_dn(I0, shape_mps[0] * shape_mps[1]), 0l, size_mps - 1);
            long RN = std::clamp(round_up(IN, shape_mps[0] * shape_mps[1]), 0l, size_mps - 1);
            long C0 = JY;
            long CN = JY;
            if(R0 == C0) {
                // We are at a tensor dimension boundary, so we can actually calculate more columns when we start from here
                // There is no point in taking too many columns though, since we discard the top right triangle of the sub block
                CN = std::clamp(round_up(C0 + 1, shape_mps[0] * shape_mps[1]) - 1, C0, JN);
            }

            auto R0_ijk   = get_offset(R0, shape_mps);
            auto RN_ijk   = get_offset(RN, shape_mps);
            auto C0_ijk   = get_offset(C0, shape_mps);
            auto CN_ijk   = get_offset(CN, shape_mps);
            auto R_ext    = get_extent(R0_ijk, RN_ijk, shape_mps);
            auto C_ext    = get_extent(C0_ijk, CN_ijk, shape_mps);
            auto ext_blk2 = std::array<long, 2>{R_ext[0] * R_ext[1] * R_ext[2], C_ext[0] * C_ext[1] * C_ext[2]};
            auto off_res  = std::array<long, 2>{I0 - offset, JY - offset};
            auto ext_res  = std::array<long, 2>{IN - I0 + 1, CN - C0 + 1};

            auto get_tile2 = [&](const std::vector<Eigen::Tensor<T, 4>> &MPOS, const Eigen::Tensor<T, 3> &ENVL, const Eigen::Tensor<T, 3> &ENVR) -> MatrixType {
                if(MPOS.empty()) {
                    auto block = MatrixType(ext_res[0], ext_res[1]);
                    block.setZero();
                    // Set 1's along the diagonal of the supermatrix
                    for(long R = 0; R < block.rows(); ++R) {
                        for(long C = 0; C < block.cols(); ++C) {
                            long I = I0 - R;
                            long J = JY - C;
                            if(I == J) block(R, C) = 1;
                        }
                    }

                    return block;
                }
                std::array<long, 3> off_envl = {R0_ijk[1], C0_ijk[1], 0};
                std::array<long, 3> ext_envl = {R_ext[1], C_ext[1], MPOS.front().dimension(0)};
                std::array<long, 3> off_envr = {R0_ijk[2], C0_ijk[2], 0};
                std::array<long, 3> ext_envr = {R_ext[2], C_ext[2], MPOS.front().dimension(1)};
                std::array<long, 4> off_mpos = {off_envl[2], off_envr[2], R0_ijk[0], C0_ijk[0]};
                std::array<long, 4> ext_mpos = {ext_envl[2], ext_envr[2], R_ext[0], C_ext[0]};

                auto  block   = Eigen::Tensor<T, 2>(ext_blk2);
                auto &threads = tenx::threads::get();

                block.device(*threads->dev) = ENVL.slice(off_envl, ext_envl)
                                                  .contract(MPOS.front().slice(off_mpos, ext_mpos), tenx::idx({2}, {0}))
                                                  .contract(ENVR.slice(off_envr, ext_envr), tenx::idx({2}, {2}))
                                                  .shuffle(tenx::array6{2, 0, 4, 3, 1, 5})
                                                  .reshape(ext_blk2);
                return Eigen::Map<MatrixType>(block.data(), ext_blk2[0], ext_blk2[1]).block(I0 - R0, JY - C0, ext_res[0], ext_res[1]);
            };

            res.block(off_res[0], off_res[1], ext_res[0], ext_res[1]) = get_tile2(MPOS_A, ENVL_A, ENVR_A) - get_tile2(MPOS_B, ENVL_B, ENVR_B) * shift;

            // MatrixType blkres = res.block(off_res[0], off_res[1], ext_res[0], ext_res[1]);
            // MatrixType blkdbg = dbg.block(off_res[0], off_res[1], ext_res[0], ext_res[1]);
            // if(!blkres.isApprox(blkdbg)) {
            //     eig::log->info("I0       {}", I0);
            //     eig::log->info("J0       {}", J0);
            //     eig::log->info("IN       {}", IN);
            //     eig::log->info("JN       {}", JN);
            //     eig::log->info("R0       {}", R0);
            //     eig::log->info("RN       {}", RN);
            //     eig::log->info("C0       {}", C0);
            //     eig::log->info("CN       {}", CN);
            //     eig::log->info("R0_ijk   {}", R0_ijk);
            //     eig::log->info("RN_ijk   {}", RN_ijk);
            //     eig::log->info("C0_ijk   {}", C0_ijk);
            //     eig::log->info("CN_ijk   {}", CN_ijk);
            //     eig::log->info("R_ext    {}", R_ext);
            //     eig::log->info("C_ext    {}", C_ext);
            //     eig::log->info("ext_blk2 {}", ext_blk2);
            //     // eig::log->info("ext_blk6 {}", ext_blk6);
            //     // eig::log->info("off_envl {}", off_envl);
            //     // eig::log->info("ext_envl {}", ext_envl);
            //     // eig::log->info("off_envr {}", off_envr);
            //     // eig::log->info("ext_envr {}", ext_envr);
            //     // eig::log->info("off_mpos {}", off_mpos);
            //     // eig::log->info("ext_mpos {}", ext_mpos);
            //     // eig::log->info("off_resC {}", off_resC);
            //     // eig::log->info("ext_resC {}", ext_resC);
            //     for(long r = 0; r < blkres.rows(); ++r) {
            //         VectorType vecC = blkres.row(r);
            //         VectorType vecF = blkdbg.row(r);
            //         eig::log->info("({},{}:{}):  {} | {} ", r, JY, JY + CN - C0 + 1, vecC, vecF);
            //     }
            //     throw except::logic_error("matC and matF mismatch");
            // }
            JY += (CN - C0 + 1);
        }

        // t_new.toc();
        return res.template selfadjointView<Eigen::Lower>();
    }
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
//             rowindices = get_offset(diagidx, shape_mps);
//             colindices = rowindices;
//         } else {
//             rowindices = get_offset(diagidx, shape_mps);
//             colindices = get_offset(diagidx + offset, shape_mps);
//         }
//         auto ir = rowindices[0];
//         auto jr = rowindices[1];
//         auto kr = rowindices[2];
//         auto ic = colindices[0];
//         auto jc = colindices[1];
//         auto kc = colindices[2];
//
//         // Index i is special, since it is composed of all the mpo physical indices.
//         // auto idxs  = get_offset(i, mpos.size(), mpos.front().dimension(2)); // Maps i to tensor indices to select the physical indices in the
//         // mpos
//         auto irxs   = get_offset(ir, mpos_A.size(), spindims); // Maps ir to tensor indices to select the upper phsical indices in the mpos
//         auto icxs   = get_offset(ic, mpos_A.size(), spindims); // Maps ic to tensor indices to select the lower phsical indices in the mpos
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

// template<typename T>
// void MatVecMPOS<T>::FactorOP() {
//     throw except::runtime_error("MatVecMPOS<T>::FactorOP(): not implemented");
// }

template<typename T>
void MatVecMPOS<T>::FactorOP() {
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
void MatVecMPOS<T>::MultOPv([[maybe_unused]] T *mps_in_, [[maybe_unused]] T *mps_out_) {
    throw except::runtime_error("void MatVecMPOS<T>::MultOPv(...) not implemented");
}

// template<typename T>
// void MatVecMPOS<T>::MultOPv([[maybe_unused]] void *x, [[maybe_unused]] int *ldx, [[maybe_unused]] void *y, [[maybe_unused]] int *ldy,
//                             [[maybe_unused]] int *blockSize, [[maybe_unused]] primme_params *primme, [[maybe_unused]] int *err) {
//     throw except::runtime_error("MatVecMPOS<T>::MultOPv(...): not implemented");
// }
template<typename T>
void MatVecMPOS<T>::MultOPv(void *x, int *ldx, void *y, int *ldy, int *blockSize, [[maybe_unused]] primme_params *primme, int *err) {
    auto token = t_multOPv->tic_token();
    switch(side) {
        case eig::Side::R: {
            for(int i = 0; i < *blockSize; i++) {
                T *x_ptr = static_cast<T *>(x) + *ldx * i;
                T *y_ptr = static_cast<T *>(y) + *ldy * i;
                if(factorization == eig::Factorization::LDLT) {
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
                } else {
                    throw except::runtime_error("Invalid factorization: {}", eig::FactorizationToString(factorization));
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

// template<typename T>
// std::vector<long> MatVecMPOS<T>::get_k_smallest(const VectorType &vec, size_t k) const {
//     std::priority_queue<T> pq;
//     for(auto d : vec) {
//         if(pq.size() >= k && pq.top() > d) {
//             pq.push(d);
//             pq.pop();
//         } else if(pq.size() < k) {
//             pq.push(d);
//         }
//     }
//     Scalar            kth_element = pq.top();
//     std::vector<long> result;
//     for(long i = 0; i < vec.size(); i++)
//         if(vec[i] <= kth_element) { result.emplace_back(i); }
//     return result;
// }
//
// template<typename T>
// std::vector<long> MatVecMPOS<T>::get_k_largest(const VectorType &vec, size_t k) const {
//     using idx_pair = std::pair<Scalar, long>;
//     std::priority_queue<idx_pair, std::vector<idx_pair>, std::greater<idx_pair>> q;
//     for(long i = 0; i < vec.size(); ++i) {
//         if(q.size() < k)
//             q.emplace(vec[i], i);
//         else if(q.top().first < vec[i]) {
//             q.pop();
//             q.emplace(vec[i], i);
//         }
//     }
//     k = q.size();
//     std::vector<long> res(k);
//     for(size_t i = 0; i < k; ++i) {
//         res[k - i - 1] = q.top().second;
//         q.pop();
//     }
//     return res;
// }

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
        eig::log->trace("MatVecMPOS<T>::CalcPc(): calculating the jacobi preconditioner ... (shift = {:.16f})", shift);
        jcbDiagA = get_diagonal_new(0, mpos_A, envL_A, envR_A);
        jcbDiagB = get_diagonal_new(0, mpos_B, envL_B, envR_B);
        eig::log->debug("MatVecMPOS<T>::CalcPc(): calculating the jacobi preconditioner ... done (shift = {:.16f})", shift);
    } else if(jcbBlockSize > 1) {
        long nblocks = 1 + ((size_mps - 1) / jcbBlockSize); // ceil: note that the last block may be smaller than blocksize!

        eig::log->debug("MatVecMPOS<T>::CalcPc(): calculating the block jacobi preconditioner | size {} | diagonal blocksize {} | nblocks {} ...", size_mps,
                        jcbBlockSize, nblocks);
        std::vector<double> sparsity;
        auto                m_rss = debug::mem_hwm_in_mb();
        auto                t_jcb = tid::ur("jcb");
        t_jcb.tic();

#pragma omp parallel for ordered schedule(dynamic, 1)
        for(long blkidx = 0; blkidx < nblocks; ++blkidx) {
            long offset = blkidx * jcbBlockSize;
            long extent = std::min((blkidx + 1) * jcbBlockSize - offset, size_mps - offset);
            eig::log->trace("calculating block {}/{}", blkidx, nblocks);
            MatrixType blockI = get_diagonal_block(offset, extent, shift, mpos_A, envL_A, envR_A, mpos_B, envL_B, envR_B);
            double     sp     = static_cast<double>(blockI.cwiseAbs().count()) / static_cast<double>(blockI.size());

            // auto       blockA    = get_diagonal_block_old(offset, extent, mpos_A, envL_A, envR_A);
            // auto       blockB    = get_diagonal_block_old(offset, extent, mpos_B, envL_B, envR_B);
            // MatrixType blockIdbg = blockA - blockB * shift; // Reset
            //
            // if(!blockI.isApprox(blockIdbg, 1e-12)) { throw except::runtime_error("blockI != blockIdbg"); }

            eig::log->trace("calculating block {}/{} ... done", blkidx, nblocks);
            auto bicg         = std::make_unique<BICGType>();
            auto cg           = std::make_unique<CGType>();
            auto llt          = std::make_unique<LLTType>();
            auto ldlt         = std::make_unique<LDLTType>();
            auto lu           = std::make_unique<LUType>();
            auto sparseRM     = std::make_unique<SparseRowM>();
            auto sparseCM     = std::make_unique<SparseType>();
            auto blockPtr     = std::make_unique<MatrixType>();
            bool lltSuccess   = false;
            bool ldltSuccess  = false;
            bool luSuccess    = false;
            bool qrSuccess    = false;
            bool bicgSuccess  = false;
            bool cgSuccess    = false;
            bool solveSuccess = false;
            switch(factorization) {
                case eig::Factorization::NONE: {
                    eig::log->warn("MatvecMPOS::CalcPc(): No factorization has been set for the preconditioner");
                    break;
                }
                case eig::Factorization::LLT: {
                    eig::log->trace("llt factorizing block {}/{} ... ", blkidx, nblocks);
                    if(std::real(blockI(0, 0)) > 0.0) {
                        llt->compute(blockI);
                    } else {
                        llt->compute(-blockI); // Can sometimes be negative definite in the generalized problem
                    }
                    lltSuccess = llt->info() == Eigen::Success;
                    // if(lltSuccess) {
                    //     blockI = llt->solve(MatrixType::Identity(extent, extent));
                    //     solveSuccess = llt->info() == Eigen::Success;
                    //     if(!solveSuccess) {
                    //         blockI = get_diagonal_block(offset, extent, shift, mpos_A, envL_A, envR_A, mpos_B, envL_B, envR_B); // Reset...
                    //     }
                    // }
                    eig::log->trace("llt factorizing block {}/{} ... info {}", blkidx, nblocks, static_cast<int>(llt->info()));

                    // if(!mpos_B.empty()) {
                    //     auto sol = Eigen::SelfAdjointEigenSolver<MatrixType>(blockI, Eigen::EigenvaluesOnly);
                    //     eig::log->info("llt factorization block {}/{} ... info {} | diag min {:.16f} max {:.16f} | eig min {:.16f} max {:.16f}",
                    //                    blkidx, nblocks, static_cast<int>(llt->info()), blockI.diagonal().minCoeff(), blockI.diagonal().maxCoeff(),
                    //                    sol.eigenvalues().minCoeff(), sol.eigenvalues().maxCoeff());
                    // }
                    if(lltSuccess) break;

                    // blockI = get_diagonal_block(offset, extent, shift, mpos_A, envL_A, envR_A, mpos_B, envL_B, envR_B);
                    // auto sol = Eigen::SelfAdjointEigenSolver<MatrixType>(blockI, Eigen::EigenvaluesOnly);
                    // eig::log->info(
                    //     "llt factorization failed on block {}/{} ... info {} | diag min {:.16f} max {:.16f} | eig min {:.16f} max {:.16f} | trying ldlt",
                    //     blkidx, nblocks, static_cast<int>(llt->info()), blockI.diagonal().real().minCoeff(), blockI.diagonal().real().maxCoeff(),
                    //     sol.eigenvalues().minCoeff(), sol.eigenvalues().maxCoeff());
                    [[fallthrough]];

                    // if(auto llt = Eigen::LLT<MatrixType, Eigen::Lower>(blockI); llt.info() == Eigen::Success) {
                    //     blockI     = llt.solve(MatrixType::Identity(extent, extent));
                    //     lltSuccess = llt.info() == Eigen::Success;
                    //     if(lltSuccess) break;
                    //     blockI = blockA - blockB * shift; // Reset
                    //     eig::log->info("llt solve failed on block {}/{}", blkidx, nblocks);
                    // }
                    // [[fallthrough]];
                }
                case eig::Factorization::LDLT: {
                    eig::log->trace("ldlt factorizing block {}/{}", blkidx, nblocks);

                    ldlt->compute(blockI);
                    ldltSuccess = ldlt->info() == Eigen::Success;
                    if(ldltSuccess) break;
                    // blockI = blockA - blockB * shift; // Reset
                    // blockI = get_diagonal_block(offset, extent, shift, mpos_A, envL_A, envR_A, mpos_B, envL_B, envR_B);

                    eig::log->debug("ldlt factorization failed on block {}/{}: info {}", blkidx, nblocks, static_cast<int>(ldlt->info()));
                    [[fallthrough]];
                    // if(auto ldlt = Eigen::LDLT<MatrixType, Eigen::Lower>(blockI); ldlt.info() == Eigen::Success) {
                    //     blockI      = ldlt.solve(MatrixType::Identity(extent, extent));
                    //     ldltSuccess = ldlt.info() == Eigen::Success;
                    //     if(ldltSuccess) break;
                    //     blockI = blockA - blockB * shift; // Reset
                    //     eig::log->info("ldlt solve failed on block {}/{}", blkidx, nblocks);
                    // }
                    // [[fallthrough]];
                }
                case eig::Factorization::LU: {
                    eig::log->trace("lu factorizing block {}/{}", blkidx, nblocks);
                    lu->compute(blockI);
                    luSuccess = true;
                    break;

                    // for(auto cidx = 0; cidx < blockI.cols(); ++cidx) {
                    //     VectorType col   = blockI.col(cidx).cwiseAbs();
                    //     blockI.col(cidx) = (col.array() < 2.0 * col.mean()).select(0, blockI.col(cidx)).eval();
                    // }
                    // blockI = blockI.template selfadjointView<Eigen::Lower>();

                    // auto lu   = Eigen::PartialPivLU<MatrixType>(blockI);
                    // blockI    = lu.solve(MatrixType::Identity(extent, extent));
                    // luSuccess = true;
                    // break;
                }
                case eig::Factorization::QR: {
                    auto qr   = Eigen::HouseholderQR<MatrixType>(blockI);
                    blockI    = qr.solve(MatrixType::Identity(extent, extent));
                    qrSuccess = true;
                    break;
                }
                case eig::Factorization::ILUT: {
                    double mean       = blockI.cwiseAbs().mean();
                    double droptol    = 1e-2; // Drops elements smaller than droptol*rowwise().cwiseAbs().mean()
                    int    fillfactor = 10;   // if the original matrix has nnz nonzeros, LU matrix has at most nnz * fillfactor nonseros
                    // sparseI = blockI.sparseView(mean, 1e-6);
                    *sparseRM = blockI.sparseView(mean, 1e-3);
                    sparseRM->makeCompressed();
                    eig::log->trace("bf sparseI block {}/{}: nnz: {:.3e}", blkidx, nblocks,
                                    static_cast<double>(sparseRM->nonZeros()) / static_cast<double>(sparseRM->size()));
                    bicg->setMaxIterations(1000);
                    bicg->setTolerance(1e-3);
                    bicg->preconditioner().setDroptol(droptol);
                    bicg->preconditioner().setFillfactor(fillfactor);
                    bicg->compute(*sparseRM);
                    bicgSuccess = bicg->info() == Eigen::Success;
                    eig::log->trace("bicg compute block {}/{}: info {}", blkidx, nblocks, static_cast<int>(bicg->info()));
                    // if(ilutSuccess) {
                    //     sparseI       = bicg.solve(id);
                    //     // ilutSuccess = bicg.info() == Eigen::Success;
                    //     // sparseI.makeCompressed();
                    // eig::log->info("af sparse block {}/{}: {:.3e}", blkidx, nblocks, static_cast<double>(sparseI.nonZeros()) / sparseI.size());
                    // }
                    if(!bicgSuccess) {
                        eig::log->info("cg solve failed on block {}/{}: {}", blkidx, nblocks, static_cast<int>(bicg->info()));
                        // Take the diagonal of blockI instead
                        sparseRM->resize(0, 0);
                        blockI = MatrixType(blockI.diagonal().cwiseInverse().asDiagonal());
                        break;
                    }
                    break;
                }
                case eig::Factorization::ILDLT: {
                    double mean = blockI.cwiseAbs().mean();
                    // double droptol    = 1e-2; // Drops elements smaller than droptol*rowwise().cwiseAbs().mean()
                    // int    fillfactor = 10;   // if the original matrix has nnz nonzeros, LDLT matrix has at most nnz * fillfactor nonseros
                    *sparseCM = blockI.sparseView(mean, 0.5);
                    // sparseCM = blockI.sparseView();
                    sparseCM->makeCompressed();
                    // *blockPtr = blockI;
                    eig::log->info("bf sparseI block {}/{}: nnz: {:.3e}", blkidx, nblocks,
                                   static_cast<double>(sparseCM->nonZeros()) / static_cast<double>(sparseCM->size()));
                    // eig::log->info("cg compute block {}/{} ", blkidx, nblocks);
                    // cg->setMaxIterations(100);
                    cg->setTolerance(1e-3);
                    // cg->preconditioner().setInitialShift(droptol);
                    // cg->preconditioner().setFillfactor(fillfactor);
                    cg->compute(*sparseCM);
                    cgSuccess = cg->info() == Eigen::Success;
                    eig::log->trace("cg compute block {}/{}: info {}", blkidx, nblocks, static_cast<int>(cg->info()));
                    if(!cgSuccess) {
                        eig::log->info("cg solve failed on block {}/{}: {}", blkidx, nblocks, static_cast<int>(cg->info()));
                        // Take the diagonal of blockI instead
                        sparseCM->resize(0, 0);
                        blockI = MatrixType(blockI.diagonal().cwiseInverse().asDiagonal());
                        break;
                    }
                    break;
                }
            }

            if(!lltSuccess and !ldltSuccess and !luSuccess and !qrSuccess and factorization != eig::Factorization::ILUT and
               factorization != eig::Factorization::ILDLT) {
                eig::log->warn("factorization {} (and others) failed on block {}/{} ... resorting to Eigen::ColPivHouseholderQR",
                               eig::FactorizationToString(factorization), blkidx, nblocks);
                blockI       = Eigen::ColPivHouseholderQR<MatrixType>(blockI).inverse(); // Should work on any matrix
                solveSuccess = true;
            }
#pragma omp ordered
            {
                sparsity.emplace_back(sp);

                if(solveSuccess) {
                    denseJcbBlocks.emplace_back(offset, blockI); //
                } else if(lltSuccess and llt->rows() > 0) {
                    lltJcbBlocks.emplace_back(offset, std::move(llt));
                } else if(ldltSuccess and ldlt->rows() > 0) {
                    ldltJcbBlocks.emplace_back(offset, std::move(ldlt));
                } else if(luSuccess and lu->rows() > 0) {
                    luJcbBlocks.emplace_back(offset, std::move(lu));
                } else if(bicgSuccess and bicg->rows() > 0) {
                    bicgstabJcbBlocks.emplace_back(offset, std::move(sparseRM), std::move(bicg));
                } else if(cgSuccess and cg->rows() > 0) {
                    cgJcbBlocks.emplace_back(offset, std::move(sparseCM), std::move(cg));
                } else if(sparseRM->size() > 0) {
                    sparseJcbBlocks.emplace_back(offset, *sparseRM);
                }
            }
        }
        t_jcb.toc();
        auto spavg = std::accumulate(sparsity.begin(), sparsity.end(), 0.0) / static_cast<double>(sparsity.size());
        eig::log->debug(
            "MatVecMPOS<T>::CalcPc(): calculating the block jacobi preconditioner | size {} | diagonal blocksize {} | nblocks {} ... done | t {:.3e} s | avg "
            "sparsity {:.3e} | mem +{:.3e} MB",
            size_mps, jcbBlockSize, nblocks, t_jcb.get_last_interval(), spavg, debug::mem_hwm_in_mb() - m_rss);
    }

    readyCalcPc = true;
}

template<typename T>
void MatVecMPOS<T>::MultPc([[maybe_unused]] void *x, [[maybe_unused]] int *ldx, [[maybe_unused]] void *y, [[maybe_unused]] int *ldy,
                           [[maybe_unused]] int *blockSize, [[maybe_unused]] primme_params *primme, [[maybe_unused]] int *err) {
    for(int i = 0; i < *blockSize; i++) {
        T shift = primme->ShiftsForPreconditioner[i];
        if(!mpos_B.empty())
            shift = std::abs(primme->stats.estimateMaxEVal) > std::abs(primme->stats.estimateMinEVal) ? primme->stats.estimateMaxEVal
                                                                                                      : primme->stats.estimateMinEVal;
        T *mps_in_ptr  = static_cast<T *>(x) + *ldx * i;
        T *mps_out_ptr = static_cast<T *>(y) + *ldy * i;
        MultPc(mps_in_ptr, mps_out_ptr, shift);
        // MultPc(mps_in_ptr, mps_out_ptr, primme->ShiftsForPreconditioner[i]);
    }
    *err = 0;
}

template<typename T>
void MatVecMPOS<T>::MultPc([[maybe_unused]] T *mps_in_, [[maybe_unused]] T *mps_out_, T shift) {
    if(preconditioner == eig::Preconditioner::NONE) return;
    auto mps_in  = Eigen::Map<VectorType>(mps_in_, size_mps);
    auto mps_out = Eigen::Map<VectorType>(mps_out_, size_mps);

    if(jcbMaxBlockSize == 1) {
        if(not readyCalcPc) { CalcPc(shift); }
        auto token = t_multPc->tic_token();
        // Diagonal jacobi preconditioner
        denseJcbDiagonal = (jcbDiagA.array() - shift * jcbDiagB.array()).matrix();
        mps_out          = denseJcbDiagonal.array().cwiseInverse().cwiseProduct(mps_in.array());
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
#pragma omp parallel for
        for(size_t idx = 0; idx < lltJcbBlocks.size(); ++idx) {
            const auto &[offset, solver] = lltJcbBlocks[idx];
            long extent                  = solver->rows();
            auto mps_out_segment         = Eigen::Map<VectorType>(mps_out_ + offset, extent);
            auto mps_in_segment          = Eigen::Map<const VectorType>(mps_in_ + offset, extent);
            mps_out_segment.noalias()    = solver->solve(mps_in_segment);
        }
#pragma omp parallel for
        for(size_t idx = 0; idx < ldltJcbBlocks.size(); ++idx) {
            const auto &[offset, solver] = ldltJcbBlocks[idx];
            long extent                  = solver->rows();
            auto mps_out_segment         = Eigen::Map<VectorType>(mps_out_ + offset, extent);
            auto mps_in_segment          = Eigen::Map<const VectorType>(mps_in_ + offset, extent);
            mps_out_segment.noalias()    = solver->solve(mps_in_segment);
        }
#pragma omp parallel for
        for(size_t idx = 0; idx < luJcbBlocks.size(); ++idx) {
            const auto &[offset, solver] = luJcbBlocks[idx];
            long extent                  = solver->rows();
            auto mps_out_segment         = Eigen::Map<VectorType>(mps_out_ + offset, extent);
            auto mps_in_segment          = Eigen::Map<const VectorType>(mps_in_ + offset, extent);
            mps_out_segment.noalias()    = solver->solve(mps_in_segment);
        }
        // #pragma omp parallel for
        for(size_t idx = 0; idx < bicgstabJcbBlocks.size(); ++idx) {
            const auto &[offset, sparseI, solver]     = bicgstabJcbBlocks[idx];
            long extent                               = solver->rows();
            mps_out.segment(offset, extent).noalias() = solver->solveWithGuess(mps_in.segment(offset, extent), mps_out.segment(offset, extent));
        }
        // #pragma omp parallel for
        for(size_t idx = 0; idx < cgJcbBlocks.size(); ++idx) {
            const auto &[offset, sparseI, solver]     = cgJcbBlocks[idx];
            long extent                               = solver->rows();
            mps_out.segment(offset, extent).noalias() = solver->solveWithGuess(mps_in.segment(offset, extent), mps_out.segment(offset, extent));
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
        if constexpr(std::is_same_v<T, fp64>)
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
    if(size.has_value()) {
        // We want the block sizes to be roughly equal, so we reduce the block size until the remainder is zero or larger than 80% of the block size
        // This ensures that the last block isn't too much smaller than the other ones
        jcbMaxBlockSize = std::clamp(size.value(), 1l, size_mps);
        long newsize    = jcbMaxBlockSize;
        long rem        = num::mod(size_mps, newsize);
        auto bestsize   = std::pair<long, long>{rem, newsize};
        while(newsize >= jcbMaxBlockSize / 2) {
            rem = num::mod(size_mps, newsize);
            if(rem > bestsize.first or rem == 0) { bestsize = std::pair<long, long>{rem, newsize}; }
            if(rem == 0) break;              // All equal size
            if(rem > 4 * newsize / 5) break; // The last is at least 80% the size of the others
            newsize -= 2;
        }
        if(bestsize.second != jcbMaxBlockSize) {
            eig::log->trace("Adjusted block size to {}", bestsize.second);
            jcbMaxBlockSize = std::clamp(bestsize.second, jcbMaxBlockSize / 2, size_mps);
        }
        // jcbMaxBlockSize = std::clamp(size.value(), 1l, size_mps);
    }
}

template<typename T>
T MatVecMPOS<T>::get_shift() const {
    if constexpr(std::is_same_v<T, fp64>)
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
    if constexpr(std::is_same_v<T, fp64>)
        return eig::Type::REAL;
    else if constexpr(std::is_same_v<T, cx64>)
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
    trip.reserve(static_cast<size_t>(size_mps));
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
long MatVecMPOS<T>::get_jcbMaxBlockSize() const {
    return jcbMaxBlockSize;
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
