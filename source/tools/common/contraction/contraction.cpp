#include "../contraction.h"
#include <Eigen/IterativeLinearSolvers>
#include <io/fmt.h>
#include <math/tenx.h>
#include <tid/tid.h>
#include <unsupported/Eigen/IterativeSolvers>

#if defined(DMRG_SAVE_CONTRACTION)
    #include <h5pp/h5pp.h>
#endif
#if defined(DMRG_BENCH_CONTRACTION)
    #include <math/num.h>
    #include <tid/tid.h>
namespace settings {
    constexpr static bool bench_expval = true;
}
#else
namespace settings {
    constexpr static bool bench_expval = false;
}
#endif

using namespace tools::common::contraction;

/* clang-format off */
template<typename Scalar>
double tools::common::contraction::expectation_value(
                           const Scalar * const mps_ptr, std::array<long,3> mps_dims,
                           const Scalar * const mpo_ptr, std::array<long,4> mpo_dims,
                           const Scalar * const envL_ptr, std::array<long,3> envL_dims,
                           const Scalar * const envR_ptr, std::array<long,3> envR_dims){

    std::string bench_suffix;
#if defined(DMRG_BENCH_CONTRACTION)
    if constexpr (settings::bench_expval){
        auto spin = mps_dims[0];
        auto chiL = num::next_multiple<long>(mps_dims[1], 5l);
        auto chiR = num::next_multiple<long>(mps_dims[2], 5l);
        auto mdim = mpo_dims[0];
        bench_suffix = fmt::format("_mps-[{},{},{}]_mpo-[{}]", spin, chiL, chiR, mdim);
    }
#endif
    auto t_expval = tid::tic_token(fmt::format("expval{}",bench_suffix), tid::level::detailed);

    // This measures the expectation value of some multisite mps with respect to some mpo operator and corresponding environments.
    // This is usually the energy E = <psi|H|psi> or variance V = <psi|(H-E)²|psi>
    // Note that the environments must contain the correct type of mpos
    auto mps = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(mps_ptr,mps_dims);
    auto mpo = Eigen::TensorMap<const Eigen::Tensor<const Scalar,4>>(mpo_ptr,mpo_dims);
    auto envL = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(envL_ptr,envL_dims);
    auto envR = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(envR_ptr,envR_dims);

    if(mps.dimension(1) != envL.dimension(0)) throw except::runtime_error("Dimension mismatch mps {} and envL {}", mps.dimensions(), envL.dimensions());
    if(mps.dimension(2) != envR.dimension(0)) throw except::runtime_error("Dimension mismatch mps {} and envR {}", mps.dimensions(), envR.dimensions());
    if(mps.dimension(0) != mpo.dimension(2)) throw except::runtime_error("Dimension mismatch mps {} and mpo {}", mps.dimensions(), mpo.dimensions());
    if(envL.dimension(2) != mpo.dimension(0)) throw except::runtime_error("Dimension mismatch envL {} and mpo {}", envL.dimensions(), mpo.dimensions());
    if(envR.dimension(2) != mpo.dimension(1)) throw except::runtime_error("Dimension mismatch envR {} and mpo {}", envR.dimensions(), mpo.dimensions());

    Eigen::Tensor<Scalar, 0> expval;
    expval.device(tenx::omp::getDevice()) =
        envL
            .contract(mps,             tenx::idx({0}, {1}))
            .contract(mpo,             tenx::idx({2, 1}, {2, 0}))
            .contract(mps.conjugate(), tenx::idx({3, 0}, {0, 1}))
            .contract(envR,            tenx::idx({0, 2, 1}, {0, 1, 2}));

    double moment = 0;
    if constexpr(std::is_same_v<Scalar,cplx>){
        moment = std::real(expval(0));
        if(abs(std::imag(expval(0))) > 1e-10)
            fmt::print("Expectation value has an imaginary part: {:.16f} + i {:.16f}\n", std::real(expval(0)), std::imag(expval(0)));
        //        throw except::runtime_error("Expectation value has an imaginary part: {:.16f} + i {:.16f}", std::real(expval(0)), std::imag(expval(0))));
    }else
        moment = expval(0);
    if(std::isnan(moment) or std::isinf(moment)) throw except::runtime_error("First moment is invalid: {}", moment);

    #if defined(DMRG_SAVE_CONTRACTION)
    {
        t_expval.toc();
        auto file       = h5pp::File("dmrg-contractions.h5", h5pp::FilePermission::READWRITE);
        auto group_num  = 0;
        auto group_name = fmt::format("contraction_{}", group_num);
        while(file.linkExists(group_name)) group_name = fmt::format("contraction_{}", ++group_num);
        file.writeDataset(mps, fmt::format("{}/mps", group_name));
        file.writeDataset(mpo, fmt::format("{}/mpo", group_name));
        file.writeDataset(envL, fmt::format("{}/envL", group_name));
        file.writeDataset(envR, fmt::format("{}/envR", group_name));
    }
    #endif

    return moment;
}

template double tools::common::contraction::expectation_value(const real * const mps_ptr,  std::array<long,3> mps_dims,
                                                              const real * const mpo_ptr,  std::array<long,4> mpo_dims,
                                                              const real * const envL_ptr, std::array<long,3> envL_dims,
                                                              const real * const envR_ptr, std::array<long,3> envR_dims);
template double tools::common::contraction::expectation_value(const cplx * const mps_ptr,  std::array<long,3> mps_dims,
                                                              const cplx * const mpo_ptr,  std::array<long,4> mpo_dims,
                                                              const cplx * const envL_ptr, std::array<long,3> envL_dims,
                                                              const cplx * const envR_ptr, std::array<long,3> envR_dims);

/* clang-format on */

template<typename Scalar_>
class MatrixReplacement;
template<typename T>
using DenseMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

namespace Eigen {
    namespace internal {
        // MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
        template<typename Scalar_>
        struct traits<MatrixReplacement<Scalar_>> : public Eigen::internal::traits<DenseMatrix<Scalar_>> {};
    }
}

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::Matrix.
template<typename Scalar_>
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement<Scalar_>> {
    public:
    // Required typedefs, constants, and method:
    typedef Scalar_ Scalar;
    typedef double  RealScalar;
    typedef int     StorageIndex;
    enum { ColsAtCompileTime = Eigen::Dynamic, MaxColsAtCompileTime = Eigen::Dynamic, IsRowMajor = false };

    const Scalar_      *envL = nullptr;
    const Scalar_      *envR = nullptr;
    const Scalar_      *mpo  = nullptr;
    std::array<long, 3> shape_mps;
    std::array<long, 4> shape_mpo;
    std::array<long, 3> shape_envL;
    std::array<long, 3> shape_envR;
    std::vector<Scalar> shift_mpo;
    long                mps_size;
    // Timers
    mutable int                 counter = 0;
    mutable std::vector<Scalar> tmp;

    [[nodiscard]] Eigen::Index rows() const { return static_cast<int>(mps_size); }; /*!< Linear size\f$d^2 \times \chi_L \times \chi_R \f$  */
    [[nodiscard]] Eigen::Index cols() const { return static_cast<int>(mps_size); }; /*!< Linear size\f$d^2 \times \chi_L \times \chi_R \f$  */

    template<typename Rhs>
    Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs> &x) const {
        return Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
    }

    // Custom API:
    MatrixReplacement() = default;

    void attachTensors(const Scalar_      *envL_,      /*!< The left block tensor.  */
                       const Scalar_      *envR_,      /*!< The right block tensor.  */
                       const Scalar_      *mpo_,       /*!< The Hamiltonian MPO's  */
                       std::array<long, 3> shape_mps_, /*!< An array containing the shapes of the mps  */
                       std::array<long, 4> shape_mpo_  /*!< An array containing the shapes of the mpo  */
    ) {
        envL       = envL_;
        envR       = envR_;
        mpo        = mpo_;
        shape_mps  = shape_mps_;
        shape_mpo  = shape_mpo_;
        shape_envL = {shape_mps_[1], shape_mps_[1], shape_mpo_[0]};
        shape_envR = {shape_mps_[2], shape_mps_[2], shape_mpo_[1]};
        if(envL == nullptr) throw std::runtime_error("Lblock is a nullptr!");
        if(envR == nullptr) throw std::runtime_error("Rblock is a nullptr!");
        if(mpo == nullptr) throw std::runtime_error("mpo is a nullptr!");
        mps_size = shape_mps[0] * shape_mps[1] * shape_mps[2];
        //        t_multAx = std::make_unique<class_tic_toc>(true, 5, "Time MultAx");
    }
};

// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of init::generic_product_impl:
namespace Eigen {
    namespace internal {

        template<typename Rhs, typename ReplScalar>
        struct generic_product_impl<MatrixReplacement<ReplScalar>, Rhs, DenseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
            : generic_product_impl_base<MatrixReplacement<ReplScalar>, Rhs, generic_product_impl<MatrixReplacement<ReplScalar>, Rhs>> {
            typedef typename Product<MatrixReplacement<ReplScalar>, Rhs>::Scalar Scalar;

            template<typename Dest>
            static void scaleAndAddTo(Dest &dst, const MatrixReplacement<ReplScalar> &mat, const Rhs &rhs, const Scalar &alpha) {
                // This method should implement "dst += alpha * lhs * rhs" inplace,
                // however, for iterative solvers, alpha is always equal to 1, so let's not worry about it.
                assert(alpha == Scalar(1) && "scaling is not implemented");
                EIGEN_ONLY_USED_FOR_DEBUG(alpha);

                //                auto token = mat.t_multAx->tic_token();
                mat.tmp.resize(static_cast<size_t>(dst.size()));
                Eigen::Map<Dest> tmp_map(mat.tmp.data(), dst.size());
                tools::common::contraction::matrix_vector_product(tmp_map.data(), rhs.data(), mat.shape_mps, mat.mpo, mat.shape_mpo, mat.envL, mat.shape_envL,
                                                                  mat.envR, mat.shape_envR);

                dst.noalias() += tmp_map;
                mat.counter++;
                //
                //                // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
                //                // but let's do something fancier (and less efficient):
                //                for(Index i=0; i<lhs.cols(); ++i)
                //                    dst += rhs(i) * lhs.my_matrix().col(i);
            }
        };

    }
}

/* clang-format off */
template<typename Scalar>
void tools::common::contraction::matrix_inverse_vector_product(Scalar * res_ptr,
                                                       const Scalar * const mps_ptr, std::array<long,3> mps_dims,
                                                       const Scalar * const mpo_ptr, std::array<long,4> mpo_dims,
                                                       const Scalar * const envL_ptr, std::array<long,3> envL_dims,
                                                       const Scalar * const envR_ptr, std::array<long,3> envR_dims){


        // Here we return x <-- A^-1 * b
        // Where A^-1 * b is obtained by solving
        //       A*x = b
        // using an iterative matrix-free solver.
        {
            auto mps = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(mps_ptr,mps_dims);
            auto mpo = Eigen::TensorMap<const Eigen::Tensor<const Scalar,4>>(mpo_ptr,mpo_dims);
            auto envL = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(envL_ptr,envL_dims);
            auto envR = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(envR_ptr,envR_dims);

            if(mps.dimension(1) != envL.dimension(0))
                throw except::runtime_error("Dimension mismatch mps {} and envL {}", mps.dimensions(), envL.dimensions());
            if(mps.dimension(2) != envR.dimension(0))
                throw except::runtime_error("Dimension mismatch mps {} and envR {}", mps.dimensions(), envR.dimensions());
            if(mps.dimension(0) != mpo.dimension(2)) throw except::runtime_error("Dimension mismatch mps {} and mpo {}", mps.dimensions(), mpo.dimensions());
            if(envL.dimension(2) != mpo.dimension(0))
                throw except::runtime_error("Dimension mismatch envL {} and mpo {}", envL.dimensions(), mpo.dimensions());
            if(envR.dimension(2) != mpo.dimension(1))
                throw except::runtime_error("Dimension mismatch envR {} and mpo {}", envR.dimensions(), mpo.dimensions());
        }

        // Define the "matrix-free" matrix replacement.
        MatrixReplacement<Scalar> matRepl;
        matRepl.attachTensors(envL_ptr, envR_ptr, mpo_ptr, mps_dims, mpo_dims);

        Eigen::Index MaxIters = 200000;
        double tolerance = 1e-14;
        Eigen::Map<tenx::VectorType<Scalar>> res(res_ptr, matRepl.rows());
        Eigen::Map<const tenx::VectorType<Scalar>> mps(mps_ptr, matRepl.rows());
        if constexpr (std::is_same_v<Scalar,std::complex<double>>){

            Eigen::BiCGSTAB<MatrixReplacement<Scalar>, Eigen::IdentityPreconditioner> bicg;
            bicg.compute(matRepl);
            bicg.setMaxIterations(MaxIters);
            bicg.setTolerance(tolerance);
            res = bicg.solve(mps);
        }
        if constexpr (std::is_same_v<Scalar,double>){
            Eigen::MINRES<MatrixReplacement<Scalar>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> minres;
            minres.setMaxIterations(MaxIters);
            minres.setTolerance(tolerance);
            minres.compute(matRepl);
            res = minres.solve(mps);
        }




}

using namespace tools::common::contraction;
template void tools::common::contraction::matrix_inverse_vector_product(
    real * res_ptr,
    const real * const mps_ptr, std::array<long,3> mps_dims,
    const real * const mpo_ptr, std::array<long,4> mpo_dims,
    const real * const envL_ptr, std::array<long,3> envL_dims,
    const real * const envR_ptr, std::array<long,3> envR_dims);
template void tools::common::contraction::matrix_inverse_vector_product(      cplx *       res_ptr,
                                                                        const cplx * const mps_ptr, std::array<long,3> mps_dims,
                                                                        const cplx * const mpo_ptr, std::array<long,4> mpo_dims,
                                                                        const cplx * const envL_ptr, std::array<long,3> envL_dims,
                                                                        const cplx * const envR_ptr, std::array<long,3> envR_dims);
/* clang-format on */

/* clang-format off */
template<typename Scalar>
void tools::common::contraction::matrix_vector_product(      Scalar * res_ptr,
                                                       const Scalar * const mps_ptr, std::array<long,3> mps_dims,
                                                       const Scalar * const mpo_ptr, std::array<long,4> mpo_dims,
                                                       const Scalar * const envL_ptr, std::array<long,3> envL_dims,
                                                       const Scalar * const envR_ptr, std::array<long,3> envR_dims){

    auto t_matvec = tid::tic_token("matrix_vector_product", tid::level::detailed);

    // This applies the mpo's with corresponding environments to local multisite mps
    // This is usually the operation H|psi>  or H²|psi>
    auto res = Eigen::TensorMap<Eigen::Tensor<Scalar,3>>(res_ptr,mps_dims);
    auto mps = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(mps_ptr,mps_dims);
    auto mpo = Eigen::TensorMap<const Eigen::Tensor<const Scalar,4>>(mpo_ptr,mpo_dims);
    auto envL = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(envL_ptr,envL_dims);
    auto envR = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(envR_ptr,envR_dims);

    if(mps.dimension(1) != envL.dimension(0)) throw except::runtime_error("Dimension mismatch mps {} and envL {}", mps.dimensions(), envL.dimensions());
    if(mps.dimension(2) != envR.dimension(0)) throw except::runtime_error("Dimension mismatch mps {} and envR {}", mps.dimensions(), envR.dimensions());
    if(mps.dimension(0) != mpo.dimension(2))  throw except::runtime_error("Dimension mismatch mps {} and mpo {}", mps.dimensions(), mpo.dimensions());
    if(envL.dimension(2) != mpo.dimension(0)) throw except::runtime_error("Dimension mismatch envL {} and mpo {}", envL.dimensions(), mpo.dimensions());
    if(envR.dimension(2) != mpo.dimension(1)) throw except::runtime_error("Dimension mismatch envR {} and mpo {}", envR.dimensions(), mpo.dimensions());

    res.device(tenx::omp::getDevice()) =
             mps
            .contract(envL,     tenx::idx({1}, {0}))
            .contract(mpo,      tenx::idx({0, 3}, {2, 0}))
            .contract(envR,     tenx::idx({0, 2}, {0, 2}))
            .shuffle(tenx::array3{1, 0, 2});

}
using namespace tools::common::contraction;
template void tools::common::contraction::matrix_vector_product(      cplx *       res_ptr,
                                                                const cplx * const mps_ptr, std::array<long,3> mps_dims,
                                                                const cplx * const mpo_ptr, std::array<long,4> mpo_dims,
                                                                const cplx * const envL_ptr, std::array<long,3> envL_dims,
                                                                const cplx * const envR_ptr, std::array<long,3> envR_dims);
template void tools::common::contraction::matrix_vector_product(      real *       res_ptr,
                                                                const real * const mps_ptr, std::array<long,3> mps_dims,
                                                                const real * const mpo_ptr, std::array<long,4> mpo_dims,
                                                                const real * const envL_ptr, std::array<long,3> envL_dims,
                                                                const real * const envR_ptr, std::array<long,3> envR_dims);

template<typename Scalar>
void  tools::common::contraction::contract_mps_bnd(      Scalar * res_ptr      , std::array<long,3> res_dims,
                                                   const Scalar * const mps_ptr, std::array<long,3> mps_dims,
                                                   const Scalar * const bnd_ptr, std::array<long,1> bnd_dims){
    auto t_con = tid::tic_token("contract_mps_bnd", tid::level::detailed);
    auto res = Eigen::TensorMap<Eigen::Tensor<Scalar,3>>(res_ptr,res_dims);
    auto mps = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(mps_ptr,mps_dims);
    auto bnd = Eigen::TensorMap<const Eigen::Tensor<const Scalar,1>>(bnd_ptr,bnd_dims);
    if(mps.dimension(2) != bnd.dimension(0)) throw except::runtime_error("Dimension mismatch mps {} (idx 2) and bnd {} (idx 0)", mps.dimensions(), bnd.dimensions());
    if(mps.dimensions() != res.dimensions()) throw except::runtime_error("Dimension mismatch mps {} and res {}", mps.dimensions(), res.dimensions());
    res.device(tenx::omp::getDevice()) = mps.contract(tenx::asDiagonal(bnd), tenx::idx({2}, {0}));
}
template void tools::common::contraction::contract_mps_bnd(      cplx *       res_ptr, std::array<long,3> res_dims,
                                                           const cplx * const mps_ptr, std::array<long,3> mps_dims,
                                                           const cplx * const bnd_ptr, std::array<long,1> bnd_dims);

template void tools::common::contraction::contract_mps_bnd(      real *       res_ptr, std::array<long,3> res_dims,
                                                           const real * const mps_ptr, std::array<long,3> mps_dims,
                                                           const real * const bnd_ptr, std::array<long,1> bnd_dims);


template<typename Scalar>
void  tools::common::contraction::contract_bnd_mps(
          Scalar * res_ptr      , std::array<long,3> res_dims,
    const Scalar * const bnd_ptr, std::array<long,1> bnd_dims,
    const Scalar * const mps_ptr, std::array<long,3> mps_dims){
    auto t_con = tid::tic_token("contract_bnd_mps", tid::level::detailed);
    auto res = Eigen::TensorMap<Eigen::Tensor<Scalar,3>>(res_ptr,res_dims);
    auto mps = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(mps_ptr,mps_dims);
    auto bnd = Eigen::TensorMap<const Eigen::Tensor<const Scalar,1>>(bnd_ptr,bnd_dims);
    if(mps.dimension(1) != bnd.dimension(0)) throw except::runtime_error("Dimension mismatch mps {} (idx 1) and bnd {} (idx 0)", mps.dimensions(), bnd.dimensions());
    if(mps.dimensions() != res.dimensions()) throw except::runtime_error("Dimension mismatch mps {} and res {}", mps.dimensions(), res.dimensions());
    res.device(tenx::omp::getDevice()) = tenx::asDiagonal(bnd).contract(mps, tenx::idx({1}, {1})).shuffle(tenx::array3{1, 0, 2});
}

template void tools::common::contraction::contract_bnd_mps(      cplx *       res_ptr, std::array<long,3> res_dims,
                                                           const cplx * const bnd_ptr, std::array<long,1> bnd_dims,
                                                           const cplx * const mps_ptr, std::array<long,3> mps_dims);

template void tools::common::contraction::contract_bnd_mps(      real *       res_ptr, std::array<long,3> res_dims,
                                                           const real * const bnd_ptr, std::array<long,1> bnd_dims,
                                                           const real * const mps_ptr, std::array<long,3> mps_dims);


template<typename Scalar>
void tools::common::contraction::contract_mps_mps(      Scalar * res_ptr       , std::array<long,3> res_dims,
                                                  const Scalar * const mpsL_ptr, std::array<long,3> mpsL_dims,
                                                  const Scalar * const mpsR_ptr, std::array<long,3> mpsR_dims){
    auto t_con = tid::tic_token("contract_mps_mps", tid::level::detailed);
    auto res  = Eigen::TensorMap<Eigen::Tensor<Scalar,3>>(res_ptr,res_dims);
    auto mpsL = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(mpsL_ptr, mpsL_dims);
    auto mpsR = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(mpsR_ptr, mpsR_dims);
    constexpr auto shuffle_idx  = std::array<long,4>{0, 2, 1, 3};
    constexpr auto contract_idx = tenx::idx({2}, {1});
    auto check_dims = std::array<long,3>{mpsL.dimension(0) * mpsR.dimension(0), mpsL.dimension(1), mpsR.dimension(2)};
    if(res_dims != check_dims) throw except::runtime_error("res dimension mismatch: dims {} | expected dims {}", res_dims, check_dims);
    if(mpsL.dimension(2) != mpsR.dimension(1)) throw except::runtime_error("Dimension mismatch mpsL {} (idx 2) and mpsR {} (idx 1)", mpsL.dimensions(), mpsR.dimensions());
    res.device(tenx::omp::getDevice()) = mpsL.contract(mpsR, contract_idx).shuffle(shuffle_idx).reshape(res_dims);
}


template void tools::common::contraction::contract_mps_mps(      cplx * res_ptr       , std::array<long,3> res_dims,
                                                           const cplx * const mpsL_ptr, std::array<long,3> mpsL_dims,
                                                           const cplx * const mpsR_ptr, std::array<long,3> mpsR_dims);

//template void tools::common::contraction::contract_mps_mps(      real * res_ptr       , std::array<long,3> res_dims,
//                                                           const real * const mpsL_ptr, std::array<long,3> mpsL_dims,
//                                                           const real * const mpsR_ptr, std::array<long,3> mpsR_dims);



template<typename Scalar>
double tools::common::contraction::contract_mps_mps_overlap(const Scalar * const mps1_ptr, std::array<long,3> mps1_dims,
                                                            const Scalar * const mps2_ptr, std::array<long,3> mps2_dims){
    auto t_con = tid::tic_token("contract_mps_mps_overlap", tid::level::detailed);
    auto mps1 = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(mps1_ptr, mps1_dims);
    auto mps2 = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(mps2_ptr, mps2_dims);
    if(mps1.dimensions() != mps2.dimensions()) throw except::runtime_error("Dimension mismatch mps1 {} and mps2 {}", mps1.dimensions(), mps2.dimensions());
    Eigen::Tensor<Scalar,0> res;
    constexpr auto idxs = tenx::idx({0,1,2},{0,1,2});
    res.device(tenx::omp::getDevice()) = mps1.contract(mps2.conjugate(), idxs);
    return std::abs(res.coeff(0));
}

template double tools::common::contraction::contract_mps_mps_overlap(const cplx * const mps1_ptr, std::array<long,3> mps1_dims,
                                                                     const cplx * const mps2_ptr, std::array<long,3> mps2_dims);
//template double tools::common::contraction::contract_mps_mps_overlap(const real * const mps1_ptr, std::array<long,3> mps1_dims,
//                                                                     const real * const mps2_ptr, std::array<long,3> mps2_dims);


template<typename Scalar>
void tools::common::contraction::contract_mps_mps_partial(      Scalar *       res_ptr , std::array<long,2> res_dims,
                                                          const Scalar * const mps1_ptr, std::array<long,3> mps1_dims,
                                                          const Scalar * const mps2_ptr, std::array<long,3> mps2_dims,
                                                          std::array<long,2> idx){
    auto res  = Eigen::TensorMap<Eigen::Tensor<Scalar,2>>(res_ptr,res_dims);
    auto mps1 = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(mps1_ptr, mps1_dims);
    auto mps2 = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(mps2_ptr, mps2_dims);
    auto idxs = tenx::idx(idx,idx);
    res.device(tenx::omp::getDevice()) = mps1.contract(mps2.conjugate(), idxs);
}

template void tools::common::contraction::contract_mps_mps_partial(      cplx *       res_ptr , std::array<long,2> res_dims,
                                                                   const cplx * const mps1_ptr, std::array<long,3> mps1_dims,
                                                                   const cplx * const mps2_ptr, std::array<long,3> mps2_dims,
                                                                   std::array<long,2> idx);
//template void tools::common::contraction::contract_mps_mps_partial(      real *       res_ptr , std::array<long,2> res_dims,
//                                                                   const real * const mps1_ptr, std::array<long,3> mps1_dims,
//                                                                   const real * const mps2_ptr, std::array<long,3> mps2_dims,
//                                                                   std::array<long,2> idx);


template<typename Scalar>
void tools::common::contraction::contract_env_mps_mpo(      Scalar *      res_ptr, std::array<long, 2> res_dims,
                                                      const Scalar *const env_ptr, std::array<long, 2> env_dims,
                                                      const Scalar *const mps_ptr, std::array<long, 3> mps_dims,
                                                      const Scalar *const mpo_ptr, std::array<long, 2> mpo_dims) {
    auto res = Eigen::TensorMap<Eigen::Tensor<Scalar, 2>>(res_ptr, res_dims);
    auto env = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 2>>(env_ptr, env_dims);
    auto mps = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(mps_ptr, mps_dims);
    auto mpo = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 2>>(mpo_ptr, mpo_dims);
    res.device(tenx::omp::getDevice()) = env.contract(mps,             tenx::idx({0}, {1}))
                                            .contract(mpo,             tenx::idx({1}, {0}))
                                            .contract(mps.conjugate(), tenx::idx({0, 2}, {1, 0}));
}

template void tools::common::contraction::contract_env_mps_mpo(      cplx *       res_ptr , std::array<long,2> res_dims,
                                                               const cplx * const env_ptr , std::array<long,2> env_dims,
                                                               const cplx * const mps_ptr , std::array<long,3> mps_dims,
                                                               const cplx * const mpo_ptr , std::array<long,2> mpo_dims);

template<typename Scalar>
void tools::common::contraction::contract_env_mps_mpo(      Scalar *      res_ptr, std::array<long, 3> res_dims,
                                                      const Scalar *const env_ptr, std::array<long, 3> env_dims,
                                                      const Scalar *const mps_ptr, std::array<long, 3> mps_dims,
                                                      const Scalar *const mpo_ptr, std::array<long, 4> mpo_dims) {
    auto res                           = Eigen::TensorMap<Eigen::Tensor<Scalar, 3>>(res_ptr, res_dims);
    auto env                           = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(env_ptr, env_dims);
    auto mps                           = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(mps_ptr, mps_dims);
    auto mpo                           = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 4>>(mpo_ptr, mpo_dims);
    res.device(tenx::omp::getDevice()) = env.contract(mps,             tenx::idx({0}, {1}))
                                            .contract(mpo,             tenx::idx({1, 2}, {0, 2}))
                                            .contract(mps.conjugate(), tenx::idx({0, 3}, {1, 0}))
                                            .shuffle(                  tenx::array3{0, 2, 1});
}

template void tools::common::contraction::contract_env_mps_mpo(      cplx *       res_ptr , std::array<long,3> res_dims,
                                                               const cplx * const env_ptr , std::array<long,3> env_dims,
                                                               const cplx * const mps_ptr , std::array<long,3> mps_dims,
                                                               const cplx * const mpo_ptr , std::array<long,4> mpo_dims);

template<typename Scalar>
void tools::common::contraction::contract_mps_mpo_env(      Scalar *      res_ptr, std::array<long, 2> res_dims,
                                                      const Scalar *const env_ptr, std::array<long, 2> env_dims,
                                                      const Scalar *const mps_ptr, std::array<long, 3> mps_dims,
                                                      const Scalar *const mpo_ptr, std::array<long, 2> mpo_dims) {
    auto res = Eigen::TensorMap<Eigen::Tensor<Scalar, 2>>(res_ptr, res_dims);
    auto env = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 2>>(env_ptr, env_dims);
    auto mps = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(mps_ptr, mps_dims);
    auto mpo = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 2>>(mpo_ptr, mpo_dims);
    res.device(tenx::omp::getDevice()) =
        env.contract(mps,             tenx::idx({0}, {2}))
           .contract(mpo,             tenx::idx({1}, {0}))
           .contract(mps.conjugate(), tenx::idx({0, 2}, {2, 0}));
}
template void tools::common::contraction::contract_mps_mpo_env(      cplx *       res_ptr , std::array<long,2> res_dims,
                                                               const cplx * const env_ptr , std::array<long,2> env_dims,
                                                               const cplx * const mps_ptr , std::array<long,3> mps_dims,
                                                               const cplx * const mpo_ptr , std::array<long,2> mpo_dims);
template<typename Scalar>
void tools::common::contraction::contract_mps_mpo_env(      Scalar *      res_ptr, std::array<long, 3> res_dims,
                                                      const Scalar *const env_ptr, std::array<long, 3> env_dims,
                                                      const Scalar *const mps_ptr, std::array<long, 3> mps_dims,
                                                      const Scalar *const mpo_ptr, std::array<long, 4> mpo_dims) {
    auto res                           = Eigen::TensorMap<Eigen::Tensor<Scalar, 3>>(res_ptr, res_dims);
    auto env                           = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(env_ptr, env_dims);
    auto mps                           = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(mps_ptr, mps_dims);
    auto mpo                           = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 4>>(mpo_ptr, mpo_dims);
    res.device(tenx::omp::getDevice()) = env.contract(mps,             tenx::idx({0}, {2}))
                                            .contract(mpo,             tenx::idx({1, 2}, {1, 2}))
                                            .contract(mps.conjugate(), tenx::idx({0, 3}, {2, 0}))
                                            .shuffle(                  tenx::array3{0, 2, 1});
}

template void tools::common::contraction::contract_mps_mpo_env(      cplx *       res_ptr , std::array<long,3> res_dims,
                                                               const cplx * const env_ptr , std::array<long,3> env_dims,
                                                               const cplx * const mps_ptr , std::array<long,3> mps_dims,
                                                               const cplx * const mpo_ptr , std::array<long,4> mpo_dims);
/* clang-format on */
