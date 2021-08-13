
#include <Eigen/IterativeLinearSolvers>
#include <io/fmt.h>
#include <iostream>
#include <math/tenx.h>
#include <tools/common/contraction.h>
#include <unsupported/Eigen/IterativeSolvers>

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
    //    eig::Form           form = eig::Form::SYMM;
    //    eig::Side           side = eig::Side::R;
    // Profiling
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
                //                std::cout << "dst size " << dst.size() << " | rhs size " << rhs.size() << std::endl;
                //                if(dst.size() != rhs.size())
                //                    dst.conservativeResize(rhs.size());

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
                throw std::runtime_error(fmt::format("Dimension mismatch mps {} and envL {}", mps.dimensions(), envL.dimensions()));
            if(mps.dimension(2) != envR.dimension(0))
                throw std::runtime_error(fmt::format("Dimension mismatch mps {} and envR {}", mps.dimensions(), envR.dimensions()));
            if(mps.dimension(0) != mpo.dimension(2)) throw std::runtime_error(fmt::format("Dimension mismatch mps {} and mpo {}", mps.dimensions(), mpo.dimensions()));
            if(envL.dimension(2) != mpo.dimension(0))
                throw std::runtime_error(fmt::format("Dimension mismatch envL {} and mpo {}", envL.dimensions(), mpo.dimensions()));
            if(envR.dimension(2) != mpo.dimension(1))
                throw std::runtime_error(fmt::format("Dimension mismatch envR {} and mpo {}", envR.dimensions(), mpo.dimensions()));
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
//            std::cout << "BiCGSTAB: #iterations: " << bicg.iterations()
//                << ", #count: " << matRepl.matvecs
//                << ", estimated error: " << bicg.error()
//                << std::endl;

        }
//        {
//            Eigen::ConjugateGradient<MatrixReplacement<Scalar>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg;
//            cg.setMaxIterations(20000);
//            cg.compute(matRepl);
//            res = cg.solve(mps);
//            std::cout << "CG      : #iterations: " << cg.iterations()
//                      << ", #count: " << matRepl.matvecs
//                      << ", estimated error: " << cg.error()
//                      << std::endl;
//        }
        if constexpr (std::is_same_v<Scalar,double>){
            Eigen::MINRES<MatrixReplacement<Scalar>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> minres;
            minres.setMaxIterations(MaxIters);
            minres.setTolerance(tolerance);
            minres.compute(matRepl);
            res = minres.solve(mps);
//            std::cout << "MINRES:   #iterations: " << minres.iterations()
//                      << ", #count: " << matRepl.matvecs
//                      << ", norm: " << res.norm()
//                      << ", estimated error: " << minres.error()
//                      << std::endl;
    //    std::cout << "x: \n" << x << std::endl;
        }

//        {
//            Eigen::ConjugateGradient<MatrixReplacement<Scalar>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg;
//            cg.setMaxIterations(20000);
//            cg.compute(matRepl);
//            res = cg.solve(mps);
//            std::cout << "CG      : #iterations: " << cg.iterations()
//                      << ", #count: " << matRepl.matvecs
//                      << ", estimated error: " << cg.error()
//                      << std::endl;
//        }




}
using namespace tools::common::contraction;
template void tools::common::contraction::matrix_inverse_vector_product(
    real * res_ptr,
    const real * const mps_ptr, std::array<long,3> mps_dims,
    const real * const mpo_ptr, std::array<long,4> mpo_dims,
    const real * const envL_ptr, std::array<long,3> envL_dims,
    const real * const envR_ptr, std::array<long,3> envR_dims);
template void tools::common::contraction::matrix_inverse_vector_product(
    cplx * res_ptr,
    const cplx * const mps_ptr, std::array<long,3> mps_dims,
    const cplx * const mpo_ptr, std::array<long,4> mpo_dims,
    const cplx * const envL_ptr, std::array<long,3> envL_dims,
    const cplx * const envR_ptr, std::array<long,3> envR_dims);
/* clang-format on */
