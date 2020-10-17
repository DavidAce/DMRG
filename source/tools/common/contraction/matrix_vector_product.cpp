
#include <general/nmspc_tensor_extra.h>
#include <tools/common/contraction.h>
#include <tools/common/fmt.h>

/* clang-format off */
template<typename Scalar>
void tools::common::contraction::matrix_vector_product(Scalar * res_ptr,
                           const Scalar * const mps_ptr, std::array<long,3> mps_dims,
                           const Scalar * const mpo_ptr, std::array<long,4> mpo_dims,
                           const Scalar * const envL_ptr, std::array<long,3> envL_dims,
                           const Scalar * const envR_ptr, std::array<long,3> envR_dims){

    // This applies the mpo's with corresponding environments to local multisite mps
    // This is usually the operation H|psi>  or H²|psi>
    auto res = Eigen::TensorMap<Eigen::Tensor<Scalar,3>>(res_ptr,mps_dims);
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

    auto dsizes = mps.dimensions();
    auto chiL   = dsizes[1];
    auto chiR   = dsizes[2];

//    res.resize(mps.dimensions());
    if(chiL > chiR) {
        Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(Textra::array3{1, 0, 2});
        res.device(*Textra::omp::dev) =
             mps_shuffled
            .contract(envL, Textra::idx({0}, {0}))
            .contract(mpo, Textra::idx({0, 3}, {2, 0}))
            .contract(envR, Textra::idx({0, 2}, {0, 2}))
            .shuffle(Textra::array3{1, 0, 2});
    } else {
        Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(Textra::array3{2, 0, 1});
        res.device(*Textra::omp::dev) =
             mps_shuffled
            .contract(envR, Textra::idx({0}, {0}))
            .contract(mpo, Textra::idx({0, 3}, {2, 1}))
            .contract(envL, Textra::idx({0, 2}, {0, 2}))
            .shuffle(Textra::array3{1, 2, 0});
    }
}
using namespace tools::common::contraction;
template void tools::common::contraction::matrix_vector_product(
                           real * res_ptr,
                           const real * const mps_ptr, std::array<long,3> mps_dims,
                           const real * const mpo_ptr, std::array<long,4> mpo_dims,
                           const real * const envL_ptr, std::array<long,3> envL_dims,
                           const real * const envR_ptr, std::array<long,3> envR_dims);
template void tools::common::contraction::matrix_vector_product(
                            cplx * res_ptr,
                            const cplx * const mps_ptr, std::array<long,3> mps_dims,
                            const cplx * const mpo_ptr, std::array<long,4> mpo_dims,
                            const cplx * const envL_ptr, std::array<long,3> envL_dims,
                            const cplx * const envR_ptr, std::array<long,3> envR_dims);
/* clang-format on */


//
//template void tools::common::contraction::matrix_vector_product(real *, const real *const, const real *const, const real *const, const real *const);
//template void tools::common::contraction::matrix_vector_product(cplx *, const cplx *const, const cplx *const, const cplx *const, const cplx *const);
//template void tools::common::contraction::matrix_vector_product(real *, const TM<T3<const real>> &, const T4<real> &, const T3<real> &, const T3<real> &);
//template void tools::common::contraction::matrix_vector_product(cplx *, const TM<T3<const cplx>> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);





//template<typename res_type, typename mps_type, typename mpo_type, typename env_type>
//void tools::common::contraction::matrix_vector_product(res_type &result, const mps_type &mps, const mpo_type &mpo, const env_type &envL, const env_type &envR) {
//    // This applies the local hamiltonian and environment to a local multisite mps
//    // This is usually the operation H|psi>  or H²|psi>
//    // Note that the environments must contain the correct type of mpos
//    static_assert(mpo_type::NumIndices == 4 and "Wrong mps tensor rank != 4 passed to calculation of matrix_vector_product");
//    static_assert(mpo_type::NumIndices == 4 and "Wrong mpo tensor rank != 4 passed to calculation of matrix_vector_product");
//    static_assert(env_type::NumIndices == 3 and "Wrong env tensor rank != 3 passed to calculation of matrix_vector_product");
//    if(mps.dimension(1) != envL.dimension(0))
//        throw std::runtime_error(fmt::format("Dimension mismatch mps {} and envL {}", mps.dimensions(), envL.dimensions()));
//    if(mps.dimension(2) != envR.dimension(0))
//        throw std::runtime_error(fmt::format("Dimension mismatch mps {} and envR {}", mps.dimensions(), envR.dimensions()));
//    if(mps.dimension(0) != mpo.dimension(2)) throw std::runtime_error(fmt::format("Dimension mismatch mps {} and mpo {}", mps.dimensions(), mpo.dimensions()));
//    if(envL.dimension(2) != mpo.dimension(0))
//        throw std::runtime_error(fmt::format("Dimension mismatch envL {} and mpo {}", envL.dimensions(), mpo.dimensions()));
//    if(envR.dimension(2) != mpo.dimension(1))
//        throw std::runtime_error(fmt::format("Dimension mismatch envR {} and mpo {}", envR.dimensions(), mpo.dimensions()));
//
//    using Scalar = typename mps_type::Scalar;
//    auto dsizes = mps.dimensions();
//    auto chiL   = dsizes[1];
//    auto chiR   = dsizes[2];
//
//    result.resize(mps.dimensions());
//    if(chiL > chiR) {
//        Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(Textra::array3{1, 0, 2});
//        result.device(*Textra::omp::dev)      = mps_shuffled.contract(envL, Textra::idx({0}, {0}))
//                                               .contract(mpo, Textra::idx({0, 3}, {2, 0}))
//                                               .contract(envR, Textra::idx({0, 2}, {0, 2}))
//                                               .shuffle(Textra::array3{1, 0, 2});
//    } else {
//        Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(Textra::array3{2, 0, 1});
//        result.device(*Textra::omp::dev)      = mps_shuffled.contract(envR, Textra::idx({0}, {0}))
//                                               .contract(mpo, Textra::idx({0, 3}, {2, 1}))
//                                               .contract(envL, Textra::idx({0, 2}, {0, 2}))
//                                               .shuffle(Textra::array3{1, 2, 0});
//    }
//}
//
//template<typename mps_type, typename mpo_type, typename env_type>
//mps_type tools::common::contraction::matrix_vector_product(const mps_type &mps, const mpo_type &mpo, const env_type &envL, const env_type &envR) {
//    mps_type result;
//    matrix_vector_product(result, mps, envL, envR);
//    return result;
//}
//
//using namespace tools::common::contraction;
//
//template T3<real> tools::common::contraction::matrix_vector_product(const T3<real> &, const T4<real> &, const T3<real> &, const T3<real> &);
//template T3<cplx> tools::common::contraction::matrix_vector_product(const T3<cplx> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);
//template void tools::common::contraction::matrix_vector_product(T3<real> &, const T3<real> &, const T4<real> &, const T3<real> &, const T3<real> &);
//template void tools::common::contraction::matrix_vector_product(T3<cplx> &, const T3<cplx> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);
//template void tools::common::contraction::matrix_vector_product(T3<real> &, const TM<T3<const real>> &, const T4<real> &, const T3<real> &, const T3<real> &);
//template void tools::common::contraction::matrix_vector_product(T3<cplx> &, const TM<T3<const cplx>> &, const T4<cplx> &, const T3<cplx> &, const T3<cplx> &);


//
// template<typename Scalar>
// void ceres_direct_functor<Scalar>::get_Hv(const VectorType &v) const {
//    t_vH->tic();
//    auto log2chiL = std::log2(dsizes[1]);
//    auto log2chiR = std::log2(dsizes[2]);
//    //            size_t log2spin  = std::log2(multiComponents.dsizes[0]);
//    if(log2chiL > log2chiR) {
//        if(print_path) tools::log->trace("get_Hv path: log2chiL > log2chiR ");
//
//        Eigen::Tensor<Scalar, 3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(v.derived().data(), dsizes).shuffle(Textra::array3{1, 0, 2});
//        Hv_tensor.device(*Textra::omp::dev) = theta.contract(envL, Textra::idx({0}, {0}))
//            .contract(mpo, Textra::idx({0, 3}, {2, 0}))
//            .contract(envR, Textra::idx({0, 2}, {0, 2}))
//            .shuffle(Textra::array3{1, 0, 2});
//    } else {
//        if(print_path) tools::log->trace("get_Hv path: log2chiL <= log2chiR ");
//
//        Eigen::Tensor<Scalar, 3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(v.derived().data(), dsizes).shuffle(Textra::array3{2, 0, 1});
//        Hv_tensor.device(*Textra::omp::dev) = theta.contract(envR, Textra::idx({0}, {0}))
//            .contract(mpo, Textra::idx({0, 3}, {2, 1}))
//            .contract(envL, Textra::idx({0, 2}, {0, 2}))
//            .shuffle(Textra::array3{1, 2, 0});
//    }
//
//    t_vH->toc();
//}