
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

    res.device(*Textra::omp::dev) =
             mps
            .contract(envL, Textra::idx({1}, {0}))
            .contract(mpo, Textra::idx({0, 3}, {2, 0}))
            .contract(envR, Textra::idx({0, 2}, {0, 2}))
            .shuffle(Textra::array3{1, 0, 2});

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
