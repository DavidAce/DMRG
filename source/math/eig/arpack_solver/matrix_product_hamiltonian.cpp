#include "matrix_product_hamiltonian.h"
#include <tools/common/contraction.h>
#include <unsupported/Eigen/CXX11/Tensor>

template<typename Scalar_>
MatrixProductHamiltonian<Scalar_>::MatrixProductHamiltonian(const Scalar_ *     envL_,      /*!< The left block tensor.  */
                                                            const Scalar_ *     envR_,      /*!< The right block tensor.  */
                                                            const Scalar_ *     mpo_,       /*!< The Hamiltonian MPO's  */
                                                            std::array<long, 3> shape_mps_, /*!< An array containing the shapes of theta  */
                                                            std::array<long, 4> shape_mpo_  /*!< An array containing the shapes of the MPO  */
                                                            )
    : envL(envL_), envR(envR_), mpo(mpo_), shape_mps(shape_mps_), shape_mpo(shape_mpo_), shape_envL({shape_mps_[1], shape_mps_[1], shape_mpo_[0]}),
      shape_envR({shape_mps_[2], shape_mps_[2], shape_mpo_[1]}) {
    if(envL == nullptr) throw std::runtime_error("Lblock is a nullptr!");
    if(envR == nullptr) throw std::runtime_error("Rblock is a nullptr!");
    if(mpo == nullptr) throw std::runtime_error("mpo is a nullptr!");
    mps_size = shape_mps[0] * shape_mps[1] * shape_mps[2];
    init_profiling();
}

template<typename T>
void MatrixProductHamiltonian<T>::MultAx(T *mps_in_, T *mps_out_) {
    t_multAx->tic();
    tools::common::contraction::matrix_vector_product(mps_out_, mps_in_, shape_mps, mpo, shape_mpo, envL, shape_envL, envR, shape_envR);
    counter++;
    t_multAx->toc();
}

template<typename Scalar>
void MatrixProductHamiltonian<Scalar>::print() const {}

template<typename T>
void MatrixProductHamiltonian<T>::set_shift(std::complex<double> sigma_) {
    if(readyShift) return; // This only happens once!!
    sigma = sigma_;
    // The MPO is a rank4 tensor ijkl where the first 2 ij indices draw a simple
    // rank2 matrix, where each element is also a matrix with the size
    // determined by the last 2 indices kl.
    // When we shift an MPO, all we do is subtract a diagonal matrix from
    // the botton left corner of the ij-matrix.
    auto mpo_size = shape_mpo[0] * shape_mpo[1] * shape_mpo[2] * shape_mpo[3];
    // Setup the start and shifted mpo tensors
    shift_mpo.resize(static_cast<size_t>(mpo_size));
    Eigen::TensorMap<Eigen::Tensor<const T, 4>> start_mpo_map(mpo, shape_mpo);
    Eigen::TensorMap<Eigen::Tensor<T, 4>>       shift_mpo_map(shift_mpo.data(), shape_mpo);

    // Copy the previous mpo into the new one
    shift_mpo_map = start_mpo_map;

    // Setup extents and handy objects
    Eigen::array<long, 4>                            offset4{shape_mpo[0] - 1, 0, 0, 0};
    Eigen::array<long, 4>                            extent4{1, 1, shape_mpo[2], shape_mpo[3]};
    Eigen::array<long, 2>                            extent2{shape_mpo[2], shape_mpo[3]};
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> sigma_Id;
    if constexpr(std::is_same_v<T, eig::real>)
        sigma_Id = std::real(sigma) * Eigen::MatrixXd::Identity(extent2[0], extent2[1]);
    else
        sigma_Id = sigma * Eigen::MatrixXcd::Identity(extent2[0], extent2[1]);
    Eigen::TensorMap<Eigen::Tensor<T, 2>> sigma_Id_map(sigma_Id.data(), sigma_Id.rows(), sigma_Id.cols());

    shift_mpo_map.slice(offset4, extent4).reshape(extent2) -= sigma_Id_map;

    // Reseat the mpo pointer
    mpo        = shift_mpo.data();
    readyShift = true;
}

template<typename Scalar>
void MatrixProductHamiltonian<Scalar>::set_mode(const eig::Form form_) {
    form = form_;
}
template<typename Scalar>
void MatrixProductHamiltonian<Scalar>::set_side(const eig::Side side_) {
    side = side_;
}
template<typename Scalar>
const eig::Form &MatrixProductHamiltonian<Scalar>::get_form() const {
    return form;
}
template<typename Scalar>
const eig::Side &MatrixProductHamiltonian<Scalar>::get_side() const {
    return side;
}

template<typename Scalar>
void MatrixProductHamiltonian<Scalar>::init_profiling() {
    t_multAx = std::make_unique<class_tic_toc>(profile_matrix_product_hamiltonian, 5, "Time MultAx");
}

// Explicit instantiations
template class MatrixProductHamiltonian<double>;
template class MatrixProductHamiltonian<std::complex<double>>;
