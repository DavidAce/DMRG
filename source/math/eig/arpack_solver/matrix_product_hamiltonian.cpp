#include "matrix_product_hamiltonian.h"
#include <tools/common/contraction.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <general/class_tic_toc.h>
#include <math/svd.h>

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
    t_factorOP = std::make_unique<class_tic_toc>(true, 5, "Time FactorOp");
    t_multOPv  = std::make_unique<class_tic_toc>(true, 5, "Time MultOpv");
    t_multAx   = std::make_unique<class_tic_toc>(true, 5, "Time MultAx");
}



template<typename Scalar>
void MatrixProductHamiltonian<Scalar>::FactorOP()
/* We don't actually invert a matrix here: we let an iterative matrix-free solver apply OP^-1 x */
{
    if(readyFactorOp) return; // happens only once
    if(not readyShift) throw std::runtime_error("Cannot FactorOP: Shift value sigma has not been set.");
    t_factorOP->tic();
    readyFactorOp = true;
    t_factorOP->toc();
}


template<typename T>
void MatrixProductHamiltonian<T>::MultAx(T *mps_in_, T *mps_out_) {
    auto token = t_multAx->tic_token();
    tools::common::contraction::matrix_vector_product(mps_out_, mps_in_, shape_mps, mpo, shape_mpo, envL, shape_envL, envR, shape_envR);
    counter++;
}





template<typename T>
void MatrixProductHamiltonian<T>::MultOPv(T *mps_in_, T *mps_out_) {
    t_multOPv->tic();
    switch(side) {
        case eig::Side::R: {
            tools::common::contraction::matrix_inverse_vector_product(mps_out_, mps_in_, shape_mps, mpo, shape_mpo, envL, shape_envL, envR, shape_envR);
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
void MatrixProductHamiltonian<Scalar>::print() const {}

template<typename T>
void MatrixProductHamiltonian<T>::set_shift(std::complex<double> sigma_) {
    if(readyShift) return; // This only happens once!!
    if(readyCompress) throw std::runtime_error("Cannot shift the matrix: it is already compressed!");
    sigma = sigma_;

    // The MPO is a rank4 tensor ijkl where the first 2 ij indices draw a simple
    // rank2 matrix, where each element is also a matrix with the size
    // determined by the last 2 indices kl.
    // When we shift an MPO, all we do is subtract a diagonal matrix from
    // the botton left corner of the ij-matrix.
    auto mpo_size = shape_mpo[0] * shape_mpo[1] * shape_mpo[2] * shape_mpo[3];
    // Setup the start and shifted mpo tensors
    mpo_internal.resize(static_cast<size_t>(mpo_size));
    Eigen::TensorMap<Eigen::Tensor<const T, 4>> start_mpo_map(mpo, shape_mpo);
    Eigen::TensorMap<Eigen::Tensor<T, 4>>       shift_mpo_map(mpo_internal.data(), shape_mpo);

    // Copy the previous mpo into the new one
    shift_mpo_map = start_mpo_map;

    // Setup extents and handy objects
    std::array<long, 4>                              offset4{shape_mpo[0] - 1, 0, 0, 0};
    std::array<long, 4>                              extent4{1, 1, shape_mpo[2], shape_mpo[3]};
    std::array<long, 2>                              extent2{shape_mpo[2], shape_mpo[3]};
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> sigma_Id;
    if constexpr(std::is_same_v<T, eig::real>)
        sigma_Id = std::real(sigma) * Eigen::MatrixXd::Identity(extent2[0], extent2[1]);
    else
        sigma_Id = sigma * Eigen::MatrixXcd::Identity(extent2[0], extent2[1]);
    Eigen::TensorMap<Eigen::Tensor<T, 2>> sigma_Id_map(sigma_Id.data(), sigma_Id.rows(), sigma_Id.cols());

    shift_mpo_map.slice(offset4, extent4).reshape(extent2) -= sigma_Id_map;

    // Reseat the mpo pointer
    mpo        = mpo_internal.data();
    readyShift = true;
}

//template<typename T>
//void MatrixProductHamiltonian<T>::set_shifted_mpo(const T* mpo_, std::array<long, 4> shape_mpo_){
//    if(readyShift) return;
//    mpo = mpo_;
//    shape_mpo = shape_mpo_;
//    readyShift = true;
//}


template<typename T>
void MatrixProductHamiltonian<T>::compress(){
    if(readyCompress) return;

    svd::settings svd_settings;
    svd_settings.use_lapacke = true;
    svd_settings.use_bdc = false;
    svd_settings.threshold = 1e-16;
    svd_settings.switchsize = 4096;
    svd::solver svd(svd_settings);

    Eigen::Tensor<T, 4> mpo_tmp  = Eigen::TensorMap<Eigen::Tensor<const T, 4>> (mpo, shape_mpo);
    Eigen::Tensor<T, 4> mpo_l2r, mpo_r2l;
    {
        // Compress left to right
        auto [U_l2r,S,V_l2r] = svd.split_mpo_l2r(mpo_tmp);
        mpo_l2r = U_l2r.contract(Textra::asDiagonal(S), Textra::idx({1},{0})).shuffle(std::array<long,4>{0,3,1,2});
        // Resize the new buffer for the compressed version of envR
        std::array<long,3> shape_newR = {shape_envR[0], shape_envR[1], S.size()};
        size_t newR_size = shape_newR[0]* shape_newR[1]* shape_newR[2];
        envR_internal.resize(newR_size);

        // Define maps to work with the pointers
        Eigen::TensorMap<Eigen::Tensor<const T, 3>> envR_map(envR, shape_envR);
        Eigen::TensorMap<Eigen::Tensor<T,3>> envR_internal_map (envR_internal.data(), shape_newR);
        envR_internal_map  = envR_map.contract(V_l2r, Textra::idx({2},{1}));

        // Reseat the pointer and update the size
        envR = envR_internal.data();
        shape_envR = shape_newR;
    }
    {
        // Compress right to left
        auto [U_r2l,S ,V_r2l] = svd.split_mpo_r2l(mpo_l2r);
        mpo_r2l = Textra::asDiagonal(S).contract(V_r2l, Textra::idx({1},{0}));

        // Resize the new buffer for the compressed version of envR
        std::array<long,3> shape_newL = {shape_envL[0], shape_envL[1], S.size()};
        size_t newL_size = shape_newL[0] * shape_newL[1]* shape_newL[2];
        envL_internal.resize(newL_size);

        // Define maps to work with the pointers
        Eigen::TensorMap<Eigen::Tensor<const T, 3>> envL_map(envL, shape_envL);
        Eigen::TensorMap<Eigen::Tensor<T,3>> envL_internal_map (envL_internal.data(), shape_newL);
        envL_internal_map  = envL_map.contract(U_r2l, Textra::idx({2},{0}));

        // Reseat the pointer and update the size
        envL = envL_internal.data();
        shape_envL = shape_newL;
    }
    {

        // Copy the new mpo
        shape_mpo = mpo_r2l.dimensions();
        mpo_internal.resize(mpo_r2l.size());
        Eigen::TensorMap<Eigen::Tensor<T,4>> mpo_internal_map (mpo_internal.data(), shape_mpo);
        mpo_internal_map = mpo_r2l;

        // Reseat the pointer
        mpo = mpo_internal.data();
    }
    readyCompress = true;
    eig::log->trace("Compressed mpo² dimensions {}", shape_mpo);
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

// Explicit instantiations
template class MatrixProductHamiltonian<double>;
template class MatrixProductHamiltonian<std::complex<double>>;
