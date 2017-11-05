//
// Created by david on 2017-10-04.
//

#ifndef DMRG_CLASS_EIGENSOLVER_H
#define DMRG_CLASS_EIGENSOLVER_H


#include "n_tensor_extra.h"


class class_SVD{
private:
    double SVDThreshold         = 1e-12;
    double truncation_error     = 0;
    int chi                     = 0;
    template <long rankU, long rankS, long rankV>
    using TensorTuple = std::tuple<Textra::Tensor<rankU,double>, // U
                                   Textra::Tensor<rankS,double>, // S
                                   Textra::Tensor<rankV,double>>; // V^T
public:
    double get_truncation_error();
    void setThreshold(double newThreshold);
    TensorTuple<2,1,2> decompose(const Textra::Tensor2d &tensor);
    TensorTuple<2,1,2> decompose(const Textra::Tensor2d &tensor, const long chi_max);
    TensorTuple<3,1,3> schmidt(const Textra::Tensor2d &tensor, long d, long chiL, long chi_max, long chiR);
};







//enum class Handedness {RGHT, LEFT};
//enum class Solver {SYM, GEN};
//
//
//
//template <typename Scalar, Handedness hand = Handedness::RGHT, int Uplo = Eigen::Lower>
//class DenseSymMatProd
//{
//private:
//    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
//    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
//    typedef Eigen::Map<const Matrix> MapConstMat;
//    typedef Eigen::Map<const Vector> MapConstVec;
//    typedef Eigen::Map<Vector> MapVec;
//
//    typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;
//
//    const MapConstMat m_mat;
//
//public:
//    DenseSymMatProd(ConstGenericMatrix& mat_) :
//            m_mat(mat_.data(), mat_.rows(), mat_.cols())
//    {}
//
//    int rows() const { return m_mat.rows(); }
//    int cols() const { return m_mat.cols(); }
//
//// y_out = A * x_in
//    void perform_op(const Scalar* x_in, Scalar* y_out) const
//    {
//        MapConstVec x(x_in,  m_mat.cols());
//        MapVec      y(y_out, m_mat.rows());
//        y.noalias() = m_mat.template selfadjointView<Uplo>() * x;
//        switch (hand){
//            case Handedness::RGHT:
//                y.noalias() = m_mat.template selfadjointView<Uplo>() * x;
//                break;
//            case Handedness::LEFT:
//                y.noalias() = x.transpose() * m_mat.template selfadjointView<Uplo>();
//                break;
//        }
//
//    }
//};
//
//
//template <typename Scalar, Handedness hand = Handedness::RGHT>
//class DenseGenMatProd
//{
//private:
//    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
//    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
//    typedef Eigen::Map<const Matrix> MapConstMat;
//    typedef Eigen::Map<const Vector> MapConstVec;
//    typedef Eigen::Map<Vector> MapVec;
//    typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;
//    const MapConstMat m_mat;
//
//public:
//    DenseGenMatProd(ConstGenericMatrix& mat_) :
//            m_mat(mat_.data(), mat_.rows(), mat_.cols())
//    {}
//
//    int rows() const { return m_mat.rows(); }
//    int cols() const { return m_mat.cols(); }
//
//    // y_out = A * x_in
//    void perform_op(const Scalar* x_in, Scalar* y_out) const
//    {
//        MapConstVec x(x_in,  m_mat.cols());
//        MapVec      y(y_out, m_mat.rows());
//        switch (hand){
//            case Handedness::RGHT:
//                y.noalias() = m_mat * x;
//                break;
//            case Handedness::LEFT:
//                y.noalias() = x.transpose() * m_mat;
//                break;
//        }
//    }
//};

//
///*! \brief Returns the largest eigenvector of a matrix */
//template<Handedness hand = Handedness::RGHT>
//Textra::Tensor<1,Textra::cdouble> dominant_eigenvector(const Textra::Tensor<2, double> &tensor){
//    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
//    DenseGenMatProd<double, hand> op(mat);
//    int ncv = std::min(settings::precision::eig_max_ncv, 3);
//    Spectra::GenEigsSolver<double, Spectra::LARGEST_REAL, DenseGenMatProd<double,hand>> eigs(&op, 1, ncv);
//    eigs.init();
//    eigs.compute(50000, 1e-16, Spectra::LARGEST_REAL);
//    if(eigs.info() != Spectra::SUCCESSFUL){
//        std::cout << "Eigenvalue solver failed." << '\n';
////            std::cout << tensor << '\n';
////            exit(1);
//    }
//    std::cout << "Eigenvalue: " << eigs.eigenvalues()<< std::endl;
//    return Textra::Matrix_to_Tensor1<Textra::cdouble>(eigs.eigenvectors().real());
//}
//
//
//template<Solver solver = Solver::GEN, Handedness hand = Handedness::RGHT>
//class class_eigensolver {
//private:
//    DenseGenMatProd<double, hand> op;
//public:
//
//    class_eigensolver(){
//
//    }
//
//
//};


#endif //DMRG_CLASS_EIGENSOLVER_H
