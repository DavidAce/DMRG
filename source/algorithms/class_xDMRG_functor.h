//
// Created by david on 2018-10-19.
//

#ifndef DMRG_CLASS_XDMRG_FUNCTOR_H
#define DMRG_CLASS_XDMRG_FUNCTOR_H
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <complex>

class class_xDMRG_functor {
public:
    template <typename T>
    int sgn(const T val) const {
        return (T(0) < val) - (val < T(0));
    }
    using cScalar = std::complex<double>;
    using  Scalar = double;
    using cMatrixType_ = Eigen::Matrix<cScalar,Eigen::Dynamic, Eigen::Dynamic>;
    using  MatrixType_ = Eigen::Matrix< Scalar,Eigen::Dynamic, Eigen::Dynamic>;
    using cVectorType_ = Eigen::Matrix<cScalar,Eigen::Dynamic, 1>;

    const size_t shape;
    const size_t nev;
//    const cScalar *H_local_ptr = nullptr;
//    const cScalar *H_local_sq_ptr = nullptr;
    const cScalar *eigvecs_ptr;
    const cScalar *eigvals_ptr;

//    cMatrixType_ H;
    MatrixType_ H2;

    class_xDMRG_functor(
            const size_t shape_,
            const size_t nev_,
//            const cScalar *H_local_ptr_,
            const cScalar *H_local_sq_ptr_,
            const cScalar *eigvecs_ptr_,
            const cScalar *eigvals_ptr_);

    template <typename Derived>
    class_xDMRG_functor(
            const size_t shape_,
            const size_t nev_,
//            const Eigen::EigenBase<Derived> &H_,
            const Eigen::EigenBase<Derived> &H_local_sq,
            const cScalar *eigvecs_ptr_,
            const cScalar *eigvals_ptr_)
            :
            shape(shape_),
            nev(nev_),
//            H(H_),
//            H2(H2_),
            eigvecs_ptr(eigvecs_ptr_),
            eigvals_ptr(eigvals_ptr_)
    {
        auto eigvecs    = Eigen::Map<const cMatrixType_> (eigvecs_ptr   ,shape,nev);
        H2 = (eigvecs.adjoint() * H_local_sq.derived().template selfadjointView<Eigen::Upper>() * eigvecs).real();

//        bool H2isHermitian;
//        H2isHermitian = H2.isApprox(H2.adjoint(), 1e-14);
//        std::cout << "H2 is hermitian: " << std::boolalpha << H2isHermitian << std::endl;


    }


    double operator()(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, Eigen::Matrix<double,Eigen::Dynamic,1> &grad) const;

};


#endif //DMRG_CLASS_XDMRG_FUNCTOR_H
