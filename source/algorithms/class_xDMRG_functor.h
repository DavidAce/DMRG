//
// Created by david on 2018-10-19.
//

#ifndef DMRG_CLASS_XDMRG_FUNCTOR_H
#define DMRG_CLASS_XDMRG_FUNCTOR_H
#include <Eigen/Core>
#include <Eigen/Dense>
#include <complex>

class class_xDMRG_functor {
public:
    template <typename T>
    int sgn(const T val) const {
        return (T(0) < val) - (val < T(0));
    }
    using cScalar = std::complex<double>;
    using cMatrixType_ = Eigen::Matrix<cScalar,Eigen::Dynamic, Eigen::Dynamic>;
    using cVectorType_ = Eigen::Matrix<cScalar,Eigen::Dynamic, 1>;

    const size_t shape;
    const size_t nev;
    const cScalar *H_local_ptr;
    const cScalar *H_local_sq_ptr;
    const cScalar *eigvecs_ptr;
    const cScalar *eigvals_ptr;

    cMatrixType_ H;
    cMatrixType_ H2;

    class_xDMRG_functor(
            const size_t shape_,
            const size_t nev_,
            const cScalar *H_local_ptr_,
            const cScalar *H_local_sq_ptr_,
            const cScalar *eigvecs_ptr_,
            const cScalar *eigvals_ptr_);

    double operator()(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, Eigen::Matrix<double,Eigen::Dynamic,1> &grad) const;

};


#endif //DMRG_CLASS_XDMRG_FUNCTOR_H
