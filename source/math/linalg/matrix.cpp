
#include "matrix.h"
#include <unsupported/Eigen/KroneckerProduct>

template<typename DerivedA,typename DerivedB>
linalg::matrix::KroneckerResultType<DerivedA,DerivedB>
linalg::matrix::kronecker(const Eigen::PlainObjectBase<DerivedA> & A, const Eigen::PlainObjectBase<DerivedB> & B, bool mirror){
    if (mirror) return kroneckerProduct(B.derived(),A.derived());
    else return kroneckerProduct(A.derived(),B.derived());
}

template linalg::matrix::KroneckerResultType<Eigen::MatrixXcd,Eigen::MatrixXcd>
linalg::matrix::kronecker(const Eigen::PlainObjectBase<Eigen::MatrixXcd> & A, const Eigen::PlainObjectBase<Eigen::MatrixXcd> & B, bool mirror);

template linalg::matrix::KroneckerResultType<Eigen::MatrixXd,Eigen::MatrixXd>
linalg::matrix::kronecker(const Eigen::PlainObjectBase<Eigen::MatrixXd> & A, const Eigen::PlainObjectBase<Eigen::MatrixXd> & B, bool mirror);

template linalg::matrix::KroneckerResultType<Eigen::Matrix2cd,Eigen::Matrix2cd>
linalg::matrix::kronecker(const Eigen::PlainObjectBase<Eigen::Matrix2cd> & A, const Eigen::PlainObjectBase<Eigen::Matrix2cd> & B, bool mirror);
