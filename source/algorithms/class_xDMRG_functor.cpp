//
// Created by david on 2018-10-19.
//

#include "class_xDMRG_functor.h"
class_xDMRG_functor::class_xDMRG_functor(
        const size_t shape_,
        const size_t nev_,
        const cScalar *H_local_ptr_,
        const cScalar *H_local_sq_ptr_,
        const cScalar *eigvecs_ptr_,
        const cScalar *eigvals_ptr_)
        :
        shape(shape_),
        nev(nev_),
        H_local_ptr(H_local_ptr_),
        H_local_sq_ptr(H_local_sq_ptr_),
        eigvecs_ptr(eigvecs_ptr_),
        eigvals_ptr(eigvals_ptr_)
{
    auto H_local    = Eigen::Map<const cMatrixType_> (H_local_ptr   ,shape,shape);
    auto H_local_sq = Eigen::Map<const cMatrixType_> (H_local_sq_ptr,shape,shape);
    auto eigvecs    = Eigen::Map<const cMatrixType_> (eigvecs_ptr   ,shape,nev);
    H               = eigvecs.adjoint() * H_local * eigvecs;
    H2              = eigvecs.adjoint() * H_local_sq * eigvecs;
};


double class_xDMRG_functor::operator()(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, Eigen::Matrix<double,Eigen::Dynamic,1> &grad) const {
    auto eigvals  = Eigen::Map<const cVectorType_> (eigvals_ptr,nev);
    double lambda = 1.0;
    double vH2v   = (v.adjoint() * H2 * v).real().sum();
    double vEv    = v.cwiseAbs2().cwiseProduct(eigvals).real().sum();
    double vv     = v.cwiseAbs2().sum();
    double var    = vH2v/vv - vEv*vEv/vv/vv;
    double log10var = std::log10(var);
//        double fx = logvar  + lambda * std::abs(vv-1.0);
    double fx = log10var  + lambda * std::pow(vv-1.0,2);
    for (int k = 0; k < v.size(); k++){
        double vi2H2ik         = 2.0*v.dot(H2.col(k)).real();             // 2 * sum_i x_i H2_ik
        double vk4EkvEv        = 4.0*v(k)*eigvals(k).real() * vEv ;       // 4 * x_k * E_k * (sum_i x_i^2 E_i)
//            grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vEv*vEv*4.0*v(k)*vv)/(std::pow(vv,4)))/var
//                      + lambda * 2.0 * v(k) * sgn(vv - 1);
        grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vEv*vEv*4.0*v(k)*vv)/(std::pow(vv,4)))/var/std::log(10)
                  + lambda * 4.0 * v(k) * (vv - 1);
    }
    return fx;
}