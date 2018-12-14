//
// Created by david on 2018-10-03.
//

#ifndef DMRG_CLASS_BFGS_OPTIMIZATION_H
#define DMRG_CLASS_BFGS_OPTIMIZATION_H
#include <Eigen/Core>
#include <cppoptlib/meta.h>
#include <cppoptlib/problem.h>
#include <cppoptlib/solver/bfgssolver.h>
#include <cppoptlib/solver/lbfgsbsolver.h>
#include <cppoptlib/solver/conjugatedgradientdescentsolver.h>
#include <cppoptlib/solver/newtondescentsolver.h>
#include <cppoptlib/solver/gradientdescentsolver.h>
#include <cppoptlib/solver/neldermeadsolver.h>
#include <cppoptlib/solver/cmaessolver.h>
#include <memory>

template<typename dScalar>
class class_bfgs_optimization  : public cppoptlib::Problem<dScalar> {
private:
    template <typename T> int sgn(T val) const {
        return (T(0) < val) - (val < T(0));
    }
public:
    using cScalar = std::complex<double>;
    using MatrixType = Eigen::Matrix<cScalar,Eigen::Dynamic, Eigen::Dynamic>;
    using VectorType = Eigen::Matrix<cScalar,Eigen::Dynamic, 1>;
    using MatrixMap  = Eigen::Map<const MatrixType>;
    using VectorMap  = Eigen::Map<const VectorType>;

    using typename cppoptlib::Problem<dScalar>::TVector;
    using typename cppoptlib::Problem<dScalar>::THessian;
    // this is just the objective (NOT optional)
    dScalar value(const TVector &x) {
        auto   v      = Eigen::Map<const Eigen::VectorXd>(x.data(), x.size());
        double lambda = 1.0;//x[v.size()];
        double vH2v   = (v.transpose() * H2 * v).real().sum(); // sum  x_i x_j H2
        double vEv    = v.cwiseAbs2().cwiseProduct(eigvals.real()).sum(); // sum  x_i x_i E
        double vv     = v.cwiseAbs2().sum();
        double var    = vH2v/vv - vEv*vEv/vv/vv;
        double logvar = std::log(var);
//        return logvar  + lambda * std::abs(vv-1.0);
        return logvar  + lambda * std::pow(vv-1.0,2);
    }

    void gradient(const TVector &x, TVector &grad) {
        auto   v      = Eigen::Map<const Eigen::VectorXd>(x.data(), x.size());
        double lambda = 1;
        double vH2v   = (v.transpose() * H2 * v).real().sum();                // sum  x_i x_j H2
        double vEv    = v.cwiseAbs2().cwiseProduct(eigvals.real()).sum();     // sum  x_i x_i E
        double vv     = v.cwiseAbs2().sum();
        double var    = vH2v/vv - vEv*vEv/vv/vv;
        for (int k = 0; k < v.size(); k++){
            double vi2H2ik         = 2.0*v.dot(H2.col(k)).real();             // 2 * sum_i x_i H2_ik
            double vk4EkvEv        = 4.0*v(k)*eigvals(k).real() * vEv ;       // 4 * x_k * E_k * (sum_i x_i^2 E_i)
//            grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vEv*vEv*4.0*v(k)*vv)/(std::pow(vv,4)))/var
//                        + lambda * 2.0 * v(k) * sgn(vv - 1);
            grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vEv*vEv*4.0*v(k)*vv)/(std::pow(vv,4)))/var
                      + lambda * 4.0 * v(k) * (vv - 1);
        }
    }

public:
    const size_t shape;
    const size_t nev;
    const MatrixType &H_local;
    const MatrixType &H_local_sq;
    const VectorType &thetavec;
    const MatrixType &eigvecs;
    const VectorType &eigvals;

//    MatrixMap H_map;

    MatrixType H;
    MatrixType H2;
    class_bfgs_optimization(
            const size_t shape_,
            const size_t nev_,
            const MatrixType &H_local_,
            const MatrixType &H_local_sq_,
            const VectorType &thetavec_,
            const MatrixType &eigvecs_,
            const VectorType &eigvals_)
    :
            shape(shape_),
            nev(nev_),
            H_local(H_local_),
            H_local_sq(H_local_sq_),
            thetavec(thetavec_),
            eigvecs(eigvecs_),
            eigvals(eigvals_)
    {
//        MatrixMap H_map (H_local_ptr.get(),shape,shape);
//        MatrixMap H2_map (H_local_sq.get(),shape,shape);
//        MatrixMap eigvecs_map (eigvecs.get(),shape,nev);

        H  = eigvecs.adjoint() * H_local    * eigvecs;
        H2 = eigvecs.adjoint() * H_local_sq * eigvecs;

    };





};


#endif //DMRG_CLASS_BFGS_OPTIMIZATION_H
