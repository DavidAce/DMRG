//
// Created by david on 2018-10-07.
//

#ifndef DMRG_CLASS_LM_OPTIMIZATION_H
#define DMRG_CLASS_LM_OPTIMIZATION_H
#include <Eigen/Eigen>
#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>
#include <iostream>

struct lm_optimization {
    template <typename T> int sgn(const T val) const {
        return (T(0) < val) - (val < T(0));
    }
    using cScalar = std::complex<double>;
    using cMatrixType_ = Eigen::Matrix<cScalar,Eigen::Dynamic, Eigen::Dynamic>;
    using cVectorType_ = Eigen::Matrix<cScalar,Eigen::Dynamic, 1>;


    int operator()(const Eigen::Matrix<double,Eigen::Dynamic,1> &x, Eigen::Matrix<double,Eigen::Dynamic,1> &fvec) const
    {
        // 'x' has dimensions n x 1
        // It contains the current estimates for the parameters.

        // 'fvec' has dimensions m x 1
        // It will contain the fitness.
        auto eigvals  = Eigen::Map<const cVectorType_> (eigvals_ptr   ,nev);
        auto   v      = Eigen::Map<const Eigen::VectorXd>(x.data(), x.size());
        double lambda = 1.0;//x[v.size()];
        double vH2v   = (v.transpose() * H2 * v).real().sum(); // sum  x_i x_j H2
        double vEv    = v.cwiseAbs2().cwiseProduct(eigvals.real()).sum(); // sum  x_i x_i E
        double vv     = v.cwiseAbs2().sum();
        double var    = vH2v/vv - vEv*vEv/vv/vv;
        double logvar = std::log(var);
        fvec.setZero();
        fvec(0) = logvar  + lambda * std::abs(vv-1.0);

        return 0;
    }


    // Compute the jacobian
    int df(const Eigen::Matrix<double,Eigen::Dynamic,1> &x, Eigen::MatrixXd &fjac) const {



        // 'fjac' has dimensions m x n
        // It will contain the jacobian
        auto eigvals  = Eigen::Map<const cVectorType_> (eigvals_ptr   ,nev);
        auto   v      = Eigen::Map<const Eigen::VectorXd>(x.data(), x.size());
        double lambda = 1;
        double vH2v   = (v.transpose() * H2 * v).real().sum();                // sum  x_i x_j H2
        double vEv    = v.cwiseAbs2().cwiseProduct(eigvals.real()).sum();     // sum  x_i x_i E
        double vv     = v.cwiseAbs2().sum();
        double var    = vH2v/vv - vEv*vEv/vv/vv;
        fjac.setZero();
        for (int k = 0; k < v.size(); k++){
            double vi2H2ik         = 2.0*v.dot(H2.col(k)).real();             // 2 * sum_i x_i H2_ik
            double vk4EkvEv        = 4.0*v(k)*eigvals(k).real() * vEv ;       // 4 * x_k * E_k * (sum_i x_i^2 E_i)
            fjac(0,k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vEv*vEv*4.0*v(k)*vv)/(std::pow(vv,4)))/var
                      + lambda * 2.0 * v(k) * sgn(vv - 1);
        }
        return 0;
    }


    // Number of data points, i.e. values.
    int m;

    // Returns 'm', the number of values.
    int values() const { return m; }

    // The number of parameters, i.e. inputs.
    int n;

    // Returns 'n', the number of inputs.
    int inputs() const { return n; }


    const size_t shape;
    const size_t nev;
    const cScalar *H_local_ptr;
    const cScalar *H_local_sq_ptr;
    const cScalar *eigvecs_ptr;
    const cScalar *eigvals_ptr;

    cMatrixType_ H;
    cMatrixType_ H2;

    lm_optimization(
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
        H  = eigvecs.adjoint() * H_local    * eigvecs;
        H2 = eigvecs.adjoint() * H_local_sq * eigvecs;
        n = (int)H.rows();
        m = n;

    };



};


#endif //DMRG_CLASS_LM_OPTIMIZATION_H
