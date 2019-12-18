//
// Created by david on 2019-12-17.
//
#pragma once
//#define EIGEN_MALLOC_ALREADY_ALIGNED 0

#include <complex.h>
#undef I
#include <ceres/ceres.h>
#include <glog/logging.h>
//#include <general/nmspc_omp.h>

namespace opt{
    void SolveRosenbrock();
    inline ceres::GradientProblemSolver::Options ceres_default_options;

    class RosenbrockBase : public ceres::FirstOrderFunction {
    public:
    protected:
        Eigen::MatrixXd H;
        Eigen::MatrixXd H2;
        mutable double variance;
        mutable double energy  ;
//        OMP omp;
    public:
//        explicit RosenbrockBase(Eigen::MatrixXd & H_);
        explicit RosenbrockBase(const Eigen::MatrixXd & H_, const Eigen::MatrixXd & H2_);
        int NumParameters() const final;
        double get_variance   () const;

    };
}


namespace opt {
    template<typename T>
    class Rosenbrock : public RosenbrockBase {
    private:
        using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
        using VectorType = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    public:
//        explicit Rosenbrock(MatrixType &H_);
        explicit Rosenbrock(const MatrixType &H_, const MatrixType &H2_);
        bool Evaluate(const double *v_ptr,
                      double *fx,
                      double *grad_ptr) const final;

    };
}