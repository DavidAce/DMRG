//
// Created by david on 2019-12-17.
//

#include "ceres_pthread.h"
//opt::RosenbrockBase::RosenbrockBase(Eigen::MatrixXd & H_):H(H_), H2(H_*H_){}
opt::RosenbrockBase::RosenbrockBase(const Eigen::MatrixXd & H_, const Eigen::MatrixXd & H2_):H(H_), H2(H2_){}

int opt::RosenbrockBase::NumParameters() const { return H.rows(); }
double  opt::RosenbrockBase::get_variance() const { return variance; }
