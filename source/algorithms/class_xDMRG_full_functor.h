//
// Created by david on 2018-11-30.
//

#ifndef DMRG_CLASS_XDMRG_FULL_FUNCTOR_H
#define DMRG_CLASS_XDMRG_FULL_FUNCTOR_H
#ifdef OpenMP_AVAILABLE
#include <omp.h>
#endif

#ifdef OpenBLAS_AVAILABLE
#include <cblas.h>
#endif

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include <general/nmspc_tensor_extra.h>
#include <general/class_tic_toc.h>
template<typename Scalar>
class class_xDMRG_full_functor {
private:
    double variance;
    double energy  ;
    double energy_lower_bound;
    double energy_upper_bound;
    double energy_target;
    double energy_window;
public:
    template <typename T>
    int sgn(const T val) const {
        return (T(0) < val) - (val < T(0));
    }

    using MatrixType_ = Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic>;
    using VectorType_ = Eigen::Matrix<Scalar,Eigen::Dynamic, 1>;
    size_t   counter = 0;
//    const size_t shape;
    void set_energy_bounds(double E_lower, double E_upper);
    bool have_bounds_on_energy = false;
    double get_variance(){return variance;}
    double get_energy  (){return energy  ;}
    size_t get_count   (){return counter;}
    Eigen::Tensor<double,4> HA_MPO;
    Eigen::Tensor<double,4> HB_MPO;
    Eigen::Tensor<double,3> Lblock;
    Eigen::Tensor<double,3> Rblock;
    Eigen::Tensor<double,4> Lblock2;
    Eigen::Tensor<double,4> Rblock2;

    Eigen::Tensor<double,6> HAHB;
    Eigen::Tensor<double,8> HAHB2;
    Eigen::DSizes<long,4>   dsizes;
    class_tic_toc t_lbfgs;
    class_xDMRG_full_functor(
            const Eigen::Tensor<Scalar,4> &HA_MPO_,
            const Eigen::Tensor<Scalar,4> &HB_MPO_,
            const Eigen::Tensor<Scalar,3> &Lblock_,
            const Eigen::Tensor<Scalar,3> &Rblock_,
            const Eigen::Tensor<Scalar,4> &Lblock2_,
            const Eigen::Tensor<Scalar,4> &Rblock2_,
            const Eigen::DSizes<long,4>   &dsizes_
            );

    double get_vH2v(const Eigen::Matrix<double,Eigen::Dynamic,1> &v);
    double get_vHv(const Eigen::Matrix<double,Eigen::Dynamic,1> &v);
    Eigen::VectorXd get_vH2 (const Eigen::Matrix<double,Eigen::Dynamic,1> &v);
    Eigen::VectorXd get_vH (const Eigen::Matrix<double,Eigen::Dynamic,1> &v);
    double operator()(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, Eigen::Matrix<double,Eigen::Dynamic,1> &grad);


};


#endif //DMRG_CLASS_XDMRG_FULL_FUNCTOR_H
