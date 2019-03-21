//
// Created by david on 2018-10-19.
//

#ifndef DMRG_CLASS_XDMRG_FUNCTOR_H
#define DMRG_CLASS_XDMRG_FUNCTOR_H
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <complex>
#include <spdlog/spdlog.h>
template<typename Scalar>
class class_xDMRG_functor {
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
    int   counter = 0;
//    void set_energy_bounds(double E_lower, double E_upper);
    bool have_bounds_on_energy = false;
    void set_energy_bounds(double E_lower, double E_upper) {
        energy_lower_bound = E_lower;
        energy_upper_bound = E_upper;
        energy_target = 0.5*(energy_upper_bound + energy_lower_bound);
        energy_window = 1.0*(energy_upper_bound - energy_lower_bound);
        have_bounds_on_energy = true;
        std::cout << "energy_lower_bound    = " << energy_lower_bound << std::endl;
        std::cout << "energy_upper_bound    = " << energy_upper_bound << std::endl;
        std::cout << "energy_target         = " << energy_target << std::endl;
        std::cout << "energy_window         = " << energy_window << std::endl;
    }
    const size_t shape;
    const size_t nev;

    const Scalar *eigvecs_ptr;
    const double *eigvals_ptr;

    Eigen::MatrixXd H2;
    double get_variance(){return variance;}
    double get_energy  (){return energy  ;}
    size_t get_count   (){return counter;}

    template<typename HType>
    class_xDMRG_functor(
            const size_t shape_,
            const size_t nev_,
            const HType  *H_local_sq_ptr_,
            const Scalar *eigvecs_ptr_,
            const double *eigvals_ptr_)
            :
            shape(shape_),
            nev(nev_),
            eigvecs_ptr(eigvecs_ptr_),
            eigvals_ptr(eigvals_ptr_)
    {
        auto H_local_sq = Eigen::Map<const Eigen::Matrix<HType,Eigen::Dynamic,Eigen::Dynamic>> (H_local_sq_ptr_,shape,shape);
        auto eigvecs    = Eigen::Map<const MatrixType_> (eigvecs_ptr   ,shape,nev);
        if constexpr(std::is_same<Scalar,HType>::value){
            H2 = (eigvecs.adjoint() * H_local_sq.derived().template selfadjointView<Eigen::Upper>() * eigvecs).real();
        }else{
            H2 = (eigvecs.adjoint() * H_local_sq.derived() * eigvecs).real();
        }

    }



    template <typename Derived>
    class_xDMRG_functor(
            const size_t shape_,
            const size_t nev_,
            const Eigen::EigenBase<Derived> &H_local_sq,
            const Scalar *eigvecs_ptr_,
            const double *eigvals_ptr_)
            :
            shape(shape_),
            nev(nev_),
            eigvecs_ptr(eigvecs_ptr_),
            eigvals_ptr(eigvals_ptr_)
    {
        auto eigvecs    = Eigen::Map<const MatrixType_> (eigvecs_ptr   ,shape,nev);
        if constexpr(std::is_same<Scalar,typename Derived::Scalar>::value){
            H2 = (eigvecs.adjoint() * H_local_sq.derived().template selfadjointView<Eigen::Upper>() * eigvecs).real();
        }else{
            H2 = (eigvecs.adjoint() * H_local_sq.derived() * eigvecs).real();
        }
    }


//    double operator()(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, Eigen::Matrix<double,Eigen::Dynamic,1> &grad) const;
    double operator()(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, Eigen::Matrix<double,Eigen::Dynamic,1> &grad) {
        auto eigvals  = Eigen::Map<const VectorType_> (eigvals_ptr,nev);
        double vH2v,vEv,vv,energy_penalty;
        double lambda,lambda2,var,log10var, fx;
        #pragma omp parallel
        {
            #pragma omp sections
            {
                #pragma omp section
                { vH2v = (v.adjoint() * H2.selfadjointView<Eigen::Upper>() * v).real().sum(); }
                #pragma omp section
                { vEv = v.cwiseAbs2().cwiseProduct(eigvals).real().sum(); }
                #pragma omp section
                { vv = v.cwiseAbs2().sum(); }
            }
            #pragma omp barrier
            #pragma omp single
            {
                lambda = 1.0;
                lambda2 = 1.0;
                var = vH2v / vv - vEv * vEv / vv / vv;
                variance = var;
                energy = vEv / vv;
//                energy_penalty = std::abs(energy - energy_target) - 0.5*energy_window;
                double energy_distance = std::abs(energy - energy_target);
                energy_penalty = energy_distance > energy_window/2.0 ?  energy - energy_target : 0.0;
                log10var = std::log10(var);
//                double fx = log10var  + lambda * std::pow(vv-1.0,2);
                if(std::isnan(log10var)){
                    spdlog::warn("log10 variance is NAN");
//                std::cout << "v: \n" << v << std::endl;
                    spdlog::warn("vH2v            =  " , vH2v );
                    spdlog::warn("vEv             =  " , vEv  );
                    spdlog::warn("vv              =  " , vv   );
                    spdlog::warn("vH2v/vv         =  " , vH2v/vv    );
                    spdlog::warn("vEv*vEv/vv/vv   =  " , vEv*vEv/vv/vv    );
                    spdlog::warn("var             =  " , var);
//                exit(1);

                    log10var    = std::abs(var) ==0  ?  -20.0 : std::log10(std::abs(var));
                }
                fx = log10var  +  lambda * std::abs(vv-1.0);
                if (have_bounds_on_energy) {
                    fx += lambda2 * std::pow(energy_penalty,2);
                }
//                if (energy_penalty > 0 and have_bounds_on_energy){
//                }


            }
            #pragma omp barrier
            #pragma omp for schedule(static,1)
            for (int k = 0; k < v.size(); k++){
                double vi2H2ik         = 2.0*v.dot(H2.col(k));             // 2 * sum_i x_i H2_ik
                double vk4EkvEv        = 4.0*v(k)*eigvals(k) * vEv ;       // 4 * x_k * E_k * (sum_i x_i^2 E_i)
//            grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vEv*vEv*4.0*v(k)*vv)/(std::pow(vv,4)))/var/std::log(10)
//                      + lambda * 4.0 * v(k) * (vv - 1);
                grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vEv*vEv*4.0*v(k)*vv)/(std::pow(vv,4)))/var/std::log(10)
                          + lambda * 2.0 * v(k) * sgn(vv - 1);
                if (have_bounds_on_energy){
//                    grad(k) -= lambda2*sgn(energy - energy_target) * 2*v(k)/vv*(eigvals(k) - energy) / energy_penalty;
                    grad(k) += lambda2* 4 * v(k) * energy_penalty*(eigvals(k)/vv - energy/vv);
                }

            }
        }
        counter++;
        return fx;
    }


};


#endif //DMRG_CLASS_XDMRG_FUNCTOR_H
