//
// Created by david on 2018-11-30.
//

#include "class_xDMRG_full_functor.h"
#include <spdlog/spdlog.h>

template<typename Scalar>
class_xDMRG_full_functor<Scalar>::class_xDMRG_full_functor(
        const Eigen::Tensor<Scalar,4> &HA_MPO_,
        const Eigen::Tensor<Scalar,4> &HB_MPO_,
        const Eigen::Tensor<Scalar,3> &Lblock_,
        const Eigen::Tensor<Scalar,3> &Rblock_,
        const Eigen::Tensor<Scalar,4> &Lblock2_,
        const Eigen::Tensor<Scalar,4> &Rblock2_,
        const Eigen::DSizes<long,4>   &dsizes_
)

//            HA_MPO(HA_MPO_),
//            HB_MPO(HB_MPO_),
//            Lblock(Lblock_),
//            Rblock(Rblock_),
//            Lblock2(Lblock2_),
//            Rblock2(Rblock2_),
//            dsizes(dsizes_)
{
    t_lbfgs.set_properties(true,5,"");
    HA_MPO = HA_MPO_.real();
    HB_MPO = HB_MPO_.real();
    Lblock = Lblock_.real();
    Rblock = Rblock_.real();
    Lblock2 = Lblock2_.real();
    Rblock2 = Rblock2_.real();
    dsizes  = dsizes_;
    HAHB    = HA_MPO.contract(HB_MPO, Textra::idx({1},{0}));
    HAHB2   = HAHB.contract(HAHB, Textra::idx({2,5},{1,4}));
}

template<typename Scalar>
void class_xDMRG_full_functor<Scalar>::set_energy_bounds(double E_lower, double E_upper) {
    energy_lower_bound = E_lower;
    energy_upper_bound = E_upper;
    energy_target = 0.5*(energy_upper_bound + energy_lower_bound);
    energy_window = 1.0*(energy_upper_bound - energy_lower_bound);


    have_bounds_on_energy = true;
}


template<typename Scalar>
double class_xDMRG_full_functor<Scalar>::get_vH2v(const Eigen::Matrix<double,Eigen::Dynamic,1> &v){
//        t_opt.tic();

    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), dsizes);
    Eigen::Tensor<double, 0> H2 =
            Lblock2
                    .contract(theta,                Textra::idx({0}  ,{1}))
                    .contract(HAHB2,                Textra::idx({2,1,3,4},{4,0,1,3}))
                    .contract(theta.conjugate(),    Textra::idx({0,3,5},{1,0,2}))
                    .contract(Rblock2,              Textra::idx({0,3,1,2},{0,1,2,3})).real();
//        t_opt.toc();
    return H2(0);
}

template<typename Scalar>
double class_xDMRG_full_functor<Scalar>::get_vHv(const Eigen::Matrix<double,Eigen::Dynamic,1> &v){
    t_lbfgs.tic();

    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), dsizes);
    Eigen::Tensor<double, 0>  E =
            Lblock
                    .contract(theta,                Textra::idx({0},{1}))
                    .contract(HAHB,                Textra::idx({1,2,3},{0,1,4}))
                    .contract(theta.conjugate(),    Textra::idx({0,2,4},{1,0,2}))
                    .contract(Rblock,               Textra::idx({0,2,1},{0,1,2}));

    t_lbfgs.toc();
    return E(0);
}


template<typename Scalar>
Eigen::VectorXd class_xDMRG_full_functor<Scalar>::get_vH2 (const Eigen::Matrix<double,Eigen::Dynamic,1> &v){
//        t_opt.tic();
    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), dsizes);
    Eigen::Tensor<double, 4> vH2 =
            Lblock2
                    .contract(theta,                 Textra::idx({0}  ,{1}))
                    .contract(HAHB2,                Textra::idx({2,1,3,4},{4,0,1,3}))
                    .contract(Rblock2,               Textra::idx({1,2,4},{0,2,3}))
                    .shuffle(Textra::array4{1,0,2,3});
//        t_opt.toc();
    return Eigen::Map<Eigen::VectorXd>(vH2.data(),vH2.size());
}

template<typename Scalar>
Eigen::VectorXd class_xDMRG_full_functor<Scalar>::get_vH (const Eigen::Matrix<double,Eigen::Dynamic,1> &v){
//        t_opt.tic();
    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), dsizes);
    Eigen::Tensor<double, 4> vH = Lblock
            .contract(theta,                    Textra::idx({0},{1}))
            .contract(HAHB,                    Textra::idx({1,2,3},{0,1,4}))
            .contract(Rblock ,                  Textra::idx({1,3},{0,2}))
            .shuffle(Textra::array4{1,0,2,3});
//        t_opt.toc();
    return Eigen::Map<Eigen::VectorXd>(vH.data(),vH.size());
}



template<typename Scalar>
double class_xDMRG_full_functor<Scalar>::operator()(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, Eigen::Matrix<double,Eigen::Dynamic,1> &grad) {
    double vH2v,vHv,vv,energy_penalty;
    double lambda,var,log10var, fx;
    Eigen::VectorXd vH, vH2;
//    num_threads(2)
    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {vH2v   = get_vH2v(v);}
            #pragma omp section
            {vHv    = get_vHv(v);}
            #pragma omp section
            {vv     = v.cwiseAbs2().sum();}
            #pragma omp section
            {vH2    = get_vH2(v);}
            #pragma omp section
            {vH     = get_vH(v);}
        }
        #pragma omp barrier
        #pragma omp single
        {
            lambda         = 1.0;
            var            = vH2v/vv - vHv*vHv/vv/vv;
            variance       = var;
            energy         = vHv/vv;
            energy_penalty = std::abs(energy - energy_target) - 0.5*energy_window;
            log10var    = std::log10(var);
//            log10norm   = std::log10(std::abs(vv));
//            fx = log10var  +  lambda * std::abs(vv-1.0);
//            fx = log10var  +  lambda * log10norm;
            if(std::isnan(log10var)){
                spdlog::warn("log10 variance is NAN");
//                std::cout << "v: \n" << v << std::endl;
                spdlog::warn("vH2v            =  " , vH2v );
                spdlog::warn("vHv             =  " , vHv  );
                spdlog::warn("vv              =  " , vv   );
                spdlog::warn("vH2v/vv         =  " , vH2v/vv    );
                spdlog::warn("vEv*vEv/vv/vv   =  " , vHv*vHv/vv/vv    );
                spdlog::warn("var             =  " , var);
//                exit(1);

                log10var    = std::abs(var) ==0  ?  -20.0 : std::log10(std::abs(var));
            }
            fx = log10var  + lambda * std::pow(vv-1.0,2);
            if (energy_penalty > 0 and have_bounds_on_energy){
                fx += energy_penalty;
            }
        }
        #pragma omp barrier
        #pragma omp for schedule(static,1)
        for (int k = 0; k < v.size(); k++){
            double vi2H2ik         = 2.0*vH2(k);             // 2 * sum_i x_i H2_ik
            double vk4EkvEv        = 4.0*vH(k) * vHv ;       // 4 * x_k * E_k * (sum_i x_i^2 E_i)
            grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vHv*vHv*4.0*v(k)*vv)/(std::pow(vv,4)))/var/std::log(10)
                      + lambda * 4.0 * v(k) * (vv - 1);
            if (energy_penalty > 0 and have_bounds_on_energy){
//                grad(k) +=  (2*vH(k)/vv - 2*v(k)*vHv/vv/vv);
                grad(k) +=  sgn(energy-energy_target) *2/vv*(vH(k) - v(k)*energy);
            }

//            grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vHv*vHv*4.0*v(k)*vv)/(std::pow(vv,4)))/var/std::log(10)
//                      + lambda * 2.0 * v(k) * sgn(vv - 1);//
//            grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vHv*vHv*4.0*v(k)*vv)/(std::pow(vv,4)))/var/std::log(10)
//                      + lambda * 2.0 * v(k) / vv / std::log(10);
        }
    }
    counter++;
    return fx;
}

template class class_xDMRG_full_functor<double>;
template class class_xDMRG_full_functor<std::complex<double>>;