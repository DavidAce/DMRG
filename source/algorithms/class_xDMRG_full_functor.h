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
public:
    template <typename T>
    int sgn(const T val) const {
        return (T(0) < val) - (val < T(0));
    }

    using MatrixType_ = Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic>;
    using VectorType_ = Eigen::Matrix<Scalar,Eigen::Dynamic, 1>;
    size_t   counter = 0;
//    const size_t shape;

    double get_variance(){return variance;}
    double get_energy  (){return energy  ;}
    size_t get_count   (){return counter;}
    Eigen::Tensor<double,4> HA_MPO;
    Eigen::Tensor<double,4> HB_MPO;
    Eigen::Tensor<double,3> Lblock;
    Eigen::Tensor<double,3> Rblock;
    Eigen::Tensor<double,4> Lblock2;
    Eigen::Tensor<double,4> Rblock2;
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
    }





    double get_vH2v(const Eigen::Matrix<double,Eigen::Dynamic,1> &v){
//        t_lbfgs.tic();

        auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), dsizes);
        Eigen::Tensor<double, 0> H2 =
                        Lblock2
                        .contract(theta,                Textra::idx({0}  ,{1}))
                        .contract(HA_MPO,               Textra::idx({1,3},{0,2}))
                        .contract(HB_MPO,               Textra::idx({4,2},{0,2}))
                        .contract(HA_MPO,               Textra::idx({1,3},{0,2}))
                        .contract(HB_MPO,               Textra::idx({4,3},{0,2}))
                        .contract(theta.conjugate(),    Textra::idx({0,3,5},{1,0,2}))
                        .contract(Rblock2,              Textra::idx({0,3,1,2},{0,1,2,3})).real();
//        t_lbfgs.toc();
        return H2(0);
    }

    double get_vHv(const Eigen::Matrix<double,Eigen::Dynamic,1> &v){
        t_lbfgs.tic();

        auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), dsizes);
        Eigen::Tensor<double, 0>  E =
                         Lblock
                        .contract(theta,                Textra::idx({0},{1}))
                        .contract(HA_MPO,               Textra::idx({1,2},{0,2}))
                        .contract(HB_MPO,               Textra::idx({3,1},{0,2}))
                        .contract(theta.conjugate(),    Textra::idx({0,2,4},{1,0,2}))
                        .contract(Rblock,               Textra::idx({0,2,1},{0,1,2}));

        t_lbfgs.toc();
        return E(0);
    }


    Eigen::VectorXd get_vH2 (const Eigen::Matrix<double,Eigen::Dynamic,1> &v){
//        t_lbfgs.tic();
        auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), dsizes);
        Eigen::Tensor<double, 4> vH2 =
                        Lblock2
                        .contract(theta,                 Textra::idx({0}  ,{1}))
                        .contract(HA_MPO,                Textra::idx({1,3},{0,2}))
                        .contract(HB_MPO,                Textra::idx({4,2},{0,2}))
                        .contract(HA_MPO,                Textra::idx({1,3},{0,2}))
                        .contract(HB_MPO,                Textra::idx({4,3},{0,2}))
                        .contract(Rblock2,               Textra::idx({1,2,4},{0,2,3}))
                        .shuffle(Textra::array4{1,0,2,3});
//        t_lbfgs.toc();
        return Eigen::Map<Eigen::VectorXd>(vH2.data(),vH2.size());
    }

    Eigen::VectorXd get_vH (const Eigen::Matrix<double,Eigen::Dynamic,1> &v){
//        t_lbfgs.tic();
        auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), dsizes);
        Eigen::Tensor<double, 4> vH = Lblock
                .contract(theta,                    Textra::idx({0},{1}))
                .contract(HA_MPO ,                  Textra::idx({1,2},{0,2}))//  idx({1,2,3},{0,4,5}))
                .contract(HB_MPO ,                  Textra::idx({3,1},{0,2}))//  idx({1,2,3},{0,4,5}))
                .contract(Rblock ,                  Textra::idx({1,3},{0,2}))
                .shuffle(Textra::array4{1,0,2,3});
//        t_lbfgs.toc();
        return Eigen::Map<Eigen::VectorXd>(vH.data(),vH.size());
    }


    //    double operator()(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, Eigen::Matrix<double,Eigen::Dynamic,1> &grad) const;
    double operator()(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, Eigen::Matrix<double,Eigen::Dynamic,1> &grad) {
        double vH2v,vHv,vv;
        double lambda,var,log10var, fx;
        Eigen::VectorXd vH, vH2;
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
                lambda      = 1.0;
                var         = vH2v/vv - vHv*vHv/vv/vv;
                variance    = var;
                energy      = vHv/vv;
                log10var    = std::log10(var);
                // fx = log10var  + lambda * std::pow(vv-1.0,2);
                fx = log10var  +  lambda * std::abs(vv-1.0);
                if(std::isnan(log10var)){
                    std::cout << "v: \n" << v << std::endl;
                    std::cout << "vH2v            : " << vH2v << std::endl;
                    std::cout << "vHv             : " << vHv  << std::endl;
                    std::cout << "vv              : " << vv   << std::endl;
                    std::cout << "vH2v/vv         : " << vH2v/vv    << std::endl;
                    std::cout << "vEv*vEv/vv/vv   : " << vHv*vHv/vv/vv    << std::endl;
                    std::cout << "var             : " << var  << std::endl << std::endl;
                    exit(1);
                }
            }
            #pragma omp barrier
            #pragma omp for schedule(static,1)
            for (int k = 0; k < v.size(); k++){
                double vi2H2ik         = 2.0*vH2(k);             // 2 * sum_i x_i H2_ik
                double vk4EkvEv        = 4.0*vH(k) * vHv ;       // 4 * x_k * E_k * (sum_i x_i^2 E_i)
//            grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vHv*vHv*4.0*v(k)*vv)/(std::pow(vv,4)))/var/std::log(10)
//                      + lambda * 4.0 * v(k) * (vv - 1);
                grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vHv*vHv*4.0*v(k)*vv)/(std::pow(vv,4)))/var/std::log(10)
                          + lambda * 2.0 * v(k) * sgn(vv - 1);
            }
        }
        counter++;
        return fx;
    }


};


#endif //DMRG_CLASS_XDMRG_FULL_FUNCTOR_H
