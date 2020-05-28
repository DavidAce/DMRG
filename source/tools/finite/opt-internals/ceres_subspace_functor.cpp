//
// Created by david on 2019-07-15.
//
#include "ceres_subspace_functor.h"
#include <tensors/state/class_state_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/class_tensors_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>

using namespace tools::finite::opt::internal;

template<typename Scalar>
tools::finite::opt::internal::ceres_subspace_functor<Scalar>::ceres_subspace_functor(
        const class_tensors_finite & tensors,
        const class_algorithm_status & status,
        const MatrixType & H2_subspace,
        const Eigen::VectorXd  & eigvals_)
        :
        ceres_base_functor(tensors,status),
        H2(H2_subspace),
        eigvals(eigvals_)
{

    energy_reduced = tensors.model->get_energy_reduced();
    num_parameters = static_cast<int>(eigvals.size());
    if constexpr (std::is_same<Scalar,std::complex<double>>::value){num_parameters *= 2;}
}



template<typename Scalar>
bool tools::finite::opt::internal::ceres_subspace_functor<Scalar>::Evaluate(const double* v_double_double,
                                                                             double* fx,
                                                                             double* grad_double_double) const
{
    using namespace tools::common::profile;
    t_op->tic();
    Scalar vH2v,vHv;
    Scalar ene,ene2,var;
    double vv, log10var;
    double norm_func,norm_grad;
    VectorType Hv, H2v;
    int vecSize = NumParameters();
    if constexpr (std::is_same<Scalar,std::complex<double>>::value){vecSize = NumParameters()/2;}
    Eigen::Map<const VectorType> v (reinterpret_cast<const Scalar*>(v_double_double)   , vecSize);
    vv = v.squaredNorm();
    norm = std::sqrt(vv);

    t_vH->tic();
    Hv  = eigvals.asDiagonal() * v;
    t_vH->toc();

    t_vHv->tic();
    vHv = v.dot(Hv);
    t_vHv->toc();

    t_vH2->tic();
    H2v  = H2.template selfadjointView<Eigen::Upper>()*v;
    t_vH2->toc();

    t_vH2v->tic();
    vH2v = v.dot(H2v);
    t_vH2v->toc();


    // Do this next bit carefully to avoid negative variance when numbers are very small
    ene             = vHv/vv;
    ene2            = vH2v/vv;
    if (std::real(ene2) < 0.0 ) tools::log->debug("Counter = {}. ene2 is negative:  {:.16f} + i {:.16f}" , counter, std::real(ene2) , std::imag(ene2));
    ene2             = std::real(ene2) <  0.0 ? std::abs(ene2)                         : std::real(ene2);
    ene2             = std::real(ene2) == 0.0 ? std::numeric_limits<double>::epsilon() : std::real(ene2);

    var             = ene2 - ene*ene;
    if (std::real(var)  < 0.0 ) tools::log->debug("Counter = {}. var  is negative:  {:.16f} + i {:.16f}" , counter, std::real(var)  , std::imag(var));
    var             = std::real(var) <  0.0 ? std::abs(var)                          : std::real(var);
    var             = std::real(var) == 0.0 ? std::numeric_limits<double>::epsilon() : std::real(var);

    energy         = std::real(ene + energy_reduced) / static_cast<double>(length);
    variance       = std::abs(var)  / static_cast<double>(length);
    norm_offset    = std::abs(vv) - 1.0 ;
    std::tie(norm_func,norm_grad) = windowed_func_grad(norm_offset,0.05);
    log10var       = std::log10(variance);

    if (fx != nullptr){
        fx[0] = log10var +  norm_func;
    }

    Eigen::Map<VectorType>  grad (reinterpret_cast<Scalar*>(grad_double_double), vecSize);
    if (grad_double_double != nullptr){
        auto vv_1  = std::pow(vv,-1);
        auto var_1 = 1.0/var/std::log(10);
        grad = var_1 * vv_1 * (H2v - 2.0*ene*Hv - (ene2 - 2.0*ene*ene)*v);
        if constexpr (std::is_same<Scalar,double>::value){
            grad *= 2.0;
        }
        grad += norm_grad * v;
    }



//    tools::log->trace("log10 var: {:<24.18f} Energy: {:<24.18f} |Grad|: {:<24.18f} |Grad|_inf: {:<24.18f} SqNorm: {:<24.18f} Norm: {:<24.18f} Norm_func: {:<24.18f} |Norm_grad *v|: {:<24.18f} fx: {:<24.18f}",
//                      std::log10(std::abs(var)/length),
//                      std::real(ene + energy_reduced) / length,
//                      grad.norm(),
//                      grad.cwiseAbs().maxCoeff(),
//                      vv,
//                      norm,
//                      norm_func,
//                      (norm_grad * v).norm(),
//                      fx[0]);



    if(std::isnan(log10var) or std::isinf(log10var)){
        tools::log->warn("log₁₀ variance is invalid");
        tools::log->warn("vv              = {:.16f} + i{:.16f}" , std::real(vv)  , std::imag(vv));
        tools::log->warn("vH2v            = {:.16f} + i{:.16f}" , std::real(vH2v) ,std::imag(vH2v) );
        tools::log->warn("vHv             = {:.16f} + i{:.16f}" , std::real(vHv)  ,std::imag(vHv)  );
        tools::log->warn("var             = {:.16f} + i{:.16f}" , std::real(var)  ,std::imag(var));
        tools::log->warn("ene             = {:.16f} + i{:.16f}" , std::real(ene)  ,std::imag(ene));
        tools::log->warn("log₁₀(var/L)    = {:.16f}" , std::log10(variance/static_cast<double>(length)) );
        std::cout << "v: \n " << v << std::endl;
        throw std::runtime_error("Subspace functor failed at counter = " + std::to_string(counter) );
    }

    counter++;
    t_op->toc();
    return true;
}



template class tools::finite::opt::internal::ceres_subspace_functor<double>;
template class tools::finite::opt::internal::ceres_subspace_functor<std::complex<double>>;
