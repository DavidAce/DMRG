//
// Created by david on 2019-05-31.
//
#include <mps_routines/mps_tools/finite/opt.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <model/class_hamiltonian_base.h>
#include <algorithms/class_simulation_state.h>

MPS_Tools::Finite::Opt::internals::subspace_functor::subspace_functor(
        const class_superblock &superblock,
        const class_simulation_state &sim_state,
        const Eigen::MatrixXd & eigvecs_,
        const Eigen::VectorXd & eigvals_)
        :
        base_functor(superblock,sim_state),
        eigvecs(eigvecs_),
        eigvals(eigvals_)
{
//    eigvecs = eigvecs_;
//    eigvals = eigvals_;
    H2 = (eigvecs.adjoint() * superblock.get_H_local_sq_matrix_real().selfadjointView<Eigen::Upper>() * eigvecs);
}






double MPS_Tools::Finite::Opt::internals::subspace_functor::operator()(const Eigen::VectorXd &v, Eigen::VectorXd &grad) {
    long double vH2v,vEv,vv,var;
    double lambda,lambda2,log10var, fx, energy_penalty;

#pragma omp parallel
    {
#pragma omp sections
        {
#pragma omp section
            { vH2v = v.dot(H2.selfadjointView<Eigen::Upper>() * v); }
#pragma omp section
            { vEv = v.cwiseAbs2().dot(eigvals); }
#pragma omp section
            { vv = v.squaredNorm(); }
        }
#pragma omp barrier
#pragma omp single
        {
            lambda          = 1.0;
            lambda2         = 1.0;
            var             = std::abs(vH2v * vv - vEv * vEv) / std::pow(vv,2);
            variance        = var;
            var             = var == 0  ? std::numeric_limits<double>::epsilon() : var;
            energy          = vEv / vv;
            double energy_distance = std::abs(energy - energy_target);
            energy_penalty = energy_distance > energy_window/2.0 ?  energy - energy_target : 0.0;
            log10var = std::log10(var);
        }


#pragma omp barrier
#pragma omp for schedule(static,1)
        for (int k = 0; k < grad.size(); k++){
            double vi2H2ik         = 2.0*v.dot(H2.col(k));             // 2 * sum_i x_i H2_ik
            double vk4EkvEv        = 4.0*v(k)*eigvals(k) * vEv ;       // 4 * x_k * E_k * (sum_i x_i^2 E_i)
            grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/(std::pow(vv,2)) - (vk4EkvEv*vv*vv - vEv*vEv*4.0*v(k)*vv)/(std::pow(vv,4)))/var/std::log(10)
                      + lambda * 4.0 * v(k) * (vv - 1);
            if (have_bounds_on_energy){
                grad(k) += lambda2* 4 * v(k) * energy_penalty*(eigvals(k)/vv - energy/vv);
            }
            grad(k) = std::isinf(grad(k)) ? std::numeric_limits<double>::max() * sgn(grad(k)) : grad(k);

        }

    }

    if(std::isnan(log10var) or std::isinf(log10var)){
        MPS_Tools::log->warn("log10 variance is invalid");
        MPS_Tools::log->warn("vH2v            = {}" , vH2v );
        MPS_Tools::log->warn("vEv             = {}" , vEv  );
        MPS_Tools::log->warn("vv              = {}" , vv   );
        MPS_Tools::log->warn("vH2v/vv         = {}" , vH2v/vv    );
        MPS_Tools::log->warn("vEv*vEv/vv/vv   = {}" , vEv*vEv/vv/vv    );
        MPS_Tools::log->warn("var             = {}" , var);
        MPS_Tools::log->warn("log10var        = {}" , log10var );
//                exit(1);
        log10var    = std::abs(var) ==0  ?  -20.0 : std::log10(std::abs(var));
    }
    fx = log10var  + lambda * std::pow(vv-1.0,2);

    if (have_bounds_on_energy) {
        fx += lambda2 * std::pow(energy_penalty,2);
    }


//    grad.normalize();
//    Eigen::VectorXd m_drt;
//    m_drt.resizeLike(grad);
//    m_drt.noalias() = -grad;
//    double step  = 1.0 / grad.norm();
//    double step2 = 1.0 / m_drt.norm();
//    double norm  = grad.norm()  ;
//    double norm2 = m_drt.norm() ;
//
//    if(step < 1e-15 or step2 < 1e-15){
//        throw std::runtime_error("step too small");
//    }
//    if(step > 1e15 or step2 > 1e15){
//        throw std::runtime_error("step too large");
//    }

    counter++;
    return fx;
}
