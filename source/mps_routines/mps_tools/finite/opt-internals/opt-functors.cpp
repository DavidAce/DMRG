//
// Created by david on 2019-05-25.
//




#include <spdlog/spdlog.h>
#include <LBFGS.h>
#include <general/class_tic_toc.h>
#include <mps_routines/mps_tools/finite/opt.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <model/class_hamiltonian_base.h>
#include <algorithms/class_simulation_state.h>


LBFGSpp::LBFGSParam<double> MPS_Tools::Finite::Opt::internals::get_lbfgs_params(){
    using namespace LBFGSpp;
    LBFGSpp::LBFGSParam<double> param;
    // READ HERE http://pages.mtu.edu/~msgocken/ma5630spring2003/lectures/lines/lines/node3.html
    // I think c1 corresponds to ftol, and c2 corresponds to wolfe
    param.max_iterations = 2000;
    param.max_linesearch = 60; // Default is 20. 5 is really bad, 80 seems better.
    param.m              = 8;     // Default is 6
    param.past           = 1;     // Or perhaps it was this one that helped.
    param.epsilon        = 1e-5;  // Default is 1e-5.
    param.delta          = 1e-8;  // Default is 0. Trying this one instead of ftol.
    param.ftol           = 1e-4;  // Default is 1e-4. this really helped at threshold 1e-8. Perhaps it should be low. Ok..it didn't
    param.wolfe          = 0.9;   // Default is 0.9
    param.min_step       = 1e-40;
    param.max_step       = 1e+40;
    param.linesearch     = LINE_SEARCH_ALGORITHM::LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
    return param;
}



MPS_Tools::Finite::Opt::internals::direct_functor::direct_functor(
        const class_superblock & superblock_): base_functor()

{
    superblock.HA_MPO        = superblock_.HA->MPO().real();
    superblock.HB_MPO        = superblock_.HB->MPO().real();
    superblock.Lblock        = superblock_.Lblock->block.real();
    superblock.Rblock        = superblock_.Rblock->block.real();
    superblock.Lblock2       = superblock_.Lblock2->block.real();
    superblock.Rblock2       = superblock_.Rblock2->block.real();
    superblock.dsizes        = superblock_.dimensions();
    superblock.HAHB          = superblock.HA_MPO.contract(superblock.HB_MPO, Textra::idx({1},{0}));
    superblock.HAHB2         = superblock.HAHB.contract(superblock.HAHB, Textra::idx({2,5},{1,4}));
}


double MPS_Tools::Finite::Opt::internals::direct_functor::operator()(const Eigen::VectorXd &v, Eigen::VectorXd &grad) {
    t_op->tic();
    long double vH2v,vHv,vv,var;
    double lambda,log10var, fx;
    Eigen::VectorXd vH, vH2;

    #pragma omp parallel
    {
            #pragma omp sections
        {
        #pragma omp section
            {std::tie(vH2,vH2v)  = get_vH2_vH2v(v,superblock);}
        #pragma omp section
            {std::tie(vH,vHv)    = get_vH_vHv(v,superblock);}
        #pragma omp section
            {vv     = v.cwiseAbs2().sum();}
        }
    #pragma omp barrier

    #pragma omp single
        {
            lambda         = 1.0;
            var            = std::abs(vH2v*vv - vHv*vHv)/std::pow(vv,2);
            variance       = var;
            var            = var == 0  ? std::numeric_limits<double>::epsilon() : var;
            energy         = vHv/vv;
            log10var       = std::log10(var);

        }
    #pragma omp barrier
        auto vv2   = std::pow(vv,2);
        auto vv4   = std::pow(vv,4);
        auto varlog10 = 1.0/var/std::log(10.0);
    #pragma omp for schedule(static,1)
        for (int k = 0; k < v.size(); k++){
            double vi2H2ik         = 2.0*vH2(k);             // 2 * sum_i x_i H2_ik
            double vk4EkvEv        = 4.0*vH(k) * vHv ;       // 4 * x_k * E_k * (sum_i x_i^2 E_i)
            grad(k) = ((vi2H2ik * vv - vH2v * 2.0*v(k))/vv2 - (vk4EkvEv*vv*vv - vHv*vHv*4.0*v(k)*vv)/vv4)*varlog10
                      + lambda * 4.0 * v(k) * (vv - 1);
            grad(k) = std::isinf(grad(k)) ? std::numeric_limits<double>::max() * sgn(grad(k)) : grad(k);
        }
    }

    if(std::isnan(log10var) or std::isinf(log10var)){
        MPS_Tools::log->warn("log10 variance is invalid");
        std::cout << "v: \n" << v << std::endl;
        std::cout << "grad: \n" << grad << std::endl;
        MPS_Tools::log->warn("vH2v            = {}" , vH2v );
        MPS_Tools::log->warn("vHv             = {}" , vHv  );
        MPS_Tools::log->warn("vv              = {}" , vv   );
        MPS_Tools::log->warn("vH2v/vv         = {}" , vH2v/vv    );
        MPS_Tools::log->warn("vEv*vEv/vv/vv   = {}" , vHv*vHv/vv/vv    );
        MPS_Tools::log->warn("var             = {}" , var);
        exit(1);
//        log10var    = std::abs(var) == 0  ?  -20.0 : std::log10(std::abs(var));
    }
    fx = log10var  + lambda * std::pow(vv-1.0,2);

    counter++;
    t_op->toc();
    return fx;
}









MPS_Tools::Finite::Opt::internals::guided_functor::guided_functor(
        const class_superblock & superblock_, class_simulation_state &sim_state): base_functor()

{
    superblock.HA_MPO        = superblock_.HA->MPO().real();
    superblock.HB_MPO        = superblock_.HB->MPO().real();
    superblock.Lblock        = superblock_.Lblock->block.real();
    superblock.Rblock        = superblock_.Rblock->block.real();
    superblock.Lblock2       = superblock_.Lblock2->block.real();
    superblock.Rblock2       = superblock_.Rblock2->block.real();
    superblock.dsizes        = superblock_.dimensions();
    superblock.HAHB          = superblock.HA_MPO.contract(superblock.HB_MPO, Textra::idx({1},{0}));
    superblock.HAHB2         = superblock.HAHB.contract(superblock.HAHB, Textra::idx({2,5},{1,4}));
    energy_target            = sim_state.energy_target;
    energy_lower_bound       = sim_state.energy_lbound;
    energy_upper_bound       = sim_state.energy_ubound;


}


double MPS_Tools::Finite::Opt::internals::guided_functor::operator()(const Eigen::VectorXd &v_and_lambdas, Eigen::VectorXd &grad) {
    t_op->tic();
    long double vH2v,vHv,vv,var;
    double log10var, fx;
    Eigen::VectorXd vH, vH2;
    auto v       = Eigen::Map<const Eigen::VectorXd>(v_and_lambdas.data(), v_and_lambdas.size()-2);
//    auto lambdas = v_and_lambdas.tail(2);
    auto lambdas = Eigen::VectorXd::Ones(2);
//    lambdas.setConstant(1.0);
    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {std::tie(vH2,vH2v)  = get_vH2_vH2v(v,superblock);}
            #pragma omp section
            {std::tie(vH,vHv)    = get_vH_vHv(v,superblock);}
            #pragma omp section
            {vv     = v.dot(v);}
        }
        #pragma omp barrier

        #pragma omp single
        {

            energy         = vHv/vv;
            energy_offset  = energy - energy_target;

            var            = vH2v/vv - energy*energy;
            variance       = var;
            var            = var == 0  ? std::numeric_limits<double>::epsilon() : var;

            norm_offset    = vv - 1.0 ;
            log10var       = std::log10(var);
        }
        #pragma omp barrier
        auto vv_1  = std::pow(vv,-1);
        auto vv_2  = std::pow(vv,-2);
        auto var_1 = 1.0/var/std::log(10.0);
//        std::cout << "Variance: " << var << std::endl;
//        std::cout << "Energy  : " << energy << std::endl;
//        std::cout << "Energy o: " << energy_offset << std::endl;
//        std::cout << "norm o  : " << norm_offset << std::endl;
//        std::cout << "lambda 0: " << lambdas(0) << std::endl;
//        std::cout << "lambda 1: " << lambdas(1) << std::endl;
        #pragma omp for schedule(static,1)
        for (int k = 0; k < v.size(); k++){

            double vi2H2ikvv_1         = 2.0 * vH2(k)*vv_1;             // 2 * sum_i x_i H2_ik
            double vk2vH2vvv_2         = 2.0 * v(k) * vH2v * vv_2;

            double vi2Hikvv_1          = 2.0*vH(k)*vv_1;
            double vk2vHvvv_2          = 2.0*v(k) * vHv * vv_2;

            grad(k)  = vi2H2ikvv_1 - vk2vH2vvv_2;
            grad(k) -= 2.0*energy * (vi2Hikvv_1 -vk2vHvvv_2);
            grad(k) *= var_1;
//            grad(k) += 2.0 *lambdas(0) * sgn(norm_offset)   * v(k);
//            grad(k) += 1.0 *lambdas(1) * sgn(energy_offset) * (vi2Hikvv_1 - vk2vHvvv_2) ;
            grad(k) += 2.0 *std::abs(lambdas(0)) * norm_offset   * v(k);
            grad(k) += 2.0 *std::abs(lambdas(1)) * energy_offset * (vi2Hikvv_1 - vk2vHvvv_2) ;
//            grad(k) += 2.0 * sgn(norm_offset)   * v(k);
//            grad(k) += 2.0 * sgn(energy_offset) * (vi2Hikvv_1 - vk2vHvvv_2) ;
//            grad(k) = std::isinf(grad(k)) ? std::numeric_limits<double>::max() * sgn(grad(k)) : grad(k);
        }
//        grad(v.size())   = std::abs(norm_offset);
//        grad(v.size()+1) = std::abs(energy_offset);
        grad(v.size())   = sgn(lambdas(0)) * std::pow(norm_offset,2);
        grad(v.size()+1) = sgn(lambdas(1)) * std::pow(energy_offset,2);
//        grad(v.size())   = 0.0;
//        grad(v.size()+1) = 0.0;
    }

    if(std::isnan(log10var) or std::isinf(log10var)){
        MPS_Tools::log->warn("log10 variance is invalid");
        std::cout << "v: \n" << v << std::endl;
        std::cout << "grad: \n" << grad << std::endl;
        MPS_Tools::log->warn("vH2v            = {}" , vH2v );
        MPS_Tools::log->warn("vHv             = {}" , vHv  );
        MPS_Tools::log->warn("vv              = {}" , vv   );
        MPS_Tools::log->warn("vH2v/vv         = {}" , vH2v/vv    );
        MPS_Tools::log->warn("vEv*vEv/vv/vv   = {}" , vHv*vHv/vv/vv    );
        MPS_Tools::log->warn("var             = {}" , var);
        MPS_Tools::log->warn("lambda 0        = {}" , lambdas(0));
        MPS_Tools::log->warn("lambda 1        = {}" , lambdas(1));
        exit(1);
    }
//    fx = log10var  + lambdas(0) * std::abs(norm_offset) +  lambdas(1) * std::abs(energy_offset);
    fx = log10var  + std::abs(lambdas(0)) * std::pow(norm_offset,2) +  std::abs(lambdas(1)) * std::pow(energy_offset,2);

    counter++;
    t_op->toc();
    return fx;
}












std::pair<Eigen::VectorXd,double> MPS_Tools::Finite::Opt::internals::get_vH_vHv(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components & superblock){
    auto vH = MPS_Tools::Finite::Opt::internals::get_vH(v,superblock);
    t_vHv->tic();
    auto vHv = vH.cwiseProduct(v).sum();
    t_vHv->toc();
    return std::make_pair(vH,vHv);
}

std::pair<Eigen::VectorXd,double> MPS_Tools::Finite::Opt::internals::get_vH2_vH2v(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components & superblock){
    auto vH2 = MPS_Tools::Finite::Opt::internals::get_vH2(v,superblock);
    t_vH2v->tic();
    auto vH2v = vH2.cwiseProduct(v).sum();
    t_vH2v->toc();
    return std::make_pair(vH2,vH2v);
}




double MPS_Tools::Finite::Opt::internals::get_vH2v(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components & superblock){
    t_vH2v->tic();
    auto theta   = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), superblock.dsizes);

    Eigen::Tensor<double, 0> H2 =
            superblock.Lblock2
                    .contract(theta,                Textra::idx({0}  ,{1}))
                    .contract(superblock.HAHB2,     Textra::idx({2,1,3,4},{4,0,1,3}))
                    .contract(theta.conjugate(),    Textra::idx({0,3,5},{1,0,2}))
                    .contract(superblock.Rblock2,   Textra::idx({0,3,1,2},{0,1,2,3})).real();
    t_vH2v->toc();
    return H2(0);
}

double MPS_Tools::Finite::Opt::internals::get_vHv(const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components & superblock){
    t_vHv->tic();
    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), superblock.dsizes);
    Eigen::Tensor<double, 0>  E =
            superblock.Lblock
                    .contract(theta,                           Textra::idx({0},{1}))
                    .contract(superblock.HAHB,                 Textra::idx({1,2,3},{0,1,4}))
                    .contract(theta.conjugate(),               Textra::idx({0,2,4},{1,0,2}))
                    .contract(superblock.Rblock,               Textra::idx({0,2,1},{0,1,2}));
    t_vHv->toc();
    return E(0);
}


Eigen::VectorXd MPS_Tools::Finite::Opt::internals::get_vH2 (const Eigen::Matrix<double,Eigen::Dynamic,1> &v,const superblock_components & superblock){
    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), superblock.dsizes);
    t_vH2->tic();
    Eigen::Tensor<double, 4> vH2 =
            superblock.Lblock2
                    .contract(theta,                            Textra::idx({0},{1}))
                    .contract(superblock.HAHB2,                 Textra::idx({2,1,3,4},{4,0,1,3}))
                    .contract(superblock.Rblock2,               Textra::idx({1,2,4},{0,2,3}))
                    .shuffle(Textra::array4{1,0,2,3});
    t_vH2->toc();
    return Eigen::Map<Eigen::VectorXd>(vH2.data(),vH2.size());
}

Eigen::VectorXd MPS_Tools::Finite::Opt::internals::get_vH (const Eigen::Matrix<double,Eigen::Dynamic,1> &v, const superblock_components & superblock){
    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), superblock.dsizes);
    t_vH->tic();
    Eigen::Tensor<double, 4> vH =
            superblock.Lblock
                    .contract(theta,                               Textra::idx({0},{1}))
                    .contract(superblock.HAHB,                     Textra::idx({1,2,3},{0,1,4}))
                    .contract(superblock.Rblock ,                  Textra::idx({1,3},{0,2}))
                    .shuffle(Textra::array4{1,0,2,3});
    t_vH->toc();

    return Eigen::Map<Eigen::VectorXd>(vH.data(),vH.size());
}















// Subspace functors










MPS_Tools::Finite::Opt::internals::subspace_functor::subspace_functor(
        const class_superblock &superblock,
        const Eigen::MatrixXd & eigvecs_,
        const Eigen::VectorXd & eigvals_)
        :
        eigvecs(eigvecs_),
        eigvals(eigvals_)
{
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
            { vH2v = (v.adjoint() * H2.selfadjointView<Eigen::Upper>() * v).real().sum(); }
            #pragma omp section
            { vEv = v.cwiseAbs2().cwiseProduct(eigvals).real().sum(); }
            #pragma omp section
            { vv = v.cwiseAbs2().sum(); }
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


    //        grad.normalize();
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




