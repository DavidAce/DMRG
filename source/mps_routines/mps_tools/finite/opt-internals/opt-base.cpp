//
// Created by david on 2019-03-18.
//
#include <general/class_tic_toc.h>
#include <algorithms/class_simulation_state.h>
#include <mps_routines/mps_tools/finite/opt.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <model/class_hamiltonian_base.h>
#include <LBFGS.h>



MPS_Tools::Finite::Opt::internals::base_functor::base_functor(
        const class_superblock & superblock,
        class_simulation_state & sim_state)
{
    reset_timers();
    initialize_params();
    superComponents.HA_MPO        = superblock.HA->MPO().real();
    superComponents.HB_MPO        = superblock.HB->MPO().real();
    superComponents.Lblock        = superblock.Lblock->block.real();
    superComponents.Rblock        = superblock.Rblock->block.real();
    superComponents.Lblock2       = superblock.Lblock2->block.real();
    superComponents.Rblock2       = superblock.Rblock2->block.real();
    superComponents.dsizes        = superblock.dimensions();
    superComponents.HAHB          = superComponents.HA_MPO.contract(superComponents.HB_MPO, Textra::idx({1},{0}));
    superComponents.HAHB2         = superComponents.HAHB.contract(superComponents.HAHB, Textra::idx({2,5},{1,4}));
    length                        = superblock.get_length();

    //All energies in sim_state are per site!
    energy_target            = sim_state.energy_target;
    energy_max               = sim_state.energy_max;
    energy_min               = sim_state.energy_min;
    energy_lower_bound       = sim_state.energy_lbound;
    energy_upper_bound       = sim_state.energy_ubound;
    energy_target_dens       = (energy_target - energy_min ) / (energy_max - energy_min);
    energy_window            = 0.5*(energy_upper_bound - energy_lower_bound ) / (energy_max - energy_min);

    iteration                = sim_state.iteration;
}


double MPS_Tools::Finite::Opt::internals::base_functor::get_variance()const{return variance;}
double MPS_Tools::Finite::Opt::internals::base_functor::get_energy  ()const{return energy  ;}
size_t MPS_Tools::Finite::Opt::internals::base_functor::get_count   ()const{return counter;}



std::pair<Eigen::VectorXd,double> MPS_Tools::Finite::Opt::internals::base_functor::get_vH_vHv(const Eigen::Matrix<double,Eigen::Dynamic,1> &v){
    auto vH = MPS_Tools::Finite::Opt::internals::base_functor::get_vH(v);
    t_vHv->tic();
    auto vHv = vH.dot(v);
    t_vHv->toc();
    return std::make_pair(vH,vHv);
}

std::pair<Eigen::VectorXd,double> MPS_Tools::Finite::Opt::internals::base_functor::get_vH2_vH2v(const Eigen::Matrix<double,Eigen::Dynamic,1> &v){
    auto vH2 = MPS_Tools::Finite::Opt::internals::base_functor::get_vH2(v);
    t_vH2v->tic();
    auto vH2v = vH2.dot(v);
    t_vH2v->toc();
    return std::make_pair(vH2,vH2v);
}


Eigen::VectorXd MPS_Tools::Finite::Opt::internals::base_functor::get_vH2 (const Eigen::Matrix<double,Eigen::Dynamic,1> &v){
    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), superComponents.dsizes);
    t_vH2->tic();
    Eigen::Tensor<double, 4> vH2 =
            superComponents.Lblock2
                    .contract(theta,                            Textra::idx({0},{1}))
                    .contract(superComponents.HAHB2,            Textra::idx({2,1,3,4},{4,0,1,3}))
                    .contract(superComponents.Rblock2,          Textra::idx({1,2,4},{0,2,3}))
                    .shuffle(Textra::array4{1,0,2,3});
    t_vH2->toc();
    return Eigen::Map<Eigen::VectorXd>(vH2.data(),vH2.size());
}

Eigen::VectorXd MPS_Tools::Finite::Opt::internals::base_functor::get_vH (const Eigen::Matrix<double,Eigen::Dynamic,1> &v){
    auto theta = Eigen::TensorMap<const Eigen::Tensor<const double,4>> (v.data(), superComponents.dsizes);
    t_vH->tic();
    Eigen::Tensor<double, 4> vH =
            superComponents.Lblock
                    .contract(theta,                               Textra::idx({0},{1}))
                    .contract(superComponents.HAHB,                Textra::idx({1,2,3},{0,1,4}))
                    .contract(superComponents.Rblock ,             Textra::idx({1,3},{0,2}))
                    .shuffle(Textra::array4{1,0,2,3});
    t_vH->toc();

    return Eigen::Map<Eigen::VectorXd>(vH.data(),vH.size());
}




