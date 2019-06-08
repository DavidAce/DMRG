//
// Created by david on 2019-06-02.
//


#include <mps_tools/finite/opt.h>
#include <algorithms/class_simulation_state.h>
#include <mps_state/class_superblock.h>
#include <mps_state/class_environment.h>
#include <model/class_hamiltonian_base.h>


MPS_Tools::Finite::Opt::internals::cppoptlib_functor::cppoptlib_functor(const class_superblock &superblock,
                                                                        const class_simulation_state &sim_state){
    reset_timers();
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


double MPS_Tools::Finite::Opt::internals::cppoptlib_functor::get_variance()const{return variance;}
double MPS_Tools::Finite::Opt::internals::cppoptlib_functor::get_energy  ()const{return energy  ;}
size_t MPS_Tools::Finite::Opt::internals::cppoptlib_functor::get_count   ()const{return counter;}
size_t MPS_Tools::Finite::Opt::internals::cppoptlib_functor::get_iter    ()const{return iteration;}




void MPS_Tools::Finite::Opt::internals::cppoptlib_functor::compute_exp_vals(const Eigen::VectorXd &v){
    if(not exp_vals_computed){
        counter ++;
        #pragma omp parallel
        {
            #pragma omp sections
            {
                #pragma omp section
                { std::tie(vH2, vH2v) = get_vH2_vH2v(v,superComponents); }
                #pragma omp section
                { std::tie(vH, vHv) = get_vH_vHv(v,superComponents); }
                #pragma omp section
                { vv = v.squaredNorm(); }
            }
        }
    }
//    exp_vals_computed = true;
}

bool MPS_Tools::Finite::Opt::internals::cppoptlib_functor::callback(const cppoptlib::Criteria<double> &state, const Eigen::VectorXd &x0){
    iteration = state.iterations;
    exp_vals_computed = false;
    std::cout << "CALLBACK" << std::endl;
    return true;
}


double MPS_Tools::Finite::Opt::internals::cppoptlib_functor::value(const Eigen::VectorXd &v){
    compute_exp_vals(v);
//    exp_vals_computed = false;
    std::cout << "VALUE" << std::endl;

    auto lambda         = 1.0;
    var                 = std::abs(vH2v*vv - vHv*vHv)/std::pow(vv,2);
    variance            = var;
    var                 = var == 0  ? std::numeric_limits<double>::epsilon() : var;
    energy              = vHv/vv;
    norm_offset         = vv-1.0;
    auto log10var       = std::log10(var);
    return log10var  + lambda * std::pow(norm_offset,2);
}

void MPS_Tools::Finite::Opt::internals::cppoptlib_functor::gradient(const Eigen::VectorXd &v,  Eigen::VectorXd &grad){
    compute_exp_vals(v);
    std::cout << "GRADIENT" << std::endl;
    auto lambda         = 1.0;
    var                 = std::abs(vH2v*vv - vHv*vHv)/std::pow(vv,2);
    variance            = var;
    var                 = var == 0  ? std::numeric_limits<double>::epsilon() : var;
    energy              = vHv/vv;
    norm_offset         = vv-1.0;
    auto vv_1           = std::pow(vv,-1);
    auto vv_2           = std::pow(vv,-2);
    auto var_1          = 1.0/var/std::log(10);
    grad                = var_1 * (2.0*(vH2*vv_1 - v * vH2v * vv_2) - 4.0 * energy * (vH * vv_1 - v * vHv * vv_2))
                           + lambda * 2.0 * norm_offset * 2.0 * v;

}