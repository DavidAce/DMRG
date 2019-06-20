//
// Created by david on 2019-06-02.
//


#include <mps_tools/finite/opt.h>
#include <algorithms/class_simulation_state.h>
#include <mps_state/class_superblock.h>
#include <mps_state/class_environment.h>
#include <model/class_hamiltonian_base.h>

template<typename Scalar>
mpstools::finite::opt::internals::cppoptlib_functor<Scalar>::cppoptlib_functor(const class_superblock &superblock,
                                                                        const class_simulation_state &sim_state)
                                                                        : superComponents(superblock)
                                                                        {
    reset_timers();
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

template<typename Scalar>
void mpstools::finite::opt::internals::cppoptlib_functor<Scalar>::compute_exp_vals(const VectorType &v){
    if(not exp_vals_computed){
        counter ++;
        #pragma omp parallel
        {
            #pragma omp sections
            {
                #pragma omp section
                { std::tie(vH2, vH2v) = get_H2v_vH2v(v, superComponents); }
                #pragma omp section
                { std::tie(vH, vHv) = get_Hv_vHv(v, superComponents); }
                #pragma omp section
                { vv = v.squaredNorm(); }
            }
        }
    }
//    exp_vals_computed = true;
}

template<typename Scalar>
bool mpstools::finite::opt::internals::cppoptlib_functor<Scalar>::callback(const cppoptlib::Criteria<Scalar> &state, [[maybe_unused]] const VectorType &x0){
    iteration = state.iterations;
    exp_vals_computed = false;
    std::cout << "CALLBACK" << std::endl;
    return true;
}

template<typename Scalar>
Scalar mpstools::finite::opt::internals::cppoptlib_functor<Scalar>::value(const VectorType &v){
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

template<typename Scalar>
void mpstools::finite::opt::internals::cppoptlib_functor<Scalar>::gradient(const VectorType &v,  VectorType &grad){
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





template class mpstools::finite::opt::internals::cppoptlib_functor<double>;
template class mpstools::finite::opt::internals::cppoptlib_functor<std::complex<double>>;

