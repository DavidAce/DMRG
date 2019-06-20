//
// Created by david on 2019-06-08.
//

#include <mps_tools/nmspc_mps_tools.h>
#include <sim_parameters/nmspc_sim_settings.h>

void mpstools::common::profiling::obs:: print_profiling(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        std::cout << "\nComputing observables breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_ene_mpo.print_time_w_percent(t_parent);
        t_ene_ham.print_time_w_percent(t_parent);
        t_ene_mom.print_time_w_percent(t_parent);
        t_var_mpo.print_time_w_percent(t_parent);
        t_var_ham.print_time_w_percent(t_parent);
        t_var_mom.print_time_w_percent(t_parent);
        t_entropy.print_time_w_percent(t_parent);
        t_temp1.print_time_w_percent(t_parent);
        t_temp2.print_time_w_percent(t_parent);
        t_temp3.print_time_w_percent(t_parent);
        t_temp4.print_time_w_percent(t_parent);
    }
}
void mpstools::common::profiling::init_profiling(bool on, int precision){
    using namespace settings::profiling;
    obs::t_ene_mpo.set_properties(on, precision,"↳ Energy (MPO)           ");
    obs::t_ene_ham.set_properties(on, precision,"↳ Energy (HAM)           ");
    obs::t_ene_mom.set_properties(on, precision,"↳ Energy (MOM)           ");
    obs::t_var_mpo.set_properties(on, precision,"↳ Variance (MPO)         ");
    obs::t_var_ham.set_properties(on, precision,"↳ Variance (HAM)         ");
    obs::t_var_mom.set_properties(on, precision,"↳ Variance (MOM)         ");
    obs::t_entropy.set_properties(on, precision,"↳ Ent. Entropy           ");
    obs::t_temp1.set_properties  (on, precision,"↳ Temp1                  ");
    obs::t_temp2.set_properties  (on, precision,"↳ Temp2                  ");
    obs::t_temp3.set_properties  (on, precision,"↳ Temp3                  ");
    obs::t_temp4.set_properties  (on, precision,"↳ Temp4                  ");
}