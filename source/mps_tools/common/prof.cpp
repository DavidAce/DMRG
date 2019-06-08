//
// Created by david on 2019-06-08.
//

#include <mps_state/nmspc_mps_tools.h>
#include <sim_parameters/nmspc_sim_settings.h>

void MPS_Tools::Common::Prof::Obs:: print_profiling(class_tic_toc &t_parent){
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
void MPS_Tools::Common::Prof::enable_profiling(int precision){
    using namespace settings::profiling;
    Obs::t_ene_mpo.set_properties(on, precision,"↳ Energy (MPO)           ");
    Obs::t_ene_ham.set_properties(on, precision,"↳ Energy (HAM)           ");
    Obs::t_ene_mom.set_properties(on, precision,"↳ Energy (MOM)           ");
    Obs::t_var_mpo.set_properties(on, precision,"↳ Variance (MPO)         ");
    Obs::t_var_ham.set_properties(on, precision,"↳ Variance (HAM)         ");
    Obs::t_var_mom.set_properties(on, precision,"↳ Variance (MOM)         ");
    Obs::t_entropy.set_properties(on, precision,"↳ Ent. Entropy           ");
    Obs::t_temp1.set_properties  (on, precision,"↳ Temp1                  ");
    Obs::t_temp2.set_properties  (on, precision,"↳ Temp2                  ");
    Obs::t_temp3.set_properties  (on, precision,"↳ Temp3                  ");
    Obs::t_temp4.set_properties  (on, precision,"↳ Temp4                  ");
}