
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_common_prof.cpp:

Program Listing for File prof.cpp
=================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_common_prof.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/common/prof.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-06-08.
   //
   
   #include <tools/nmspc_tools.h>
   #include <simulation/nmspc_settings.h>
   
   
   
   
   void tools::common::profile::print_profiling(){
       if (settings::profiling::on) {
           t_tot    .print_time();
           t_pre    .print_time_w_percent_if_nonzero(t_tot);
           t_pos    .print_time_w_percent_if_nonzero(t_tot);
           t_sim    .print_time_w_percent_if_nonzero(t_tot);
           t_con    .print_time_w_percent_if_nonzero(t_sim);
           t_eig    .print_time_w_percent_if_nonzero(t_sim);
           t_svd    .print_time_w_percent_if_nonzero(t_sim);
           t_opt    .print_time_w_percent_if_nonzero(t_sim);
           t_evo    .print_time_w_percent_if_nonzero(t_sim);
           t_env    .print_time_w_percent_if_nonzero(t_sim);
           t_ent    .print_time_w_percent_if_nonzero(t_sim);
           t_ene    .print_time_w_percent_if_nonzero(t_sim);
           t_var    .print_time_w_percent_if_nonzero(t_sim);
           t_prj    .print_time_w_percent_if_nonzero(t_sim);
           t_chk    .print_time_w_percent_if_nonzero(t_sim);
           t_hdf    .print_time_w_percent_if_nonzero(t_sim);
           t_ene_ham.print_time_w_percent_if_nonzero(t_sim);
           t_ene_mom.print_time_w_percent_if_nonzero(t_sim);
           t_var_ham.print_time_w_percent_if_nonzero(t_sim);
           t_var_mom.print_time_w_percent_if_nonzero(t_sim);
       }
   }
   
   
   void tools::common::profile::init_profiling(){
   
       t_tot.set_properties    (settings::profiling::on, settings::profiling::precision,"+Total Time              ");
       t_pre.set_properties    (settings::profiling::on, settings::profiling::precision,"↳ Preprocessing          ");
       t_pos.set_properties    (settings::profiling::on, settings::profiling::precision,"↳ Postprocessing         ");
       t_sim.set_properties    (settings::profiling::on, settings::profiling::precision, "↳+Simulation             ");
       t_con.set_properties    (settings::profiling::on, settings::profiling::precision," ↳ Convergence checks    ");
       t_eig.set_properties    (settings::profiling::on, settings::profiling::precision," ↳ Eig. decomp.          ");
       t_svd.set_properties    (settings::profiling::on, settings::profiling::precision," ↳ Svd. decomp.          ");
       t_opt.set_properties    (settings::profiling::on, settings::profiling::precision," ↳ Optimization (Ceres)  ");
       t_evo.set_properties    (settings::profiling::on, settings::profiling::precision," ↳ Time evolution        ");
       t_env.set_properties    (settings::profiling::on, settings::profiling::precision," ↳ Environment upd.      ");
       t_ent.set_properties    (settings::profiling::on, settings::profiling::precision," ↳ Entanglement entropy  ");
       t_ene.set_properties    (settings::profiling::on, settings::profiling::precision," ↳ Energy                ");
       t_var.set_properties    (settings::profiling::on, settings::profiling::precision," ↳ Variance              ");
       t_prj.set_properties    (settings::profiling::on, settings::profiling::precision," ↳ Projections           ");
       t_chk.set_properties    (settings::profiling::on, settings::profiling::precision," ↳ Checks                ");
       t_hdf.set_properties    (settings::profiling::on, settings::profiling::precision," ↳ h5pp storage          ");
       t_ene_ham.set_properties(settings::profiling::on, settings::profiling::precision," ↳ Energy (HAM)          ");
       t_ene_mom.set_properties(settings::profiling::on, settings::profiling::precision," ↳ Energy (MOM)          ");
       t_var_ham.set_properties(settings::profiling::on, settings::profiling::precision," ↳ Variance (HAM)        ");
       t_var_mom.set_properties(settings::profiling::on, settings::profiling::precision," ↳ Variance (MOM)        ");
   
   }
