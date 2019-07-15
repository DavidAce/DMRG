
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_tools_finite_profile.cpp:

Program Listing for File profile.cpp
====================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_tools_finite_profile.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/state/tools/finite/profile.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-07-05.
   //
   
   #include <state/tools/nmspc_tools.h>
   #include <simulation/nmspc_settings.h>
   
   
   void tools::finite::profile::print_profiling(class_tic_toc &t_parent){
       if (settings::profiling::on) {
           std::cout << "\nComputing observables breakdown:" << std::endl;
           std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
           t_ene.print_time_w_percent(t_parent);
           t_var.print_time_w_percent(t_parent);
           t_ent.print_time_w_percent(t_parent);
   
       }
   }
   
   
   void tools::finite::profile::init_profiling(){
       t_eig.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Eig. decomp.      ");
       t_svd.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Svd. decomp.      ");
       t_ene.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Energy            ");
       t_var.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Variance          ");
       t_ent.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Ent. Entr.        ");
       t_hdf.set_properties(settings::profiling::on, settings::profiling::precision, "↳ h5pp storage      ");
       t_prj.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Projections       ");
       t_dbg.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Debugging         ");
   }
