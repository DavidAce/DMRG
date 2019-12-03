
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_common_io_h5pp_tables.cpp:

Program Listing for File h5pp_tables.cpp
========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_common_io_h5pp_tables.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/common/io/h5pp_tables.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-11-07.
   //
   #include <tools/nmspc_tools.h>
   #include <io/table_types.h>
   #include <io/class_h5table_buffer.h>
   
   void tools::common::io::h5table::write_sim_status(const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_simulation_status> &h5tbuf) {
       log->trace("Appending simulation status table: {}...",h5tbuf.get_table_name());
       h5tbuf.append_record(sim_status);
       log->trace("Appending simulation status table: {}... OK",h5tbuf.get_table_name());
   }
   
   
   void tools::common::io::h5table::write_profiling(const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_profiling> &h5tbuf) {
       log->trace("Appending profiling data to table: {}...",h5tbuf.get_table_name());
       class_h5table_profiling::data_struct profiling_entry;
       profiling_entry.step            = sim_status.step;
       profiling_entry.iteration       = sim_status.iteration;
       profiling_entry.position        = sim_status.position;
       profiling_entry.t_tot           = sim_status.wall_time;
       profiling_entry.t_run           = sim_status.simu_time;
       profiling_entry.t_eig           = tools::common::profile::t_eig.get_measured_time();
       profiling_entry.t_svd           = tools::common::profile::t_svd.get_measured_time();
       profiling_entry.t_ene           = tools::common::profile::t_ene.get_measured_time();
       profiling_entry.t_var           = tools::common::profile::t_var.get_measured_time();
       profiling_entry.t_ent           = tools::common::profile::t_ent.get_measured_time();
       profiling_entry.t_hdf           = tools::common::profile::t_hdf.get_measured_time();
       profiling_entry.t_prj           = tools::common::profile::t_prj.get_measured_time();
       profiling_entry.t_opt           = tools::common::profile::t_opt.get_measured_time();
       profiling_entry.t_chk           = tools::common::profile::t_chk.get_measured_time();
       profiling_entry.t_ene_mpo       = tools::common::profile::t_ene_mpo.get_measured_time();
       profiling_entry.t_ene_ham       = tools::common::profile::t_ene_ham.get_measured_time();
       profiling_entry.t_ene_mom       = tools::common::profile::t_ene_mom.get_measured_time();
       profiling_entry.t_var_mpo       = tools::common::profile::t_var_mpo.get_measured_time();
       profiling_entry.t_var_ham       = tools::common::profile::t_var_ham.get_measured_time();
       profiling_entry.t_var_mom       = tools::common::profile::t_var_mom.get_measured_time();
   
   
       profiling_entry.t_env = 0;
       profiling_entry.t_evo = 0;
       profiling_entry.t_udt = 0;
       profiling_entry.t_ste = 0;
       profiling_entry.t_prt = 0;
       profiling_entry.t_obs = 0;
       profiling_entry.t_mps = 0;
       profiling_entry.t_chi = 0;
   
       h5tbuf.append_record(profiling_entry);
       log->trace("Appending profiling data to table: {}... OK",h5tbuf.get_table_name());
   
   }
