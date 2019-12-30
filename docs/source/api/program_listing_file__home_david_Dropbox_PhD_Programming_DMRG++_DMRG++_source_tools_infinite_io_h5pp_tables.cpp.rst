
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_infinite_io_h5pp_tables.cpp:

Program Listing for File h5pp_tables.cpp
========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_infinite_io_h5pp_tables.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/infinite/io/h5pp_tables.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-11-07.
   //
   
   #include <tools/nmspc_tools.h>
   #include <io/class_h5table_buffer.h>
   #include <state/class_state_infinite.h>
   
   
   void tools::infinite::io::h5table::write_measurements(const class_state_infinite &state, const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_measurements_infinite> &h5tbuf) {
       log->trace("Appending measurement entry to table: {}...",h5tbuf.get_table_name());
       class_h5table_measurements_infinite::data_struct measurements_entry;
       measurements_entry.step                            = sim_status.step;
       measurements_entry.iteration                       = sim_status.iteration;
       measurements_entry.position                        = sim_status.position;
       measurements_entry.length                          = tools::infinite::measure::length(state);
       measurements_entry.bond_dimension                  = tools::infinite::measure::bond_dimension(state);
       measurements_entry.bond_dimension_limit            = state.get_chi_lim();
       measurements_entry.bond_dimension_maximum          = state.get_chi_max();
       measurements_entry.entanglement_entropy            = tools::infinite::measure::entanglement_entropy(state);
       measurements_entry.norm                            = tools::infinite::measure::norm(state);
       measurements_entry.energy_mpo                      = tools::infinite::measure::energy_mpo(state);
       measurements_entry.energy_per_site_mpo             = tools::infinite::measure::energy_per_site_mpo(state);
       measurements_entry.energy_per_site_ham             = tools::infinite::measure::energy_per_site_ham(state);
       measurements_entry.energy_per_site_mom             = tools::infinite::measure::energy_per_site_mom(state);
       measurements_entry.energy_variance_mpo             = tools::infinite::measure::energy_variance_mpo(state);
       measurements_entry.energy_variance_per_site_mpo    = tools::infinite::measure::energy_variance_per_site_mpo(state);
       measurements_entry.energy_variance_per_site_ham    = tools::infinite::measure::energy_variance_per_site_ham(state);
       measurements_entry.energy_variance_per_site_mom    = tools::infinite::measure::energy_variance_per_site_mom(state);
       measurements_entry.truncation_error                = tools::infinite::measure::truncation_error(state);
       measurements_entry.wall_time                       = sim_status.wall_time;
       measurements_entry.wall_time                       = sim_status.phys_time;
       measurements_entry.wall_time                       = sim_status.delta_t;
       h5tbuf.append_record(measurements_entry);
       log->trace("Appending measurement entry to table: {}... OK",h5tbuf.get_table_name());
   }
   
   
   
   
   void tools::infinite::io::h5table::write_sim_status(const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_simulation_status> &h5tbuf) {
       tools::common::io::h5table::write_sim_status(sim_status,h5tbuf);
   }
   
   void tools::infinite::io::h5table::write_profiling(const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_profiling> &h5tbuf) {
       tools::common::io::h5table::write_profiling(sim_status,h5tbuf);
   }
