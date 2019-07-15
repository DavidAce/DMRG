
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_simulation_class_simulation_status.h:

Program Listing for File class_simulation_status.h
==================================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_simulation_class_simulation_status.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/simulation/class_simulation_status.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-02-14.
   //
   
   #ifndef DMRG_CLASS_SIMULATION_STATUS_H
   #define DMRG_CLASS_SIMULATION_STATUS_H
   
   #include <memory>
   #include <string>
   #include <iostream>
   #include <vector>
   #include <array>
   #include <hdf5.h>
   #include <hdf5_hl.h>
   
   struct status_data{
       // common variables
       size_t iteration                      = 0; //In idmrg and itebd: iterations, in fdmrg and xdmrg: full sweeps along the state.
       size_t step                           = 0; //In fdmrg and xdmrg: how many individual moves along the state.
       size_t position                       = 0;
       long   chi_temp                       = 16;
       long   chi_max                        = 16;
       size_t min_sweeps                     = 2 ;
       double energy_min                     = 0;
       double energy_max                     = 0;
       double energy_target                  = 0;
       double energy_ubound                  = 0;
       double energy_lbound                  = 0;
       double energy_dens                    = 0;
       double energy_dens_target             = 0;
       double energy_dens_window             = 0;
       double phys_time                      = 0;
       double wall_time                      = 0;
       double simu_time                      = 0;
       double delta_t                        = 0; //Make sure this one gets initialized to delta_t0!
       bool   time_step_has_converged        = false;
       bool   simulation_has_converged       = false;
       bool   simulation_has_saturated       = false;
       bool   simulation_has_to_stop         = false;
       bool   bond_dimension_has_reached_max = false;
       bool   entanglement_has_converged     = false;
       bool   entanglement_has_saturated     = false;
       bool   variance_mpo_has_converged     = false;
       bool   variance_mpo_has_saturated     = false;
       bool   variance_ham_has_converged     = false;
       bool   variance_ham_has_saturated     = false;
       bool   variance_mom_has_converged     = false;
       bool   variance_mom_has_saturated     = false;
       size_t variance_mpo_saturated_for     = 0;
       size_t variance_ham_saturated_for     = 0;
       size_t variance_mom_saturated_for     = 0;
   };
   
   class class_simulation_status : public status_data{
   public:
   
       void clear();
   //    void get_all();
       friend std::ostream& operator <<(std::ostream& os, const class_simulation_status & sim_status);
   
   };
   
   
   
   #endif //DMRG_CLASS_SIMULATION_STATUS_H
