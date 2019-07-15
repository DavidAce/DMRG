
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_io_log_types.h:

Program Listing for File log_types.h
====================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_io_log_types.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/io/log_types.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2018-05-24.
   //
   
   #ifndef DMRG_LOG_TYPES_H
   #define DMRG_LOG_TYPES_H
   
   #include <vector>
   #include <array>
   #include <hdf5.h>
   #include <hdf5_hl.h>
   #include <simulation/class_simulation_status.h>
   
   
   
   class class_log_dmrg {
   private:
       struct data_struct {
           int     iteration;
           int     chain_length;
           int     position;
           long    chi;
           long    chi_max;
           double  energy_mpo; double  energy_ham; double  energy_mom;
           double  energy_min; double  energy_max; double  energy_tgt;
           double  variance_mpo; double  variance_ham; double  variance_mom;
           double  entanglement_entropy_midchain;
           double  truncation_error;
           double  wall_time;
   
           data_struct(
           int iteration_,
           int chain_length_,
           int position_,
           long chi_,
           long chi_max_,
           double energy_mpo_, double energy_ham_, double energy_mom_,
           double energy_min_, double energy_max_, double energy_tgt_,
           double variance_mpo_, double variance_ham_, double variance_mom_,
           double entropy_,
           double truncation_error_,
           double wall_time_) :
                   iteration(iteration_),
                   chain_length(chain_length_),
                   position(position_),
                   chi(chi_),
                   chi_max(chi_max_),
                   energy_mpo(energy_mpo_), energy_ham(energy_ham_), energy_mom(energy_mom_),
                   energy_min(energy_min_), energy_max(energy_max_), energy_tgt(energy_tgt_),
                   variance_mpo(variance_mpo_), variance_ham(variance_ham_), variance_mom(variance_mom_),
                   entanglement_entropy_midchain(entropy_),
                   truncation_error(truncation_error_),
                   wall_time(wall_time_)
           {}
       };
       struct meta_struct {
           constexpr static hsize_t NFIELDS = 17;
           size_t dst_size = sizeof(data_struct);
           std::array<size_t, NFIELDS> dst_offsets = {HOFFSET(data_struct, iteration),
                                                      HOFFSET(data_struct, chain_length),
                                                      HOFFSET(data_struct, position),
                                                      HOFFSET(data_struct, chi),
                                                      HOFFSET(data_struct, chi_max),
                                                      HOFFSET(data_struct, energy_mpo), HOFFSET(data_struct, energy_ham), HOFFSET(data_struct, energy_mom),
                                                      HOFFSET(data_struct, energy_min), HOFFSET(data_struct, energy_max), HOFFSET(data_struct, energy_tgt),
                                                      HOFFSET(data_struct, variance_mpo), HOFFSET(data_struct, variance_ham), HOFFSET(data_struct, variance_mom),
                                                      HOFFSET(data_struct, entanglement_entropy_midchain),
                                                      HOFFSET(data_struct, truncation_error),
                                                      HOFFSET(data_struct, wall_time)
           };
           std::array<size_t, NFIELDS> dst_sizes = {
                   sizeof(data_struct::iteration),
                   sizeof(data_struct::chain_length),
                   sizeof(data_struct::position),
                   sizeof(data_struct::chi),
                   sizeof(data_struct::chi_max),
                   sizeof(data_struct::energy_mpo), sizeof(data_struct::energy_ham), sizeof(data_struct::energy_mom),
                   sizeof(data_struct::energy_min), sizeof(data_struct::energy_max), sizeof(data_struct::energy_tgt),
                   sizeof(data_struct::variance_mpo), sizeof(data_struct::variance_ham), sizeof(data_struct::variance_mom),
                   sizeof(data_struct::entanglement_entropy_midchain), sizeof(data_struct::truncation_error), sizeof(data_struct::wall_time)
           };
           std::array<const char *, NFIELDS> field_names = {
                   "iteration",
                   "chain_length",
                   "position",
                   "chi",
                   "chi_max",
                   "energy","energy_per_site_ham","energy_per_site_mom",
                   "energy_min","energy_max","energy_tgt",
                   "variance_mpo","variance_ham","variance_mom",
                   "entanglement_entropy_midchain",
                   "truncation_error",
                   "wall_time"
           };
   
           std::array<hid_t, NFIELDS> field_types = {H5T_NATIVE_INT,
                                                     H5T_NATIVE_INT,
                                                     H5T_NATIVE_INT,
                                                     H5T_NATIVE_LONG,
                                                     H5T_NATIVE_LONG,
                                                     H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                     H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                     H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                     H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE};
           hsize_t chunk_size = 1;
           void *fill_data = nullptr;
           int compress = 0;
       };
   public:
       class_log_dmrg() = default;
       meta_struct meta;
       std::vector<data_struct> buffer;
   
   };
   
   
   class class_log_tebd{
   private:
       struct data_struct {
           int     iteration;
           long    chi;
           long    chi_max;
           double  time_step;
           double  energy_mpo; double  energy_ham; double  energy_mom;
           double  variance_mpo; double  variance_ham; double  variance_mom;
           double  entanglement_entropy;
           double  truncation_error;
           double  phys_time;
           double  wall_time;
   
           data_struct(
                   int    iteration_,
                   long   chi_,
                   long   chi_max_,
                   double time_step_,
                   double energy_mpo_, double energy_ham_, double energy_mom_,
                   double variance_mpo_, double variance_ham_, double variance_mom_,
                   double entropy_,
                   double truncation_error_,
                   double phys_time_,
                   double wall_time_) :
                   iteration(iteration_),
                   chi(chi_),
                   chi_max(chi_max_),
                   time_step(time_step_),
                   energy_mpo(energy_mpo_), energy_ham(energy_ham_), energy_mom(energy_mom_),
                   variance_mpo(variance_mpo_), variance_ham(variance_ham_), variance_mom(variance_mom_),
                   entanglement_entropy(entropy_),
                   truncation_error(truncation_error_),
                   phys_time(phys_time_),
                   wall_time(wall_time_)
           {}
       };
       struct meta_struct {
           constexpr static hsize_t NFIELDS = 14;
           size_t dst_size = sizeof(data_struct);
           std::array<size_t, NFIELDS> dst_offsets = {HOFFSET(data_struct, iteration),
                                                      HOFFSET(data_struct, chi),
                                                      HOFFSET(data_struct, chi_max),
                                                      HOFFSET(data_struct, time_step),
                                                      HOFFSET(data_struct, energy_mpo), HOFFSET(data_struct, energy_ham), HOFFSET(data_struct, energy_mom),
                                                      HOFFSET(data_struct, variance_mpo), HOFFSET(data_struct, variance_ham), HOFFSET(data_struct, variance_mom),
                                                      HOFFSET(data_struct, entanglement_entropy),
                                                      HOFFSET(data_struct, truncation_error),
                                                      HOFFSET(data_struct, phys_time),
                                                      HOFFSET(data_struct, wall_time)
           };
           std::array<size_t, NFIELDS> dst_sizes = {
                   sizeof(data_struct::iteration),
                   sizeof(data_struct::chi),
                   sizeof(data_struct::chi_max),
                   sizeof(data_struct::time_step),
                   sizeof(data_struct::energy_mpo), sizeof(data_struct::energy_ham), sizeof(data_struct::energy_mom),
                   sizeof(data_struct::variance_mpo), sizeof(data_struct::variance_ham), sizeof(data_struct::variance_mom),
                   sizeof(data_struct::entanglement_entropy),
                   sizeof(data_struct::truncation_error),
                   sizeof(data_struct::phys_time),
                   sizeof(data_struct::wall_time)
           };
           std::array<const char *, NFIELDS> field_names = {"iteration",
                                                            "chi",
                                                            "chi_max",
                                                            "time_step",
                                                            "energy","energy_per_site_ham","energy_per_site_mom",
                                                            "variance_mpo","variance_ham","variance_mom",
                                                            "entanglement_entropy_midchain",
                                                            "truncation_error",
                                                            "phys_time",
                                                            "wall_time"
           };
   
           std::array<hid_t, NFIELDS> field_types = {H5T_NATIVE_INT,
                                                     H5T_NATIVE_LONG,
                                                     H5T_NATIVE_LONG,
                                                     H5T_NATIVE_DOUBLE,
                                                     H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                     H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                     H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                     H5T_NATIVE_DOUBLE};
           hsize_t chunk_size = 1;
           void *fill_data = nullptr;
           int compress = 0;
       };
   public:
       class_log_tebd() = default;
       meta_struct meta;
       std::vector<data_struct> buffer;
   };
   
   
   
   class class_log_profiling{
   private:
       struct data_struct{
           int    iteration;
           int    step;
           double t_tot;
           double t_opt;
           double t_sim;
           double t_svd;
           double t_env;
           double t_evo;
           double t_udt;
           double t_sto;
           double t_ste;
           double t_prt;
           double t_obs;
           double t_mps;
           double t_chi;
   
           data_struct(
                   int    iteration_,
                   int    step_,
                double t_tot_,  double t_opt_, double t_sim_,
                double t_svd_,  double t_env_, double t_evo_,
                double t_udt_,  double t_sto_, double t_ste_,
                double t_prt_,  double t_obs_, double t_mps_,
                double t_chi_)
                   :iteration(iteration_),step(step_),
                    t_tot(t_tot_), t_opt(t_opt_), t_sim(t_sim_),
                    t_svd(t_svd_), t_env(t_env_), t_evo(t_evo_),
                    t_udt(t_udt_), t_sto(t_sto_), t_ste(t_ste_),
                    t_prt(t_prt_), t_obs(t_obs_), t_mps(t_mps_),
                    t_chi(t_chi_)
           {}
       };
       struct meta_struct{
           constexpr static hsize_t                NFIELDS     = 15;
           size_t           dst_size                           = sizeof (data_struct);
           std::array       <size_t,NFIELDS>       dst_offsets = {HOFFSET(data_struct, iteration), HOFFSET(data_struct, step),
                                                                  HOFFSET(data_struct, t_tot), HOFFSET(data_struct, t_opt), HOFFSET(data_struct, t_sim),
                                                                  HOFFSET(data_struct, t_svd), HOFFSET(data_struct, t_env), HOFFSET(data_struct, t_evo),
                                                                  HOFFSET(data_struct, t_udt), HOFFSET(data_struct, t_sto), HOFFSET(data_struct, t_ste),
                                                                  HOFFSET(data_struct, t_prt), HOFFSET(data_struct, t_obs), HOFFSET(data_struct, t_mps),
                                                                  HOFFSET(data_struct, t_chi)
           };
           std::array       <size_t,NFIELDS>       dst_sizes   = {sizeof(data_struct::iteration), sizeof(data_struct::step),
                                                                  sizeof(data_struct::t_tot), sizeof(data_struct::t_opt), sizeof(data_struct::t_sim),
                                                                  sizeof(data_struct::t_svd), sizeof(data_struct::t_env), sizeof(data_struct::t_evo),
                                                                  sizeof(data_struct::t_udt), sizeof(data_struct::t_sto), sizeof(data_struct::t_ste),
                                                                  sizeof(data_struct::t_prt), sizeof(data_struct::t_obs), sizeof(data_struct::t_mps),
                                                                  sizeof(data_struct::t_chi)
           };
           std::array       <const char*,NFIELDS>  field_names = {"iteration","step",
                                                                  "t_tot", "t_opt", "t_sim",
                                                                  "t_svd", "t_env", "t_evo",
                                                                  "t_udt", "t_sto", "t_ste",
                                                                  "t_prt", "t_obs", "t_mps",
                                                                  "t_con"
           };
   
           std::array       <hid_t,NFIELDS>        field_types = {H5T_NATIVE_INT, H5T_NATIVE_INT,
                                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                                  H5T_NATIVE_DOUBLE
           };
   
           hsize_t          chunk_size                         = 1;
           void             *fill_data                         = nullptr;
           int              compress                           = 0;
       };
   public:
       class_log_profiling() = default;
       meta_struct meta;
       std::vector<data_struct> buffer;
   };
   
   
   class class_log_simulation_status{
   private:
       struct meta_struct{
           constexpr static hsize_t                NFIELDS     = 34;
           size_t           dst_size                           = sizeof (status_data);
           std::array       <size_t,NFIELDS>       dst_offsets =
                   {
                       HOFFSET(status_data, iteration                     ),
                       HOFFSET(status_data, step                          ),
                       HOFFSET(status_data, position                      ),
                       HOFFSET(status_data, chi_temp                      ),
                       HOFFSET(status_data, chi_max                       ),
                       HOFFSET(status_data, min_sweeps                    ),
                       HOFFSET(status_data, energy_min                    ),
                       HOFFSET(status_data, energy_max                    ),
                       HOFFSET(status_data, energy_target                 ),
                       HOFFSET(status_data, energy_ubound                 ),
                       HOFFSET(status_data, energy_lbound                 ),
                       HOFFSET(status_data, energy_dens                   ),
                       HOFFSET(status_data, energy_dens_target            ),
                       HOFFSET(status_data, energy_dens_window            ),
                       HOFFSET(status_data, phys_time                     ),
                       HOFFSET(status_data, wall_time                     ),
                       HOFFSET(status_data, simu_time                     ),
                       HOFFSET(status_data, delta_t                       ),
                       HOFFSET(status_data, time_step_has_converged       ),
                       HOFFSET(status_data, simulation_has_converged      ),
                       HOFFSET(status_data, simulation_has_saturated      ),
                       HOFFSET(status_data, simulation_has_to_stop        ),
                       HOFFSET(status_data, bond_dimension_has_reached_max),
                       HOFFSET(status_data, entanglement_has_converged    ),
                       HOFFSET(status_data, entanglement_has_saturated    ),
                       HOFFSET(status_data, variance_mpo_has_converged    ),
                       HOFFSET(status_data, variance_mpo_has_saturated    ),
                       HOFFSET(status_data, variance_ham_has_converged    ),
                       HOFFSET(status_data, variance_ham_has_saturated    ),
                       HOFFSET(status_data, variance_mom_has_converged    ),
                       HOFFSET(status_data, variance_mom_has_saturated    ),
                       HOFFSET(status_data, variance_mpo_saturated_for    ),
                       HOFFSET(status_data, variance_ham_saturated_for    ),
                       HOFFSET(status_data, variance_mom_saturated_for    )
                   };
           std::array       <size_t,NFIELDS>       dst_sizes   = {
                   sizeof(status_data::iteration                     ),
                   sizeof(status_data::step                          ),
                   sizeof(status_data::position                      ),
                   sizeof(status_data::chi_temp                      ),
                   sizeof(status_data::chi_max                       ),
                   sizeof(status_data::min_sweeps                    ),
                   sizeof(status_data::energy_min                    ),
                   sizeof(status_data::energy_max                    ),
                   sizeof(status_data::energy_target                 ),
                   sizeof(status_data::energy_ubound                 ),
                   sizeof(status_data::energy_lbound                 ),
                   sizeof(status_data::energy_dens                   ),
                   sizeof(status_data::energy_dens_target            ),
                   sizeof(status_data::energy_dens_window            ),
                   sizeof(status_data::phys_time                     ),
                   sizeof(status_data::wall_time                     ),
                   sizeof(status_data::simu_time                     ),
                   sizeof(status_data::delta_t                       ),
                   sizeof(status_data::time_step_has_converged       ),
                   sizeof(status_data::simulation_has_converged      ),
                   sizeof(status_data::simulation_has_saturated      ),
                   sizeof(status_data::simulation_has_to_stop        ),
                   sizeof(status_data::bond_dimension_has_reached_max),
                   sizeof(status_data::entanglement_has_converged    ),
                   sizeof(status_data::entanglement_has_saturated    ),
                   sizeof(status_data::variance_mpo_has_converged    ),
                   sizeof(status_data::variance_mpo_has_saturated    ),
                   sizeof(status_data::variance_ham_has_converged    ),
                   sizeof(status_data::variance_ham_has_saturated    ),
                   sizeof(status_data::variance_mom_has_converged    ),
                   sizeof(status_data::variance_mom_has_saturated    ),
                   sizeof(status_data::variance_mpo_saturated_for    ),
                   sizeof(status_data::variance_ham_saturated_for    ),
                   sizeof(status_data::variance_mom_saturated_for    )
           };
   
           std::array       <const char*,NFIELDS>  field_names =
                   {
                       "iteration",
                       "step",
                       "position",
                       "chi_temp",
                       "chi_max",
                       "min_sweeps",
                       "energy_min",
                       "energy_max",
                       "energy_target",
                       "energy_ubound",
                       "energy_lbound",
                       "energy_dens",
                       "energy_dens_target",
                       "energy_dens_window",
                       "phys_time",
                       "wall_time",
                       "simu_time",
                       "delta_t",
                       "time_step_has_converged",
                       "simulation_has_converged",
                       "simulation_has_saturated",
                       "simulation_has_to_stop",
                       "bond_dimension_has_reached_max",
                       "entanglement_has_converged",
                       "entanglement_has_saturated",
                       "variance_mpo_has_converged",
                       "variance_mpo_has_saturated",
                       "variance_ham_has_converged",
                       "variance_ham_has_saturated",
                       "variance_mom_has_converged",
                       "variance_mom_has_saturated",
                       "variance_mpo_saturated_for",
                       "variance_ham_saturated_for",
                       "variance_mom_saturated_for"
                   };
   
           std::array       <hid_t,NFIELDS>        field_types =
                   {
                           H5T_NATIVE_UINT,
                           H5T_NATIVE_UINT,
                           H5T_NATIVE_UINT,
                           H5T_NATIVE_LONG,
                           H5T_NATIVE_LONG,
                           H5T_NATIVE_UINT,
                           H5T_NATIVE_DOUBLE,
                           H5T_NATIVE_DOUBLE,
                           H5T_NATIVE_DOUBLE,
                           H5T_NATIVE_DOUBLE,
                           H5T_NATIVE_DOUBLE,
                           H5T_NATIVE_DOUBLE,
                           H5T_NATIVE_DOUBLE,
                           H5T_NATIVE_DOUBLE,
                           H5T_NATIVE_DOUBLE,
                           H5T_NATIVE_DOUBLE,
                           H5T_NATIVE_DOUBLE,
                           H5T_NATIVE_DOUBLE,
                           H5T_NATIVE_HBOOL,
                           H5T_NATIVE_HBOOL,
                           H5T_NATIVE_HBOOL,
                           H5T_NATIVE_HBOOL,
                           H5T_NATIVE_HBOOL,
                           H5T_NATIVE_HBOOL,
                           H5T_NATIVE_HBOOL,
                           H5T_NATIVE_HBOOL,
                           H5T_NATIVE_HBOOL,
                           H5T_NATIVE_HBOOL,
                           H5T_NATIVE_HBOOL,
                           H5T_NATIVE_HBOOL,
                           H5T_NATIVE_HBOOL,
                           H5T_NATIVE_UINT,
                           H5T_NATIVE_UINT,
                           H5T_NATIVE_UINT
                   };
   
           hsize_t          chunk_size                         = 1;
           void             *fill_data                         = nullptr;
           int              compress                           = 0;
       };
   public:
       class_log_simulation_status() = default;
       meta_struct meta;
       std::vector<status_data> buffer;
   };
   
   
   
   
   #endif //DMRG_LOG_TYPES_H
