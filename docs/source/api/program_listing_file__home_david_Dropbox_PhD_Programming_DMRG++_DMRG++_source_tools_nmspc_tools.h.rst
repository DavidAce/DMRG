
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_nmspc_tools.h:

Program Listing for File nmspc_tools.h
======================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_nmspc_tools.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/nmspc_tools.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-01-29.
   //
   #pragma once
   
   #include <memory>
   #include <string>
   #include <complex>
   #include <list>
   #include <vector>
   #include <io/nmspc_logger.h>
   #include <general/class_tic_toc.h>
   #include <tools/finite/opt-internals/enum_classes.h>
   #include <Eigen/Core>
   #include <unsupported/Eigen/CXX11/Tensor>
   class class_state_infinite;
   class class_mps_2site;
   class class_state_finite;
   class class_model_base;
   class class_simulation_status;
   
   
   namespace h5pp{
       class File;
   }
   template <typename table_type> class class_h5table_buffer;
   class class_h5table_measurements_finite;
   class class_h5table_measurements_infinite;
   class class_h5table_profiling;
   class class_h5table_simulation_status;
   
   
   
   
   namespace tools{
       inline std::shared_ptr<spdlog::logger> log;
       namespace finite
       {
           using Scalar = std::complex<double>;
           namespace mps {
               extern void initialize                          (class_state_finite & state, size_t length);
               extern void randomize                           (class_state_finite &state, const std::string &parity_sector = "random", int seed_state = -1, bool use_pauli_eigenstates = false, bool enumeration =  false);
               extern void normalize                           (class_state_finite & state, const std::optional<size_t> chi_lim = std::nullopt);
               extern void rebuild_environments                (class_state_finite & state);
               extern int  move_center_point                   (class_state_finite & state);          
               extern void project_to_closest_parity_sector    (class_state_finite & state, std::string paulistring);
   
               namespace internals{
                   inline bool seed_state_unused = true;
                   extern void set_product_state_in_parity_sector_from_bitset(class_state_finite & state, const std::string &parity_sector, const int seed_state);
                   extern void set_product_state_in_parity_sector_randomly(class_state_finite & state, const std::string &parity_sector);
                   extern void set_product_state_randomly(class_state_finite & state, const std::string &parity_sector, bool use_pauli_eigenstates);
               }
   
           }
   
           namespace mpo {
               extern void initialize                 (class_state_finite & state, size_t length, std::string model_type);
               extern void randomize                  (class_state_finite & state, int seed_state = -1);
               extern void reduce_mpo_energy          (class_state_finite & state);
               extern void reduce_mpo_energy_multi    (class_state_finite & state);
               extern void reduce_mpo_energy_2site    (class_state_finite & state);
           }
   
           namespace ops {
               extern std::list<Eigen::Tensor<Scalar,4>>
                           make_mpo_list                 (const std::list<std::unique_ptr<class_model_base>> & mpos_L, const std::list<std::unique_ptr<class_model_base>> & mpos_R);
               extern void apply_mpo                     (class_state_finite & state, const Eigen::Tensor<Scalar,4> & mpo, const Eigen::Tensor<Scalar,3> &Ledge, const Eigen::Tensor<Scalar,3> & Redge);
               extern void apply_mpos                    (class_state_finite & state, const std::list<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge);
               extern class_state_finite
               get_projection_to_parity_sector           (const class_state_finite & state, const Eigen::MatrixXcd & paulimatrix, int sign);
               extern class_state_finite
               get_projection_to_closest_parity_sector   (const class_state_finite & state, const Eigen::MatrixXcd & paulimatrix);
               extern class_state_finite
               get_projection_to_closest_parity_sector   (const class_state_finite & state, std::string parity_sector);
               extern double overlap                     (const class_state_finite & state1, const class_state_finite & state2);
               extern double expectation_value           (const class_state_finite & state1, const class_state_finite & state2, const std::list<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge);
               extern double exp_sq_value                (const class_state_finite & state1, const class_state_finite & state2, const std::list<Eigen::Tensor<Scalar,4>> & mpos, const Eigen::Tensor<Scalar,4> & Ledge, const Eigen::Tensor<Scalar,4> & Redge);
           }
   
   
           namespace opt{
   //            enum class OptMode  {OVERLAP, VARIANCE};
   //            enum class OptSpace {SUBSPACE,DIRECT};
   //            enum class OptType  {REAL, CPLX};
   
   
               extern Eigen::Tensor<Scalar,3> find_excited_state(const class_state_finite & state, const class_simulation_status & sim_status, OptMode optMode, OptSpace optSpace, OptType optType);
               extern Eigen::Tensor<Scalar,4> find_ground_state (const class_state_finite & state, std::string ritz = "SR");
               extern void truncate_theta(Eigen::Tensor<Scalar,3> & theta, class_state_finite & state);
               extern void truncate_left (Eigen::Tensor<Scalar,3> & theta, class_state_finite & state);
               extern void truncate_right(Eigen::Tensor<Scalar,3> & theta, class_state_finite & state);
               extern void truncate_theta(Eigen::Tensor<Scalar,4> & theta, class_state_finite & state,const std::optional<size_t> chi_lim = std::nullopt);
           }
   
           namespace multisite{
               extern Eigen::DSizes<long,3> get_dimensions  (const class_state_finite &state, const std::list<size_t> &list_of_sites);
               extern size_t                get_problem_size(const class_state_finite &state, const std::list<size_t> &list_of_sites);
               extern std::list<size_t>     generate_site_list(class_state_finite &state, const size_t threshold, const size_t max_sites, const size_t min_sites = 2);
           }
   
   
           namespace measure{
   
   //            extern void do_all_measurements                           (class_state_finite & state);
               extern int length                                         (const class_state_finite & state);
               extern size_t bond_dimension_current                      (const class_state_finite & state);
               extern size_t bond_dimension_midchain                     (const class_state_finite & state);
               extern std::vector<size_t> bond_dimensions                (const class_state_finite & state);
               extern double norm                                        (const class_state_finite & state);
   
   
               extern double spin_component                              (const class_state_finite & state, const Eigen::Matrix2cd &paulimatrix);
               extern Eigen::Tensor<Scalar,1> mps_wavefn                 (const class_state_finite & state);
               extern double entanglement_entropy_current                (const class_state_finite & state);
               extern double entanglement_entropy_midchain               (const class_state_finite & state);
               extern std::vector<double> entanglement_entropies         (const class_state_finite & state);
               extern std::vector<double> spin_components                (const class_state_finite & state);
   
               namespace twosite{
                   extern double energy_minus_energy_reduced                 (const class_state_finite & state, const Eigen::Tensor<Scalar,4> & theta);
                   extern double energy                                      (const class_state_finite & state, const Eigen::Tensor<Scalar,4> & theta);
                   extern double energy_per_site                             (const class_state_finite & state, const Eigen::Tensor<Scalar,4> & theta);
                   extern double energy_variance                             (const class_state_finite & state, const Eigen::Tensor<Scalar,4> & theta);
                   extern double energy_variance_per_site                    (const class_state_finite & state, const Eigen::Tensor<Scalar,4> & theta);
               }
   
               namespace multisite{
                   namespace internal{
                       inline double digits;
                       double significant_digits(double H2, double E2);
                   }
                   extern double energy_minus_energy_reduced             (const class_state_finite & state, const Eigen::Tensor<Scalar,3> & multitheta);
                   extern double energy                                  (const class_state_finite & state, const Eigen::Tensor<Scalar,3> & multitheta);
                   extern double energy_per_site                         (const class_state_finite & state, const Eigen::Tensor<Scalar,3> & multitheta);
                   extern double energy_variance                         (const class_state_finite & state, const Eigen::Tensor<Scalar,3> & multitheta);
                   extern double energy_variance_per_site                (const class_state_finite & state, const Eigen::Tensor<Scalar,3> & multitheta);
                   extern double energy                                  (const class_state_finite & state);
                   extern double energy_per_site                         (const class_state_finite & state);
                   extern double energy_variance                         (const class_state_finite & state);
                   extern double energy_variance_per_site                (const class_state_finite & state);
               }
   
               extern double energy                                      (const class_state_finite & state);
               extern double energy_per_site                             (const class_state_finite & state);
               extern double energy_variance                             (const class_state_finite & state);
               extern double energy_variance_per_site                    (const class_state_finite & state);
               extern double energy_normalized                           (const class_state_finite & state, const class_simulation_status & sim_status);
               template<typename Derived>
               double energy_minus_energy_reduced(const class_state_finite & state, const Eigen::TensorBase<Derived,Eigen::ReadOnlyAccessors> & theta){
                   constexpr int rank = Derived::NumIndices;
                   if constexpr (rank == 4) return twosite::energy_minus_energy_reduced(state,theta);
                   if constexpr (rank == 3) return multisite::energy_minus_energy_reduced(state,theta);
                   static_assert("Wrong rank, expected 3 or 4" and (rank == 3 or rank == 4));
               }
               template<typename Derived>
               double energy(const class_state_finite & state, const Eigen::TensorBase<Derived,Eigen::ReadOnlyAccessors> & theta){
                   constexpr int rank = Derived::NumIndices;
                   if constexpr (rank == 4) return twosite::energy(state,theta);
                   if constexpr (rank == 3) return multisite::energy(state,theta);
                   static_assert("Wrong rank, expected 3 or 4" and (rank == 3 or rank == 4));
               }
               template<typename Derived>
               double energy_per_site(const class_state_finite & state, const Eigen::TensorBase<Derived,Eigen::ReadOnlyAccessors> & theta){
                   constexpr int rank = Derived::NumIndices;
                   if constexpr (rank == 4) return twosite::energy_per_site(state,theta);
                   if constexpr (rank == 3) return multisite::energy_per_site(state,theta);
                   static_assert("Wrong rank, expected 3 or 4" and (rank == 3 or rank == 4));
               }
               template<typename Derived>
               double energy_variance(const class_state_finite & state, const Eigen::TensorBase<Derived,Eigen::ReadOnlyAccessors> & theta){
                   constexpr int rank = Derived::NumIndices;
                   if constexpr (rank == 4) return twosite::energy_variance(state,theta);
                   if constexpr (rank == 3) return multisite::energy_variance(state,theta);
                   static_assert("Wrong rank, expected 3 or 4" and (rank == 3 or rank == 4));
               }
               template<typename Derived>
               double energy_variance_per_site(const class_state_finite & state, const Eigen::TensorBase<Derived,Eigen::ReadOnlyAccessors> & theta){
                   constexpr int rank = Derived::NumIndices;
                   if constexpr (rank == 4) return twosite::energy_variance_per_site(state,theta);
                   if constexpr (rank == 3) return multisite::energy_variance_per_site(state,theta);
                   static_assert("Wrong rank, expected 3 or 4" and (rank == 3 or rank == 4));
   
               }
           }
   
   
           namespace print {
               extern void print_full_state    (const class_state_finite & state);
               extern void print_state         (const class_state_finite & state);                                                
               extern void print_state_compact (const class_state_finite & state);                                                
               extern void print_hamiltonians  (const class_state_finite & state);
           }
   
           namespace io{
   
               namespace h5dset{
                   namespace internals{
                       inline bool make_extendable_dataset(const std::string & prefix_path);
                   }
                   extern void write_all_state                              (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
                   extern void write_bond_matrices                          (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
                   extern void write_bond_matrix                            (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
                   extern void write_full_mps                               (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
                   extern void write_full_mpo                               (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
                   extern void write_model                                  (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
                   extern void write_array_measurements                     (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
   
               }
   
   
               namespace h5table{
                   extern void write_measurements                       (const class_state_finite &state, const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_measurements_finite> &h5tbuf);
                   extern void write_sim_status                         (const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_simulation_status> &h5tbuf);
                   extern void write_profiling                          (const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_profiling> &h5tbuf);
               }
   
               namespace h5restore{
                   extern void load_from_hdf5                               (const h5pp::File & h5ppFile, class_state_finite & state    , class_simulation_status & sim_status, const std::string & prefix_path);
                   extern class_state_finite load_state_from_hdf5           (const h5pp::File & h5ppFile, const std::string & prefix_path);
               }
   
   
           }
   
   
           namespace debug {
               extern void check_integrity             (const class_state_finite & state);
               extern void check_integrity_of_mps      (const class_state_finite & state);
               extern void check_integrity_of_mpo      (const class_state_finite & state);
               extern void check_normalization_routine (const class_state_finite & state);
               extern void print_parity_properties     (const class_state_finite & state);
   
           }
   
       }
   
   
   
   
       namespace infinite
       {
           using Scalar = std::complex<double>;
   
           namespace mps{
               extern void initialize                          (class_state_infinite & state, std::string model_type_str);
               extern class_state_infinite set_random_state    (const class_state_infinite & state, [[maybe_unused]]  std::string parity, [[maybe_unused]] int seed_state);
           }
   
           namespace env{
               extern void initialize                          (class_state_infinite & state);
   //            extern void rebuild_environments                (class_state_infinite & state);
           }
   
           namespace mpo{
               extern void initialize                          (class_state_infinite & state, std::string model_type_str);
               extern void randomize                           (class_state_infinite &state, int seed_model);
   
               }
           namespace opt{
               extern Eigen::Tensor<Scalar,4> find_ground_state(const class_state_infinite & state, std::string ritz = "SR");
               extern Eigen::Tensor<Scalar,4> time_evolve_theta(const class_state_infinite & state, const Eigen::Tensor<Scalar, 4> &U);
               extern void truncate_theta(Eigen::Tensor<Scalar,4> &theta, class_state_infinite & state);
   
           }
   
           namespace measure{
               extern int    length                          (const class_state_infinite & state);
               extern int    bond_dimension                  (const class_state_infinite & state);
               extern double truncation_error                (const class_state_infinite & state);
               extern double norm                            (const class_state_infinite & state);
               extern double energy_mpo                      (const class_state_infinite & state);
               extern double energy_mpo                      (const class_state_infinite & state, const Eigen::Tensor<Scalar,4> &theta);
               extern double energy_per_site_mpo             (const class_state_infinite & state);
               extern double energy_per_site_ham             (const class_state_infinite & state);
               extern double energy_per_site_mom             (const class_state_infinite & state);
               extern double energy_variance_mpo             (const class_state_infinite & state, const Eigen::Tensor<Scalar,4> &theta, double &energy_mpo);
               extern double energy_variance_mpo             (const class_state_infinite & state, const Eigen::Tensor<Scalar,4> &theta);
               extern double energy_variance_mpo             (const class_state_infinite & state);
               extern double energy_variance_per_site_mpo    (const class_state_infinite & state);
               extern double energy_variance_per_site_ham    (const class_state_infinite & state);
               extern double energy_variance_per_site_mom    (const class_state_infinite & state);
               extern double entanglement_entropy    (const class_state_infinite & state);
           }
   
           namespace print {
               extern void print_state         (const class_state_infinite & state);                                                
               extern void print_state_compact (const class_state_infinite & state);                                                
               extern void print_hamiltonians  (const class_state_infinite & state);
           }
   
           namespace io{
               namespace h5dset{
                   extern void write_all_state(const class_state_infinite &state, h5pp::File &h5ppFile, std::string sim_name);
                   extern void write_2site_mps                    (const class_state_infinite & state, h5pp::File & h5ppFile, std::string sim_name);
                   extern void write_2site_mpo                    (const class_state_infinite & state, h5pp::File & h5ppFile, std::string sim_name);
                   extern void write_2site_env                    (const class_state_infinite & state, h5pp::File & h5ppFile, std::string sim_name);
                   extern void write_2site_env2                   (const class_state_infinite & state, h5pp::File & h5ppFile, std::string sim_name);
                   extern void write_hamiltonian_params           (const class_state_infinite & state, h5pp::File & h5ppFile, std::string sim_name);
                   extern void write_all_measurements             (const class_state_infinite & state, h5pp::File & h5ppFile, std::string sim_name);
               }
   
               namespace h5table{
                   extern void write_measurements                       (const class_state_infinite &state, const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_measurements_infinite> &h5tbuf);
                   extern void write_sim_status                         (const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_simulation_status> &h5tbuf);
                   extern void write_profiling                          (const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_profiling> &h5tbuf);
               }
   
               namespace h5restore{
                   extern void load_from_hdf5                     (const h5pp::File & h5ppFile, class_state_infinite & state, class_simulation_status &sim_status, std::string sim_name);
                   extern void load_superblock_from_hdf5          (const h5pp::File & h5ppFile, class_state_infinite & state, std::string sim_name);
                   extern void load_sim_status_from_hdf5          (const h5pp::File & h5ppFile, class_simulation_status & sim_status, std::string sim_name);
               }
   
   
           }
   
   
   
           namespace debug {
               extern void check_integrity             (const class_state_infinite & state);
               extern void check_integrity_of_mps      (const class_state_infinite & state);
               extern void check_normalization_routine (const class_state_infinite & state);
   
           }
   
   
       }
   
   
   
   
   
       namespace common{
           using Scalar = std::complex<double>;
   
           namespace io {
               namespace h5dset{
                   extern void write_simulation_status(const class_simulation_status &sim_status, h5pp::File &h5ppFile, std::string sim_name);
               }
               namespace h5table{
                   extern void write_sim_status                         (const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_simulation_status> &h5tbuf);
                   extern void write_profiling                          (const class_simulation_status &sim_status, class_h5table_buffer<class_h5table_profiling> &h5tbuf);
               }
   
               namespace h5restore{
                   extern class_simulation_status load_sim_status_from_hdf5(const h5pp::File &h5ppFile, std::string sim_name);
               }
   
               namespace h5tmp{
                   extern std::string set_tmp_prefix(const std::string &output_filename);
                   extern std::string unset_tmp_prefix(const std::string &output_filename);
                   extern void copy_from_tmp(const std::string & output_filename);
                   extern void create_directory(const std::string & dir);
                   extern void remove_from_temp(const std::string output_filename);
   
               }
   
           }
   
   
           namespace profile{
               inline class_tic_toc t_eig;
               inline class_tic_toc t_svd;
               inline class_tic_toc t_ene;
               inline class_tic_toc t_var;
               inline class_tic_toc t_ent;
               inline class_tic_toc t_hdf;
               inline class_tic_toc t_prj;
               inline class_tic_toc t_opt;
               inline class_tic_toc t_chk;
               inline class_tic_toc t_ene_mpo;
               inline class_tic_toc t_ene_ham;
               inline class_tic_toc t_ene_mom;
               inline class_tic_toc t_var_mpo;
               inline class_tic_toc t_var_ham;
               inline class_tic_toc t_var_mom;
   
               extern void print_profiling(class_tic_toc &t_parent);
               extern void init_profiling();
           }
   
   
           namespace views {
               extern Eigen::Tensor<Scalar,4> theta, theta_evn_normalized, theta_odd_normalized;
               extern Eigen::Tensor<Scalar,4> theta_sw ;
               extern Eigen::Tensor<Scalar,3> LAGA, LCGB;
               extern Eigen::Tensor<Scalar,2> l_evn, r_evn;
               extern Eigen::Tensor<Scalar,2> l_odd, r_odd;
               extern Eigen::Tensor<Scalar,4> transfer_matrix_LAGA;
               extern Eigen::Tensor<Scalar,4> transfer_matrix_LCGB;
               extern Eigen::Tensor<Scalar,4> transfer_matrix_evn;
               extern Eigen::Tensor<Scalar,4> transfer_matrix_odd;
               extern bool components_computed;
               extern void compute_mps_components(const class_state_infinite &state);
   
               extern Eigen::Tensor<Scalar,4> get_theta                       (const class_state_finite & state, std::complex<double> norm = 1.0);              
               extern Eigen::Tensor<Scalar,4> get_theta                       (const class_state_infinite & state, std::complex<double> norm = 1.0);              
               extern Eigen::Tensor<Scalar,4> get_theta_swapped               (const class_state_infinite & state, std::complex<double> norm = 1.0);              
               extern Eigen::Tensor<Scalar,4> get_theta_evn                   (const class_state_infinite & state, std::complex<double> norm = 1.0);              
               extern Eigen::Tensor<Scalar,4> get_theta_odd                   (const class_state_infinite & state, std::complex<double> norm = 1.0);              
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_zero        (const class_state_infinite & state);
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_LBGA        (const class_state_infinite & state, std::complex<double> norm = 1.0);
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_GALC        (const class_state_infinite & state, std::complex<double> norm = 1.0);
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_GBLB        (const class_state_infinite & state, std::complex<double> norm = 1.0);
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_LCGB        (const class_state_infinite & state, std::complex<double> norm = 1.0);
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_theta_evn   (const class_state_infinite & state, std::complex<double> norm = 1.0);
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_theta_odd   (const class_state_infinite & state, std::complex<double> norm = 1.0);
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_AB          (const class_state_infinite & state, int p);
   
               extern Eigen::Tensor<Scalar,4> get_theta                       (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);              
               extern Eigen::Tensor<Scalar,4> get_theta_swapped               (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);              
               extern Eigen::Tensor<Scalar,4> get_theta_evn                   (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);              
               extern Eigen::Tensor<Scalar,4> get_theta_odd                   (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);              
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_zero        (const class_mps_2site  &MPS);
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_LBGA        (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_GALC        (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_GBLB        (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_LCGB        (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_theta_evn   (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_theta_odd   (const class_mps_2site  &MPS, std::complex<double> norm = 1.0);
               extern Eigen::Tensor<Scalar,4> get_transfer_matrix_AB          (const class_mps_2site  &MPS, int p);
   
   
           }
       }
   
   }
   
   
   
   
   
   
