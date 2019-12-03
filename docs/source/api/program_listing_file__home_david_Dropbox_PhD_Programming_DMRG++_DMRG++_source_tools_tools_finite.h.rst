
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_tools_finite.h:

Program Listing for File tools_finite.h
=======================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_tools_finite.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/tools_finite.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-10-13.
   //
   
   #pragma once
   
   #include <memory>
   #include <string>
   #include <vector>
   #include <complex>
   
   class class_state_finite;
   class class_model_base;
   class class_simulation_status;
   
   
   namespace tools{
       namespace finite
       {
           namespace mps {
               extern void initialize                          (class_state_finite & state, size_t length);
               extern void randomize                           (class_state_finite &state, const std::string &parity_sector = "random", int seed_state = -1, bool use_pauli_eigenstates = false, bool enumeration =  false);
               extern void normalize                           (class_state_finite & state, const std::vector<size_t> & bond_dimensions = {});
               extern void rebuild_environments                (class_state_finite & state);
               extern int  move_center_point                   (class_state_finite & state);          
               extern void project_to_closest_parity_sector    (class_state_finite & state, std::string paulistring, bool keep_bond_dimensions = false);
   
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
           }
   
           namespace ops {
               extern std::list<Eigen::Tensor<Scalar,4>>
               make_mpo_list                 (const std::list<std::unique_ptr<class_model_base>> & mpos_L, const std::list<std::unique_ptr<class_model_base>> & mpos_R);
               extern void apply_mpo                     (class_state_finite & state, const Eigen::Tensor<Scalar,4> & mpo, const Eigen::Tensor<Scalar,3> &Ledge, const Eigen::Tensor<Scalar,3> & Redge);
               extern void apply_mpos                    (class_state_finite & state, const std::list<Eigen::Tensor < Scalar,4>> & mpos, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge);
               extern class_state_finite
               get_projection_to_parity_sector           (const class_state_finite & state, const Eigen::MatrixXcd & paulimatrix, int sign, bool keep_bond_dimensions = false);
               extern class_state_finite
               get_projection_to_closest_parity_sector   (const class_state_finite & state, const Eigen::MatrixXcd & paulimatrix, bool keep_bond_dimensions = false);
               extern class_state_finite
               get_projection_to_closest_parity_sector   (const class_state_finite & state, std::string parity_sector, bool keep_bond_dimensions = false);
               extern double overlap                     (const class_state_finite & state1, const class_state_finite & state2);
               extern double expectation_value           (const class_state_finite & state1, const class_state_finite & state2, const std::list<Eigen::Tensor < Scalar,4>> & mpos, const Eigen::Tensor<Scalar,3> & Ledge, const Eigen::Tensor<Scalar,3> & Redge);
               extern double exp_sq_value                (const class_state_finite & state1, const class_state_finite & state2, const std::list<Eigen::Tensor < Scalar,4>> & mpos, const Eigen::Tensor<Scalar,4> & Ledge, const Eigen::Tensor<Scalar,4> & Redge);
           }
   
   
           namespace opt{
   //            enum class OptMode  {OVERLAP, VARIANCE};
   //            enum class OptSpace {SUBSPACE,DIRECT};
   //            enum class OptType  {REAL, CPLX};
   
   
               extern Eigen::Tensor<Scalar,3> find_excited_state(const class_state_finite & state, const class_simulation_status & sim_status, OptMode optMode, OptSpace optSpace, OptType optType);
               extern Eigen::Tensor<Scalar,4> find_ground_state (const class_state_finite & state, std::string ritz = "SR");
               extern void truncate_theta(Eigen::Tensor<Scalar,3> & theta, class_state_finite & state, long chi_, double SVDThreshold);
               extern void truncate_left (Eigen::Tensor<Scalar,3> & theta, class_state_finite & state, long chi_, double SVDThreshold);
               extern void truncate_right(Eigen::Tensor<Scalar,3> & theta, class_state_finite & state, long chi_, double SVDThreshold);
               extern void truncate_theta(Eigen::Tensor<Scalar,4> & theta, class_state_finite & state, long chi_, double SVDThreshold);
           }
   
           namespace multisite{
               extern Eigen::DSizes<long,3> get_dimensions  (const class_state_finite &state, const std::list<size_t> &list_of_sites);
               extern size_t                get_problem_size(const class_state_finite &state, const std::list<size_t> &list_of_sites);
               extern std::list<size_t>     generate_site_list(class_state_finite &state, const size_t threshold, const size_t max_sites);
           }
   
   
   
   
   
           namespace print {
               extern void print_full_state    (const class_state_finite & state);
               extern void print_state         (const class_state_finite & state);                                                
               extern void print_state_compact (const class_state_finite & state);                                                
               extern void print_hamiltonians  (const class_state_finite & state);
           }
   
           namespace io{
               namespace internals{
                   inline bool make_extendable_dataset(const std::string & prefix_path);
               }
               extern void write_all_state                              (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
               extern void write_bond_matrices                          (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
               extern void write_bond_matrix                            (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
               extern void write_full_mps                               (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
               extern void write_full_mpo                               (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
               extern void write_model                                  (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
               extern void write_entanglement                           (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
               extern void write_all_measurements                       (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path);
               extern void write_projection_to_closest_parity_sector    (const class_state_finite & state, h5pp::File & h5ppFile, const std::string & prefix_path, std::string parity_sector, bool keep_bond_dimensions);
               extern void load_from_hdf5                               (const h5pp::File & h5ppFile, class_state_finite & state    , class_simulation_status & sim_status, const std::string & prefix_path);
               extern class_state_finite load_state_from_hdf5           (const h5pp::File & h5ppFile, const std::string & prefix_path);
           }
   
   
           namespace debug {
               extern void check_integrity             (const class_state_finite & state);
               extern void check_integrity_of_mps      (const class_state_finite & state);
               extern void check_integrity_of_mpo      (const class_state_finite & state);
               extern void check_normalization_routine (const class_state_finite & state);
               extern void print_parity_properties     (const class_state_finite & state);
   
           }
   
       }
   }
