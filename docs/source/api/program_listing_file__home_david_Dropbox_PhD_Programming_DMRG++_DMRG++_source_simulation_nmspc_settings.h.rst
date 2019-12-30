
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_simulation_nmspc_settings.h:

Program Listing for File nmspc_settings.h
=========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_simulation_nmspc_settings.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/simulation/nmspc_settings.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 8/7/17.
   //
   
   #pragma once
   
   #include <string>
   #include <vector>
   #include <simulation/enums.h>
   class class_settings_reader;
   namespace h5pp{
       class File;
   }
   
   
   namespace settings {
       extern void load_from_file(class_settings_reader &indata);
       extern void load_from_hdf5(h5pp::File &h5ppFile);
   
       // Parameters for OpenMP
       // Make sure just one of these is > 1 otherwise too many threads
       // may be spawned inside of already threaded parts.
       namespace threading{
           inline int num_threads_eigen  = 1;                                                        
           inline int num_threads_omp    = 1;                                                        
           inline int num_threads_blas   = 1;                                                        
       }
   
       namespace input{
           inline std::string input_filename = "input/input.cfg";
           inline std::string input_file_raw;
       }
   
       namespace output {
           inline bool         save_logs            = true;                         
           inline bool         save_profiling       = true;                         
           inline std::string  access_mode          = "READWRITE" ;                 
           inline std::string  create_mode          = "RENAME";                     
           inline std::string  output_filename      = "output/default.h5";          
           inline StorageLevel storage_level        = StorageLevel::NORMAL;         
           inline bool         use_temp_dir         = true;                         
           inline size_t       copy_from_temp_freq  = 4;                            
           inline std::string  temp_dir             = "/scratch/local";             
       }
   
       //Parameters for the model Hamiltonian
       namespace model {
           inline std::string  model_type                              = "tf_ising";         
           inline int          seed_model                              = 1;                  
           inline int          seed_state                              = -1;                 
           inline bool         use_seed_state_as_enumeration           = true;               
           inline bool         projection_when_growing_chi             = true;               
           inline bool         projection_trial_when_stuck             = true;               
           inline bool         projection_on_every_sweep               = true;               
           inline bool         use_pauli_eigvecs                       = true;               
           inline std::string  initial_parity_sector                   = "x";                
           inline std::string  target_parity_sector                    = "x";                
           //Parameters for the transverse-field Ising model
           namespace tf_ising {
               inline double       J  = 1;                         
               inline double       g  = 1;                         
               inline double       w  = 0;                         
               inline size_t       d  = 2;                         
           }
   
           //Parameters for the transvese-field next-nearest neighbor Ising model
           namespace tf_nn_ising {
               inline double       J1  = 1;                         
               inline double       J2  = 1;                         
               inline double       g   = 1;                         
               inline double       w   = 0;                         
               inline size_t       d   = 2;                         
           }
   
           //Parameters for the selfdual transvese-field random-field next-neighbor Ising model
           namespace selfdual_tf_rf_ising {
               inline double       J_log_mean    = 0;               
               inline double       h_log_mean    = 0;               
               inline double       J_sigma       = 1;               
               inline double       h_sigma       = 0;               
               inline double       lambda        = 0;               
               inline size_t       d             = 2;               
           }
       }
   
       //Parmaters that control MPS, eigensolver and SVD precision
       namespace precision {
           inline size_t   eig_max_iter                    = 1000  ;   
           inline double   eig_threshold                   = 1e-12 ;   
           inline size_t   eig_max_ncv                     = 16    ;   
           inline double   svd_threshold                   = 1e-10 ;   
           inline double   variance_convergence_threshold  = 1e-11 ;   
           inline double   variance_slope_threshold        = 5     ;   
           inline double   entropy_slope_threshold         = 0.1   ;   
           inline double   subspace_error_factor           = 1     ;   
           inline double   max_subspace_error              = 1e-8  ;   
           inline double   min_subspace_error              = 1e-12 ;   
           inline size_t   max_sites_multidmrg             = 8     ;   
           inline size_t   max_size_full_diag              = 2048  ;   
           inline size_t   min_size_part_diag              = 4096  ;   
           inline size_t   max_size_direct                 = 131072;   
           inline double   max_norm_error                  = 1e-10 ;   
           inline size_t   max_resets                      = 4     ;   
           inline bool     use_reduced_energy              = true  ;   
           inline double   overlap_high                    = 0.99;
           inline double   overlap_cat                     = 0.70710678;
       }
   
       //Parameters controlling iDMRG
       namespace idmrg {
           inline bool on           = true;                           
           inline size_t max_steps  = 5000;                           
           inline long chi_max      = 8;                              
           inline bool chi_grow     = true;                           
           inline long chi_init     = 16;                             
           inline size_t print_freq = 1000;                           
           inline size_t write_freq = 100;                            
       }
       //Parameters controlling fDMRG
       namespace fdmrg {
           inline bool     on           = true;                         
           inline size_t   num_sites    = 16;                           
           inline size_t   max_sweeps   = 10;                           
           inline size_t   min_sweeps   = 4;                            
           inline long     chi_max      = 8;                            
           inline bool     chi_grow     = true;                         
           inline long     chi_init     = 16;                           
           inline size_t   print_freq   = 100;                          
           inline size_t   write_freq   = 100;                          
           inline bool     store_wavefn = false;                        
       }
   
       //Parameters controlling xDMRG
       namespace xdmrg {
           inline bool     on                      = true;             
           inline size_t   num_sites               = 16;               
           inline size_t   max_sweeps              = 10;               
           inline size_t   min_sweeps              = 4;                
           inline long     chi_max                 = 16;               
           inline bool     chi_grow                = true;             
           inline long     chi_init                = 16;               
           inline size_t   print_freq              = 1;                
           inline size_t   write_freq              = 1;                
           inline bool     store_wavefn            = false;            
           inline double   energy_density_target   = 0.5;              
           inline double   energy_density_window   = 0.05;             
       }
   
       //Parameters controlling iTEBD
       namespace itebd {
           inline bool     on           = true;                     
           inline size_t   max_steps    = 100000;                   
           inline double   delta_t0     = 0.1;                      
           inline double   delta_tmin   = 0.00001;                  
           inline size_t   suzuki_order = 1;                        
           inline long     chi_max      = 8;                        
           inline bool     chi_grow     = true;                     
           inline long     chi_init     = 16;                       
           inline size_t   print_freq   = 5000;                     
           inline size_t   write_freq   = 100;                      
       }
   
   
       //Profiling
       namespace profiling {
           inline bool     on        = false;             
           inline size_t   precision = 5;                 
       }
       //Console settings
       namespace console {
           inline size_t verbosity  = 2;                    
           inline bool   timestamp  = false;                
       }
   }
