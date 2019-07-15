
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_simulation_nmspc_settings.h:

Program Listing for File nmspc_settings.h
=========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_simulation_nmspc_settings.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/simulation/nmspc_settings.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 8/7/17.
   //
   
   #ifndef DMRG_N_SETTINGS_H
   #define DMRG_N_SETTINGS_H
   #include <string>
   #include <unordered_set>
   #include <vector>
   
   class class_settings_reader;
   namespace h5pp{
       class File;
   }
   
   enum class SimulationType      {iDMRG,fDMRG, xDMRG, iTEBD};
   enum class StorageLevel:size_t {NONE,LIGHT,NORMAL,FULL};
   
   namespace settings {
       extern void load_from_file(class_settings_reader &indata);
       extern void load_from_hdf5(h5pp::File &h5ppFile);
   
       namespace input{
           inline std::string input_file = "input/input.cfg";
           inline std::string input_filename = "input.cfg";
       }
       //Parameters for the model Hamiltonian
       namespace model {
           inline std::string  model_type     = "tf_ising";        
           inline int          seed_init      = 1;                 
           inline std::string  initial_sector = "sx";              
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
           inline size_t   eigMaxIter                   = 1000  ;   
           inline double   eigThreshold                 = 1e-12 ;   
           inline size_t   eigMaxNcv                    = 16    ;   
           inline double   SVDThreshold                 = 1e-8  ;   
           inline double   VarConvergenceThreshold      = 1e-8  ;   
           inline double   VarSaturationThreshold       = 1e-4  ;   
           inline double   EntEntrSaturationThreshold   = 1e-4  ;   
           inline size_t   MaxSizeFullDiag              = 2048  ;   
           inline size_t   MaxSizePartDiag              = 4096  ;   
       }
   
       //Parameters controlling iDMRG
       namespace idmrg {
           inline bool on           = true;                           
           inline size_t max_steps  = 5000;                           
           inline long chi_max      = 8;                              
           inline bool chi_grow     = true;                           
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
           inline size_t   print_freq   = 5000;                     
           inline size_t   write_freq   = 100;                      
       }
   
       namespace hdf5 {
           inline bool         save_logs            = true;                         
           inline bool         save_profiling       = true;                         
           inline std::string  access_mode          = "READWRITE" ;                 
           inline std::string  create_mode          = "RENAME";                     
           inline std::string  output_filename      = "output/default.h5";          
           inline StorageLevel storage_level        = StorageLevel::NORMAL;         
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
   #endif //DMRG_N_SETTINGS_H
