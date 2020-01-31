//
// Created by david on 8/7/17.
//

#pragma once

#include <string>
#include <vector>
#include <simulation/enums.h>
/*!
 *  \namespace settings
 *  This namespace contains settings such as time-step length, number of iterations and precision parameters for
 *  the different algorithms.
 */

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
        inline int num_threads_eigen  = 1;                                                        /*!< Number of threads for Eigen operations. num_threads <= 0 will try to use as many as possible */
        inline int num_threads_omp    = 1;                                                        /*!< Number of threads for OpenMP operations. num_threads <= 0 will try to use as many as possible */
        inline int num_threads_blas   = 1;                                                        /*!< Number of threads for BLAS operations. num_threads <= 0 will try to use as many as possible */
    }

    namespace input{
        inline std::string input_filename = "input/input.cfg";
        inline std::string input_file_raw;
    }

    namespace output {
        inline bool         save_logs            = true;                         /*!< If true, saves the history of the simulation in log files, not just the end results  (only enabled on storage level NORMAL and FULL.) */
        inline bool         save_profiling       = true;                         /*!< Whether to save profiling information to file. (only enabled on storage level NORMAL and FULL.) */
        inline std::string  access_mode          = "READWRITE" ;                 /*!< Choose access mode to the file. Choose between READWRITE, READONLY */
        inline std::string  create_mode          = "RENAME";                     /*!< Choose access mode to the file. Choose between TRUNCATE, OPEN, RENAME */
        inline std::string  output_filename      = "output/default.h5";          /*!< Name of the output HDF5 file relative to the execution point  */
        inline StorageLevel storage_level        = StorageLevel::NORMAL;         /*!< Sets the storage level: choose "0=NONE,1=LIGHT,2=NORMAL,3=FULL */
        inline bool         use_temp_dir         = true;                         /*!< If true uses a temporary directory for writes in the local drive (usually /tmp) and copies the results afterwards */
        inline size_t       copy_from_temp_freq  = 4;                            /*!< How often, in units of iterations, to copy the hdf5 file in tmp dir to target destination */
        inline std::string  temp_dir             = "/scratch/local";             /*!< Local temp directory on the local system. If it does not exist we default to /tmp instead (or whatever is the default) */
    }

    //Parameters for the model Hamiltonian
    namespace model {
        inline std::string  model_type                              = "tf_ising";         /*!< Choice of model type: {tf_ising, tf_nn_ising, selfdual_tf_rf_ising selfdual_tf_rf_ising_normal}  */
        inline long         seed                                    = 1;                  /*!< Main seed for the random number generator. */
        inline long         state_number                            = -1;                 /*!< Number whose bitfield represents the initial product state in the basis given by initial_parity_sector. Only positive state numbers are used */
        inline bool         quench_chi_when_stuck                   = false;              /*!< Reduce chi during a sweep when stuck and increasing bond dimension would not help */
        inline bool         projection_when_growing_chi             = true;               /*!< Project to target parity sector when bond dimension is increased (only works if chi_grow == true). */
        inline bool         projection_trial_when_stuck             = true;               /*!< Project to target parity sector at each sweep when stuck. */
        inline bool         projection_on_every_sweep               = true;               /*!< Project to target parity sector at each sweep. This implies doing it when stuck also. */
        inline bool         use_pauli_eigvecs                       = true;               /*!< Use random pauli eigenvectors to initialize spinors in x,y or z  */
        inline std::string  initial_parity_sector                   = "x";                /*!< Initialize in a global parity sector: {x,+x,-x, y, +y,-y, z,+z,-z, randomAxis,random,none}  */
        inline std::string  target_parity_sector                    = "x";                /*!< Project to in a global parity sector upon saturation: {x,+x,-x, y, +y,-y, z,+z,-z, randomAxis,random,none}  */

        //Parameters for the transverse-field Ising model
        namespace tf_ising {
            inline double       J  = 1;                         /*!< Ferromagnetic coupling. J < 0  Gives a ferromagnet. J > 0 an antiferromagnet. */
            inline double       g  = 1;                         /*!< Transverse field strength */
            inline double       w  = 0;                         /*!< Randomness strength for the random field */
            inline size_t       d  = 2;                         /*!< Spin dimension */
        }

        //Parameters for the transvese-field next-nearest neighbor Ising model
        namespace tf_nn_ising {
            inline double       J1  = 1;                         /*!< Ferromagnetic coupling for nearest neighbors.*/
            inline double       J2  = 1;                         /*!< Ferromagnetic coupling for next-nearest neighbors.*/
            inline double       g   = 1;                         /*!< Transverse field strength */
            inline double       w   = 0;                         /*!< Randomness strength for the random field */
            inline size_t       d   = 2;                         /*!< Spin dimension */
        }

        //Parameters for the selfdual transverse-field random-field next-neighbor Ising model
        namespace selfdual_tf_rf_ising {
            inline double       J_mean        = 1;               /*!< Mean for the distribution defining random ferromagnetic coupling strength.*/
            inline double       h_mean        = 1;               /*!< Mean for the distribution defining random transverse magnetic field strength */
            inline double       J_sigma       = 1;               /*!< Standard deviation for the log-normal distribution defining ferromagnetic coupling */
            inline double       h_sigma       = 1;               /*!< Standard deviation for the log-normal distribution defining transverse magnetic field */
            inline double       lambda        = 0;               /*!< Lambda parameter related to next nearest neighbor coupling */
            inline bool         parity_sep    = false;           /*!< Separation of +-X parity sectors */
            inline std::string  distribution  = "lognormal";     /*!< Random distribution for couplings and fields */
            inline size_t       d             = 2;               /*!< Spin dimension */
        }
    }

    //Parmaters that control MPS, eigensolver and SVD precision
    namespace precision {
        inline size_t   eig_max_iter                    = 1000  ;   /*!< Maximum number of steps for eigenvalue solver. */
        inline double   eig_threshold                   = 1e-12 ;   /*!< Minimum threshold for halting eigenvalue solver. */
        inline size_t   eig_max_ncv                     = 16    ;   /*!< Parameter controlling the column space? of the Lanczos solver. */
        inline double   svd_threshold                   = 1e-10 ;   /*!< Minimum threshold value for keeping singular values. */
        inline double   variance_convergence_threshold  = 1e-11 ;   /*!< Variance convergence threshold. The MPS state is considered good enough when its variance reaches below this value */
        inline double   variance_slope_threshold        = 5     ;   /*!< Variance saturation slope threshold [0-100%]. The variance has saturated when its (absolute) slope reaches below this value. 2 would mean the data saturates when it changes less than 2% per iteration */
        inline double   entropy_slope_threshold         = 0.1   ;   /*!< Entanglement Entropy saturation slope threshold [0-100%]. The entanglement entropy has saturated when its (absolute) slope reaches below this value. 2 would mean the data saturates when it changes less than 2% per iteration*/
        inline double   subspace_error_factor           = 1     ;   /*!< The subspace quality threshold = energy_variance * SubspaceQualityFactor decides if we go ahead in variance optimization. If the subspace error is too high, direct optimization is done instead */
        inline double   max_subspace_error              = 1e-8  ;   /*!< The maximum subspace error. Never do subspace variance optimization with subspace error greater than this. */
        inline double   min_subspace_error              = 1e-12 ;   /*!< The minimum subspace error. Always do subspace variance optimization with subspace error less than this  */
        inline size_t   max_size_full_diag              = 2048  ;   /*!< Maximum linear size allowed for full diagonalization of the local hamiltonian matrix. */
        inline size_t   max_size_part_diag              = 4096  ;   /*!< Maximum linear size allowed for partial diagonalization of the local hamiltonian matrix. */
        inline size_t   max_size_direct                 = 131072;   /*!< Maximum linear size for direct multisite dmrg. If the linear size is larger than this, the algorithm prefers 2-site dmrg. */
        inline double   max_norm_error                  = 1e-10 ;   /*!< Maximum norm deviation from unity during integrity checks */
        inline size_t   max_resets                      = 4     ;   /*!< Maximum number of resets to initial state. One must be allowed for initialization */
        inline bool     use_reduced_energy              = true  ;   /*!< Whether to subtract E/L from each mpo to avoid catastrophic cancellation when computing the variance */
        inline double   overlap_high                    = 0.99;
        inline double   overlap_cat                     = 0.70710678;
        inline size_t   max_sites_multidmrg             = 8     ;   /*!< Maximum number of sites in multi-site dmrg. Too many sites (>10 or so) makes the contractions slow. */
        inline std::string move_sites_multidmrg         = "one" ;    /*!< How many sites to move after a multi-site dmrg step, choose between {one,mid,max} */

    }

    //Parameters controlling iDMRG
    namespace idmrg {
        inline bool on           = true;                           /*!< Turns iDMRG simulation on/off. */
        inline size_t max_steps  = 5000;                           /*!< Final length of 1D quantum chain. */
        inline long chi_max      = 8;                              /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool chi_grow     = true;                           /*!< Whether to increase chi slowly up to chi_lim or go up to chi_lim directly. */
        inline long chi_init     = 16;                             /*!< Initial chi limit. Only used when chi_grow == true. */
        inline size_t print_freq = 1000;                           /*!< Print frequency for console output. (0 = off). */
        inline size_t write_freq = 100;                            /*!< Write frequency,for output file buffer. (0 = off). */

    }
    //Parameters controlling fDMRG
    namespace fdmrg {
        inline bool     on           = true;                         /*!< Turns fDMRG simulation on/off. */
        inline size_t   num_sites    = 16;                           /*!< Number of sites on the chain */
        inline size_t   max_sweeps   = 10;                           /*!< Max number sweeps along the chain. */
        inline size_t   min_sweeps   = 4;                            /*!< Min number sweeps along the chain. */
        inline long     chi_max      = 8;                            /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool     chi_grow     = true;                         /*!< Whether to increase chi slowly up to chi_lim or go up to chi_lim directly. */
        inline long     chi_init     = 16;                           /*!< Initial chi limit. Only used when chi_grow == true. */
        inline size_t   print_freq   = 100;                          /*!< Print frequency for console output. In units of sweeps. (0 = off). */
        inline size_t   write_freq   = 100;                          /*!< Write frequency,for output file buffer. In units of sweeps. (0 = off). */
        inline bool     store_wavefn = false;                        /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */
    }

    //Parameters controlling xDMRG
    namespace xdmrg {
        inline bool     on                      = true;             /*!< Turns xDMRG simulation on/off. */
        inline size_t   num_sites               = 16;               /*!< Number of sites on the chain */
        inline size_t   max_sweeps              = 10;               /*!< Max number sweeps along the chain. */
        inline size_t   min_sweeps              = 4;                /*!< Min number sweeps along the chain. */
        inline long     chi_max                 = 16;               /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool     chi_grow                = true;             /*!< Whether to increase chi slowly up to chi_lim or go up to chi_lim directly. */
        inline long     chi_init                = 16;               /*!< Initial chi limit. Only used when chi_grow == true. */
        inline size_t   print_freq              = 1;                /*!< Print frequency for console output. In units of sweeps. (0 = off). */
        inline size_t   write_freq              = 1;                /*!< Write frequency,for output file buffer. In units of sweeps. (0 = off). */
        inline bool     store_wavefn            = false;            /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */
        inline double   energy_density_target   = 0.5;              /*!< Target energy in [0-1], where 0.5 means middle of spectrum. */
        inline double   energy_density_window   = 0.05;             /*!< Accept states inside of energy_target +- energy_window. */
    }

    //Parameters controlling iTEBD
    namespace itebd {
        inline bool     on           = true;                     /*!< Turns iTEBD simulation on/off. */
        inline size_t   max_steps    = 100000;                   /*!< Number of iTEBD iterations, after which the simulation terminates regardless of convergence. Set high.*/
        inline double   delta_t0     = 0.1;                      /*!< Initial time step for iTEBD time evolution.*/
        inline double   delta_tmin   = 0.00001;                  /*!< Final time step for iTEBD time evolution.*/
        inline size_t   suzuki_order = 1;                        /*!< Order of the suzuki trotter decomposition (1,2 or 4) */
        inline long     chi_max      = 8;                        /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool     chi_grow     = true;                     /*!< Whether to increase chi slowly up to chi_lim or go up to chi_lim directly. */
        inline long     chi_init     = 16;                       /*!< Initial chi limit. Only used when chi_grow == true. */
        inline size_t   print_freq   = 5000;                     /*!< Print frequency for console output. (0 = off).*/
        inline size_t   write_freq   = 100;                      /*!< Write frequency,for output file buffer. (0 = off). */

    }


    //Profiling
    namespace profiling {
        inline bool     on        = false;             /*!< If true, turns on profiling and timings will be shown on console. */
        inline size_t   precision = 5;                 /*!< Sets precision (number of decimals) of time output. */
    }
    //Console settings
    namespace console {
        inline size_t verbosity  = 2;                    /*!< Level of verbosity desired [0-6]. Level 0 prints everything, 6 nothing. Level 2 or 3 is recommended for normal use */
        inline bool   timestamp  = false;                /*!< Whether to put a timestamp on console outputs */
    }
}
