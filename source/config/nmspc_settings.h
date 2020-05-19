//
// Created by david on 8/7/17.
//

#pragma once

#include <config/enums.h>
#include <string>
#include <vector>

/*!
 *  \namespace settings
 *  This namespace contains settings such as time-step length, number of iterations and precision parameters for
 *  the different algorithms.
 */

class class_config_reader;
class class_dmrg_config;
namespace h5pp {
    class File;
}

/* clang-format off */

namespace settings {
    extern void load_config(class_config_reader &indata);
    extern void load_config(class_dmrg_config &dmrg_config);
    extern void load_config(const std::string & config_filename);

    namespace threading{
        inline int num_threads = 1;                                              /*!< Number of threads for shared memory parallelism. num_threads <= 0 will try to use as many as possible */
    }

    namespace input{
        inline long        seed                                 = 1;                            /*!< Main seed for the random number generator. */
        inline long        bitfield                             = -1;                           /*!< Number whose bitfield represents the initial product state in the basis given by initial_parity_sector. Only positive state numbers are used */
        inline std::string config_filename                      = "input/input.cfg";            /*!< Default config filename. Can either be a .cfg file or a .h5 file with a config stored as a string in /common/config_file_contents */
        inline std::string config_file_contents;
    }

    namespace output {
        inline std::string         output_filepath                 = "output/default.h5";          /*!< Name of the output HDF5 file relative to the execution point  */
        inline bool                save_profiling                  = true;                         /*!< Whether to save profiling information to file */
        inline bool                checkpoint_keep_newest_only     = true;                         /*!< If true, checkpoints on each iteration will overwrite previous snapshots on file. Otherwise, all iterations are kept (dramaticallay increases file size) */
        inline bool                checkpoint_keep_chi_updates     = true;                         /*!< If true, a snapshot is written to file before updating the bond dimension is updated */
        inline size_t              checkpoint_frequency            = 1;                            /*!< How often, in units of iterations, to make a checkpoint. 0 disables checkpoints after iterations (chi-update checkpoints can still happen) */
        inline bool                use_temp_dir                    = true;                         /*!< If true uses a temporary directory for writes in the local drive (usually /tmp) and copies the results afterwards */
        inline size_t              copy_from_temp_freq             = 4;                            /*!< How often, in units of iterations, to copy the hdf5 file in tmp dir to target destination */
        inline std::string         temp_dir                        = "/tmp/DMRG";                  /*!< Local temp directory on the local system. If it does not exist we default to /tmp instead (or whatever is the default) */
        inline unsigned            compression_level               = 0;                            /*!< Attempt to use this compression level with HDF5. Choose between [0-9] (0 = off, 9 = max compression) */
        inline FileCollisionPolicy file_collision_policy           = FileCollisionPolicy::RESUME;  /*!< What to do when a prior output file is found. Choose between RESUME,RENAME,DELETE */

        // Storage Levels.
        // NOTE 1: A simulation can only be resumed from FULL state storage or checkpoint.
        // NOTE 2: storage_level_model == NORMAL is enough to recreate MPO's when resuming, since they can be reconstructed from the Hamiltonian parameter table
        //      NONE:   no data is saved at all
        //      LIGHT:  Mainly mid-chain data (energy/variance/polarization, schmidt values, entanglement entropy, lambda matrix, truncation error) , simulation status, and profiling (if save_profiling == true)
        //      NORMAL: Same as LIGHT + whole-chain measurements like entanglement entropies, truncation errors and schmidt values (lambda-matrices), and model Hamiltonian parameters
        //      FULL:   Same as NORMAL + MPS (Gamma + Lambda matrices) + MPO at each site.
        inline StorageLevel     storage_level_model      = StorageLevel::LIGHT;  /*!< Storage level for the model realization. LIGHT stores nothing. NORMAL stores the Hamiltonian parameter table, and FULL also the MPO's */
        inline StorageLevel     storage_level_checkpoint = StorageLevel::LIGHT;  /*!< Storage level for checkpoints, which are snapshots taken at each iteration or  snapshot taken at the end of each iteration */
        inline StorageLevel     storage_level_good_state = StorageLevel::NORMAL; /*!< Storage level for final results written when a simulation terminates successfully */
        inline StorageLevel     storage_level_fail_state = StorageLevel::NORMAL; /*!< Storage level for final results written when a simulation terminates unsuccessfully */
        inline StorageLevel     storage_level_proj_state = StorageLevel::LIGHT;  /*!< Storage level for the parity projected states, a projected version of the state written when a simulation terminates */
        inline StorageLevel     storage_level_init_state = StorageLevel::LIGHT;  /*!< Storage level for the initial states (for instance when launching a simulation or starting a new state) */
        inline StorageLevel     storage_level_emin_state = StorageLevel::LIGHT;  /*!< Storage level for the minimum energy state (ground state) */
        inline StorageLevel     storage_level_emax_state = StorageLevel::LIGHT;  /*!< Storage level for the maximum energy state */
    }


    //Profiling
    namespace profiling {
        inline bool     on        = false;                         /*!< If true, turns on profiling and timings will be shown on console. */
        inline size_t   precision = 5;                             /*!< Sets precision (number of decimals) of time output. */
    }
    //Console settings
    namespace console {
        inline size_t verbosity  = 2;                              /*!< Level of verbosity desired [0-6]. Level 0 prints everything, 6 nothing. Level 2 or 3 is recommended for normal use */
        inline bool   timestamp  = false;                          /*!< Whether to put a timestamp on console outputs */
    }

    #ifdef NDEBUG
        inline constexpr bool debug = false;
    #else
        inline constexpr bool debug = true;
    #endif



    //Parameters for the model Hamiltonian
    namespace model {
        inline ModelType    model_type = ModelType::ising_tf_rf;   /*!< Choice of model type: {ising_tf_rf_nn, ising_selfdual_tf_rf_nn}  */
        inline size_t       model_size = 16;                       /*!< Number of sites on the chain. Only relevant for finite algorithms: fDMRG and xDMRG */

        //Parameters for the transvese-field next-nearest neighbor Ising model with a random field
        namespace ising_tf_rf {
            inline double       J1         = 1;                 /*!< Ferromagnetic coupling for nearest neighbors.*/
            inline double       J2         = 1;                 /*!< Ferromagnetic coupling for next-nearest neighbors.*/
            inline double       h_tran     = 1;                 /*!< Transverse field strength */
            inline double       h_mean     = 0;                 /*!< Random field mean of distribution */
            inline double       h_stdv     = 0;                 /*!< Random field standard deviation. In distribution this is N(h_mean,h_stdv) or U(h_mean-h_stdv,h_mean+h_stdv) */
            inline long         spin_dim   = 2;                 /*!< Spin dimension */
            inline std::string  distribution  = "uniform";      /*!< Random distribution for couplings and fields */
        }

        //Parameters for the selfdual transverse-field random-field next-nearest neighbor Ising model
        namespace ising_sdual {
            inline double       J_mean        = 1;              /*!< Mean for the distribution defining random ferromagnetic coupling strength.*/
            inline double       h_mean        = 1;              /*!< Mean for the distribution defining random transverse magnetic field strength */
            inline double       J_stdv        = 1;              /*!< Standard deviation for the log-normal distribution defining ferromagnetic coupling */
            inline double       h_stdv        = 1;              /*!< Standard deviation for the log-normal distribution defining transverse magnetic field */
            inline double       lambda        = 0;              /*!< Lambda parameter related to next nearest neighbor coupling */
            inline bool         parity_sep    = false;          /*!< Separation of +-X parity sectors */
            inline long         spin_dim      = 2;              /*!< Spin dimension */
            inline std::string  distribution  = "lognormal";    /*!< Random distribution for couplings and fields */
        }
    }

    // Options that affect convergence
    namespace strategy {
        inline bool         chi_quench_when_stuck                   = false;              /*!< Reduce chi during a sweep when stuck and increasing bond dimension would not help */
        inline bool         perturb_when_stuck                      = false;              /*!< Perturb MPO parameters to get unstuck from local minima */
        inline bool         damping_when_stuck                      = false;              /*!< Modify MPO parameters, e.g. by reducing disorder, to get unstuck from local minima */
        inline bool         project_when_stuck                      = true;               /*!< Project to target parity sector at each sweep when stuck. */
        inline bool         project_on_every_sweep                  = true;               /*!< Project to target parity sector at each sweep. This implies doing it when stuck also. */
        inline bool         project_on_chi_update                   = true;               /*!< Project to target parity sector when bond dimension is increased (only works if chi_grow == true). */
        inline bool         randomize_on_chi_update                 = true;               /*!< Randomize MPS by flipping random spins when growing chi */
        inline bool         randomize_early                         = true;               /*!< Randomize MPS by flipping random spins before fully converging the first attempt (because the first attempt is biased) */
        inline bool         use_pauli_eigvecs                       = true;               /*!< Use random pauli eigenvectors to initialize_state spinors in x,y or z  */
        inline std::string  initial_parity_sector                   = "x";                /*!< Initialize in a global parity sector: {x,+x,-x, y, +y,-y, z,+z,-z, randomAxis,random,none}  */
        inline std::string  target_parity_sector                    = "x";                /*!< Project to in a global parity sector upon saturation: {x,+x,-x, y, +y,-y, z,+z,-z, randomAxis,random,none}  */
    }


    //Parmaters that control MPS, eigensolver and SVD precision
    namespace precision {
        inline size_t   eig_max_iter                    = 1000  ;   /*!< Maximum number of steps for eigenvalue solver. */
        inline double   eig_threshold                   = 1e-12 ;   /*!< Minimum threshold for halting eigenvalue solver. */
        inline size_t   eig_max_ncv                     = 16    ;   /*!< Parameter controlling the column space? of the Lanczos solver. */
        inline double   svd_threshold                   = 1e-10 ;   /*!< Minimum threshold value for keeping singular values. */
        inline double   variance_convergence_threshold  = 1e-11 ;   /*!< Desired precision on total energy variance. The MPS state is considered good enough when its energy variance reaches below this value */
        inline double   variance_slope_threshold        = 5     ;   /*!< Variance saturation slope threshold [0-100%]. The variance has saturated when its (absolute) slope reaches below this value. 2 would mean the data saturates when it changes less than 2% per iteration */
        inline double   entropy_slope_threshold         = 0.1   ;   /*!< Entanglement Entropy saturation slope threshold [0-100%]. The entanglement entropy has saturated when its (absolute) slope reaches below this value. 2 would mean the data saturates when it changes less than 2% per iteration*/
        inline double   subspace_error_factor           = 1     ;   /*!< The subspace quality threshold = energy_variance * SubspaceQualityFactor decides if we go ahead in variance optimization. If the subspace error is too high, direct optimization is done instead */
        inline double   max_subspace_error              = 1e-8  ;   /*!< The maximum subspace error. Never do subspace variance optimization with subspace error greater than this. */
        inline double   min_subspace_error              = 1e-12 ;   /*!< The minimum subspace error. Always do subspace variance optimization with subspace error less than this  */
        inline long     max_size_full_diag              = 2048  ;   /*!< Maximum linear size allowed for full diagonalization of the local hamiltonian matrix. */
        inline long     max_size_part_diag              = 4096  ;   /*!< Maximum linear size allowed for partial diagonalization of the local hamiltonian matrix. */
        inline long     max_size_direct                 = 131072;   /*!< Maximum linear size for direct multisite dmrg. If the linear size is larger than this, the algorithm prefers 2-site dmrg. */
        inline double   max_norm_error                  = 1e-10 ;   /*!< Maximum norm deviation from unity during integrity checks */
        inline size_t   max_resets                      = 4     ;   /*!< Maximum number of resets to initial state. One must be allowed for initialization */
        inline bool     use_reduced_energy              = true  ;   /*!< Whether to subtract E/L from each mpo to avoid catastrophic cancellation when computing the variance */
        inline double   overlap_high                    = 0.99;
        inline double   overlap_cat                     = 0.70710678;
        inline size_t   max_sites_multidmrg             = 8     ;    /*!< Maximum number of sites in multi-site dmrg. Too many sites (>10 or so) makes the contractions slow. */
        inline std::string move_multisite         = "one" ;    /*!< How many sites to move after a multi-site dmrg step, choose between {one,mid,max} */
    }

    //Parameters controlling iDMRG
    namespace idmrg {
        inline bool on           = true;                           /*!< Turns iDMRG simulation on/off. */
        inline size_t max_iters  = 5000;                           /*!< Maximum number of iDMRG iterations before forced termination */
        inline long chi_max      = 32;                             /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool chi_grow     = true;                           /*!< Whether to increase chi slowly up to chi_lim or go up to chi_lim directly. */
        inline long chi_init     = 16;                             /*!< Initial chi limit. Only used when chi_grow == true. */
        inline size_t print_freq = 1000;                           /*!< Print frequency for console output. (0 = off). */
//        inline size_t write_freq = 100;                            /*!< Write frequency,for output file buffer. (0 = off). */

    }


    //Parameters controlling iTEBD
    namespace itebd {
        inline bool     on           = true;                     /*!< Turns iTEBD simulation on/off. */
        inline size_t   max_iters    = 100000;                   /*!< Maximum number of iTEBD iterations before forced termination */
        inline double   delta_t0     = 0.1;                      /*!< Initial time step for iTEBD time evolution.*/
        inline double   delta_tmin   = 0.00001;                  /*!< Final time step for iTEBD time evolution.*/
        inline size_t   suzuki_order = 1;                        /*!< Order of the suzuki trotter decomposition (1,2 or 4) */
        inline long     chi_max      = 8;                        /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool     chi_grow     = true;                     /*!< Whether to increase chi slowly up to chi_lim or go up to chi_lim directly. */
        inline long     chi_init     = 16;                       /*!< Initial chi limit. Only used when chi_grow == true. */
        inline size_t   print_freq   = 5000;                     /*!< Print frequency for console output. (0 = off).*/
//        inline size_t   write_freq   = 100;                      /*!< Write frequency,for output file buffer. (0 = off). */

    }

    //Parameters controlling fDMRG
    namespace fdmrg {
        inline bool     on           = true;                         /*!< Turns fDMRG simulation on/off. */
        inline size_t   max_iters    = 10;                           /*!< Max number sweeps along the chain. */
        inline size_t   min_iters    = 4;                            /*!< Min number sweeps along the chain. */
        inline long     chi_max      = 8;                            /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool     chi_grow     = true;                         /*!< Whether to increase chi slowly up to chi_lim or go up to chi_lim directly. */
        inline long     chi_init     = 16;                           /*!< Initial chi limit. Only used when chi_grow == true. */
        inline size_t   print_freq   = 100;                          /*!< Print frequency for console output. In units of sweeps. (0 = off). */
//        inline size_t   write_freq   = 100;                          /*!< Write frequency,for output file buffer. In units of sweeps. (0 = off). */
        inline bool     store_wavefn = false;                        /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */
    }

    //Parameters controlling xDMRG
    namespace xdmrg {
        inline bool     on                      = true;             /*!< Turns xDMRG simulation on/off. */
        inline size_t   max_iters               = 10;               /*!< Max number sweeps along the chain. */
        inline size_t   min_iters               = 4;                /*!< Min number sweeps along the chain. */
        inline long     chi_max                 = 16;               /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool     chi_grow                = true;             /*!< Whether to increase chi slowly up to chi_lim or go up to chi_lim directly. */
        inline long     chi_init                = 16;               /*!< Initial chi limit. Only used when chi_grow == true. */
        inline size_t   print_freq              = 1;                /*!< Print frequency for console output. In units of sweeps. (0 = off). */
//        inline size_t   write_freq              = 1;                /*!< Write frequency,for output file buffer. In units of sweeps. (0 = off). */
        inline bool     store_wavefn            = false;            /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */
        inline double   energy_density_target   = 0.5;              /*!< Target energy in [0-1], where 0.5 means middle of spectrum. */
        inline double   energy_density_window   = 0.05;             /*!< Accept states inside of energy_target +- energy_window. */
        inline size_t   max_states              = 4;                /*!< Max number of random states to find using xDMRG on a single disorder realization */
    }

}
/* clang-format on */
