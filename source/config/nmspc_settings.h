//
// Created by david on 8/7/17.
//

#pragma once

#include "debug.h"
#include "enums.h"
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
        inline int omp_threads = 1;                                              /*!< Number of threads for openmp threads used in blas/lapack and Eigen. num_threads <= 0 will try to use as many as possible */
        inline int stl_threads = 1;                                              /*!< Number of threads for c++11 threading. Used in Eigen::Tensor. stl_threads <= 0 will try to use as many as possible */
    }

    namespace input{
        inline long        seed                                 = 1;                            /*!< Main seed for the random number generator. */
        inline long        bitfield                             = -1;                           /*!< Number whose bitfield represents the initial product state in the basis given by initial_parity_sector. Only positive state numbers are used */
        inline std::string config_filename                      = "input/input.config";            /*!< Default config filename. Can either be a .config file or a .h5 file with a config stored as a string in /common/config_file_contents */
        inline std::string config_file_contents;
    }

    namespace output {
        inline std::string         output_filepath                 = "output/default.h5";          /*!< Name of the output HDF5 file relative to the execution point  */
        inline bool                save_profiling                  = true;                         /*!< Whether to save profiling information to file */
        inline bool                savepoint_keep_newest_only      = true;                         /*!< If true, a savepoint will overwrite previous savepoints on file. Otherwise, all iterations are kept (dramaticallay increases file size) */
        inline size_t              savepoint_frequency             = 1;                            /*!< How often, in units of iterations, to make a savepoint. 0 disables regular savepoints but chi-update savepoints can still happen */
        inline bool                checkpoint_keep_newest_only     = true;                         /*!< If true, a checkpoint will overwrite previous checkpoint on file. Otherwise, all iterations are kept (dramaticallay increases file size) */
        inline bool                checkpoint_keep_chi_updates     = true;                         /*!< If true, a savepoint is written to file before the bond dimension is updated */
        inline size_t              checkpoint_frequency            = 1;                            /*!< How often, in units of iterations, to make a checkpoint. 0 disables checkpoints but chi-update checkpoints can still happen */
        inline bool                use_temp_dir                    = true;                         /*!< If true uses a temporary directory for writes in the local drive (usually /tmp) and copies the results afterwards */
        inline size_t              copy_from_temp_freq             = 4;                            /*!< How often, in units of iterations, to copy the hdf5 file in tmp dir to target destination */
        inline std::string         temp_dir                        = "/tmp/DMRG";                  /*!< Local temp directory on the local system. If it does not exist we default to /tmp instead (or whatever is the default) */
        inline unsigned            compression_level               = 1;                            /*!< GZip compression level in HDF5. Choose between [0-9] (0 = off, 9 = max compression) */
        inline FileCollisionPolicy file_collision_policy           = FileCollisionPolicy::RESUME;  /*!< What to do when a prior output file is found. Choose between RESUME,RENAME,DELETE */
        inline FileResumePolicy    file_resume_policy              = FileResumePolicy::FULL;       /*!< Depends on dataset "common/finished_all=bool" FULL: Ignore bool -> Scan .cfg to add missing items. FAST: exit if true. */

        // Storage Levels.
        // NOTE 1: A simulation can only be resumed from FULL state storage or savepoint.
        // NOTE 2: storage_level_model == NORMAL is enough to recreate MPO's when resuming, since they can be reconstructed from the Hamiltonian parameter table
        //      NONE:   no data is saved at all
        //      LIGHT:  Mainly mid-chain data (energy/variance/polarization, schmidt values, entanglement entropy, lambda matrix, truncation error) , simulation status, and profiling (if save_profiling == true)
        //      NORMAL: Same as LIGHT + whole-chain measurements like entanglement entropies, truncation errors and schmidt values (lambda-matrices), and model Hamiltonian parameters
        //      FULL:   Same as NORMAL + MPS (Gamma + Lambda matrices) + MPO at each site.
        inline StorageLevel     storage_level_model      = StorageLevel::LIGHT;  /*!< Storage level for the model realization. LIGHT stores nothing. NORMAL stores the Hamiltonian parameter table, and FULL also the MPO's */
        inline StorageLevel     storage_level_savepoint  = StorageLevel::LIGHT;  /*!< Storage level for savepoints, which are snapshots used for resume (if FULL) */
        inline StorageLevel     storage_level_checkpoint = StorageLevel::LIGHT;  /*!< Storage level for checkpoints, which are mid-simulation measurements (can also be used for resume if FULL) */
        inline StorageLevel     storage_level_good_state = StorageLevel::NORMAL; /*!< Storage level for final results written when a simulation terminates successfully */
        inline StorageLevel     storage_level_fail_state = StorageLevel::NORMAL; /*!< Storage level for final results written when a simulation terminates unsuccessfully */
        inline StorageLevel     storage_level_proj_state = StorageLevel::LIGHT;  /*!< Storage level for the parity projected states, a projected version of the state written when a simulation terminates */
        inline StorageLevel     storage_level_init_state = StorageLevel::LIGHT;  /*!< Storage level for the initial states (for instance when launching a simulation or starting a new state) */
        inline StorageLevel     storage_level_emin_state = StorageLevel::LIGHT;  /*!< Storage level for the minimum energy state (ground state) */
        inline StorageLevel     storage_level_emax_state = StorageLevel::LIGHT;  /*!< Storage level for the maximum energy state */

        namespace tmp{
            inline std::string hdf5_temp_path;
            inline std::string hdf5_final_path;
        }
    }


    //Profiling
    namespace profiling {
        inline bool     on        = false;                         /*!< If true, turns on profiling and timings will be shown on console. */
        inline bool     extra     = false;                         /*!< Prints more profiling */
        inline size_t   precision = 5;                             /*!< Sets precision (number of decimals) of time output. */
    }
    //Console settings
    namespace console {
        inline size_t verbosity  = 2;                              /*!< Level of verbosity desired [0-6]. Level 0 prints everything, 6 nothing. Level 2 or 3 is recommended for normal use */
        inline bool   timestamp  = false;                          /*!< Whether to put a timestamp on console outputs */
    }

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
            inline double       lambda        = 0;              /*!< Lambda parameter related to next nearest neighbor coupling */
            inline double       delta         = 0;              /*!< Delta defined as log(J_mean) - log(h_mean). We get J_mean and h_mean by fixing max(J_mean,h_mean) = 1 */
            inline double       J_stdv        = 1;              /*!< Standard deviation for the log-normal distribution defining ferromagnetic coupling */
            inline double       h_stdv        = 1;              /*!< Standard deviation for the log-normal distribution defining transverse magnetic field */
            inline bool         parity_sep    = false;          /*!< Separation of +-X parity sectors */
            inline long         spin_dim      = 2;              /*!< Spin dimension */
            inline std::string  distribution  = "lognormal";    /*!< Random distribution for couplings and fields */
        }

        //Parameters for the l-bit hamiltonian
        namespace lbit {
            inline double       J1_mean       = 0;              /*!< Constant offset for on-site */
            inline double       J2_mean       = 0;              /*!< Constant offset for two-body interaction */
            inline double       J3_mean       = 0;              /*!< Constant offset for three-body interaction */
            inline double       J1_wdth       = 0.5;            /*!< Width of the uniform box distribution U(-w1,w1) for on-site interactions */
            inline double       J2_wdth       = 0.5;            /*!< Width of the uniform box distribution U(-J2_wdth,J2_wdth) for two-body interaction */
            inline double       J3_wdth       = 0.5;            /*!< Width of the uniform box distribution U(-J3_wdth,J3_wdth) for three-body interaction */
            inline double       f_mixer       = 0.1;            /*!< Mixing factor for unitary transformation to real-space */
            inline size_t       u_layer       = 6;              /*!< Number of unitary 2-site layers which transform lbit <-> real spaces */
            inline long         spin_dim      = 2;              /*!< Spin dimension */
            inline std::string  distribution  = "uniform";      /*!< Random distribution for interaction strengths */
        }

    }

    // Options for strategy that affect convergence and targeted state
    namespace strategy {
        inline bool          compress_mpo_squared       = true;                                   /*!< Use SVD to compress the squared mpo bond dimension */
        inline bool          chi_quench_when_stuck      = false;                                  /*!< Reduce chi for a few iterations when stuck and increasing bond dimension would not help */
        inline bool          perturb_when_stuck         = false;                                  /*!< Perturb MPO parameters to get unstuck from local minima */
        inline bool          damping_when_stuck         = false;                                  /*!< Modify MPO parameters, e.g. by reducing disorder, to get unstuck from local minima */
        inline bool          discard_schmidt_when_stuck = true;                                   /*!< Try discarding smallest schmidt values when stuck */
        inline bool          expand_subspace_when_stuck = true;                                   /*!< Use subspace expansion when stuck in local minima */
        inline size_t        project_when_stuck_freq    = 4;                                      /*!< Project to target parity sector every nth iteration when stuck. (0 = turn off) */
        inline bool          project_on_every_iter      = true;                                   /*!< Project to target parity sector at the end of every iteration. This implies doing it when stuck also. */
        inline bool          project_on_chi_update      = true;                                   /*!< Project to target parity sector when bond dimension is increased (only works if cfg_chi_lim_grow == true). */
        inline bool          project_initial_state      = false;                                  /*!< Project to target parity sector when initializing a state. */
        inline bool          randomize_on_chi_update    = true;                                   /*!< Randomize MPS by flipping random spins when growing chi */
        inline bool          randomize_early            = true;                                   /*!< Randomize MPS by flipping random spins before fully converging the first attempt (because the first attempt is biased) */
        inline bool          use_eigenspinors           = false;                                  /*!< Use random pauli-matrix eigenvectors when initializing each mps site along x,y or z  */
        inline size_t        max_resets                 = 1;                                      /*!< Maximum number of resets to product state due to saturation. One must be allowed for initialization */
        inline size_t        multisite_max_sites        = 8;                                      /*!< Maximum number of sites in multi-site dmrg. Too many sites (>10 or so) makes the contractions slow. */
        inline MultisiteMove multisite_move             = MultisiteMove::ONE;                     /*!< How many sites to move after a multi-site dmrg step, choose between {ONE, MID, MAX} */
        inline StateInitType initial_type               = StateInitType::REAL;                    /*!< Initial state can be REAL/CPLX */
        inline StateInit     initial_state              = StateInit::RANDOM_ENTANGLED_STATE;      /*!< Initial configuration for the spin chain (only for finite systems)  */
        inline StateInit     secondary_states           = StateInit::RANDOMIZE_PREVIOUS_STATE;    /*!< Spin configuration for subsequent states (only for finite systems)  */
        inline std::string   target_sector              = "x";                                    /*!< Find an eigenstate in this parity sector. Choose between {x,+x,-x, y, +y,-y, z,+z,-z, randomAxis,random,none}  */
}


    //Parmaters that control precision in MPS, eigensolver and SVD
    namespace precision {
        inline size_t   eig_max_iter                    = 1000  ;   /*!< Maximum number of steps for eigenvalue solver. */
        inline double   eig_threshold                   = 1e-12 ;   /*!< Minimum threshold for halting eigenvalue solver. */
        inline size_t   eig_max_ncv                     = 16    ;   /*!< Parameter controlling the column space? of the Lanczos solver. */
        inline double   svd_threshold                   = 1e-10 ;   /*!< Minimum threshold value for keeping singular values. */
        inline size_t   svd_switchsize                  = 16    ;   /*!< Linear size of a matrix, below which BDCSVD will use slower but more precise JacobiSVD instead (default is 16) */
        inline double   variance_convergence_threshold  = 1e-11 ;   /*!< Desired precision on total energy variance. The MPS state is considered good enough when its energy variance reaches below this value */
        inline double   variance_slope_threshold        = 5     ;   /*!< Variance saturation slope threshold [0-100%]. The variance has saturated when its (absolute) slope reaches below this value. 2 would mean the data saturates when it changes less than 2% per iteration */
        inline double   entropy_slope_threshold         = 0.1   ;   /*!< Entanglement Entropy saturation slope threshold [0-100%]. The entanglement entropy has saturated when its (absolute) slope reaches below this value. 2 would mean the data saturates when it changes less than 2% per iteration*/
        inline double   subspace_error_factor           = 1     ;   /*!< The subspace quality threshold = energy_variance * SubspaceQualityFactor decides if we go ahead in variance optimization. If the subspace error is too high, direct optimization is done instead */
        inline double   max_subspace_error              = 1e-8  ;   /*!< The maximum subspace error. Never do subspace variance optimization with subspace error greater than this. */
        inline double   min_subspace_error              = 1e-12 ;   /*!< The minimum subspace error. Always do subspace variance optimization with subspace error less than this  */
        inline size_t   max_subspace_size               = 256   ;   /*!< Maximum number of candidate eigenstates to keep for a subspace-optimization step */
        inline long     max_size_full_diag              = 2048  ;   /*!< Maximum linear size allowed for full diagonalization of the local hamiltonian matrix. */
        inline long     max_size_part_diag              = 4096  ;   /*!< Maximum linear size allowed for partial diagonalization of the local hamiltonian matrix. */
        inline long     max_size_direct                 = 131072;   /*!< Maximum linear size for direct multisite dmrg. If the linear size is larger than this, the algorithm prefers 2-site dmrg. */
        inline double   max_norm_error                  = 1e-10 ;   /*!< Maximum norm deviation from unity during integrity checks */
        inline bool     use_reduced_energy              = true  ;   /*!< Whether to subtract E/L from each mpo to avoid catastrophic cancellation when computing the variance */
        inline double   overlap_high                    = 0.99;
        inline double   overlap_cat                     = 0.70710678;
    }


    //Parameters controlling iDMRG
    namespace idmrg {
        inline bool on           = true;                           /*!< Turns iDMRG simulation on/off. */
        inline size_t max_iters  = 5000;                           /*!< Maximum number of iDMRG iterations before forced termination */
        inline long chi_lim_max  = 32;                             /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool chi_lim_grow = true;                           /*!< Whether to increase chi slowly up to chi_lim or go up to chi_lim directly. */
        inline long chi_lim_init = 16;                             /*!< Initial chi limit. Only used when cfg_chi_lim_grow == true. */
        inline size_t print_freq = 1000;                           /*!< Print frequency for console output. In units of iterations.  (0 = off). */
    }


    //Parameters controlling iTEBD
    namespace itebd {
        inline bool     on                    = true;                /*!< Turns iTEBD simulation on/off. */
        inline size_t   max_iters             = 100000;              /*!< Maximum number of iTEBD iterations before forced termination */
        inline double   time_step_init_real   = 0.0;                 /*!< Real part of initial time step delta_t */
        inline double   time_step_init_imag   = 0.1;                 /*!< Imag part of initial time step delta_t */
        inline double   time_step_min         = 0.00001;             /*!< (Absolute value) Minimum and final time step for iTEBD time evolution. */
        inline size_t   suzuki_order          = 1;                   /*!< Order of the suzuki trotter decomposition (1,2 or 4) */
        inline long     chi_lim_max           = 8;                   /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool     chi_lim_grow          = true;                /*!< Whether to increase chi slowly up to chi_lim or go up to chi_lim directly. */
        inline long     chi_lim_init          = 16;                  /*!< Initial chi limit. Only used when cfg_chi_lim_grow == true. */
        inline size_t   print_freq            = 5000;                /*!< Print frequency for console output. In units of iterations. (0 = off).*/
    }

    //Parameters controlling fdmrg
    namespace fdmrg {
        inline bool     on           = true;                         /*!< Turns fdmrg simulation on/off. */
        inline size_t   max_iters    = 10;                           /*!< Max number of iterations. One iterations moves L steps. */
        inline size_t   min_iters    = 4;                            /*!< Min number of iterations. One iterations moves L steps. */
        inline long     chi_lim_max  = 8;                            /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool     chi_lim_grow = true;                         /*!< Whether to increase chi slowly up to chi_lim or go up to chi_lim directly. */
        inline long     chi_lim_init = 16;                           /*!< Initial chi limit. Only used when cfg_chi_lim_grow == true. */
        inline size_t   print_freq   = 100;                          /*!< Print frequency for console output. In units of iterations. (0 = off). */
        inline bool     store_wavefn = false;                        /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */
    }

    //Parameters controlling flbit
    namespace flbit {
        inline bool     on                      = false;                    /*!< Turns flbit simulation on/off. */
        inline size_t   max_iters               = 10000;                    /*!< Max number of iterations. One iterations moves L steps. */
        inline size_t   min_iters               = 4;                        /*!< Min number of iterations. One iterations moves L steps. */
        inline long     chi_lim_max             = 8;                        /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool     chi_lim_grow            = true;                     /*!< Whether to increase chi slowly up to chi_lim or go up to chi_lim directly. */
        inline long     chi_lim_init            = 16;                       /*!< Initial chi limit. Only used when cfg_chi_lim_grow == true. */
        inline double   time_start_real         = 1e-1;                     /*!< Starting time point (real) */
        inline double   time_start_imag         = 0;                        /*!< Starting time point (imag) */
        inline double   time_final_real         = 1e6;                      /*!< Finishing time point (real) */
        inline double   time_final_imag         = 0;                        /*!< Finishing time point (imag) */
        inline size_t   time_num_steps          = 500;                      /*!< Number of steps from start to finish */
        inline size_t   print_freq              = 1;                        /*!< Print frequency for console output. In units of iterations. (0 = off). */
        inline bool     store_wavefn            = false;                    /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */
    }

    //Parameters controlling xDMRG
    namespace xdmrg {
        inline bool     on                              = true;             /*!< Turns xDMRG simulation on/off. */
        inline size_t   max_iters                       = 10;               /*!< Max number of iterations. One iterations moves L steps. */
        inline size_t   min_iters                       = 4;                /*!< Min number of iterations. One iterations moves L steps. */
        inline size_t   olap_iters                      = 2;                /*!< Number of initial iterations selecting the candidate state with best overlap to the current state */
        inline size_t   vsub_iters                      = 2;                /*!< Number of iterations using the subspace optimization for variance, after overlap iterations */
        inline long     chi_lim_max                     = 16;               /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool     chi_lim_grow                    = true;             /*!< Whether to increase chi slowly up to chi_lim or go up to chi_lim directly. */
        inline long     chi_lim_init                    = 16;               /*!< Initial chi limit. Only used when chi_grow == true, or starting from an entangled state. */
        inline long     chi_lim_olap                    = 16;               /*!< Chi limit during initial OVERLAP|SUBSPACE mode. set to <= 0 for unlimited */
        inline long     chi_lim_vsub                    = 32;               /*!< Chi limit during initial VARIANCE|SUBSPACE mode. set to <= 0 for unlimited */
        inline size_t   print_freq                      = 1;                /*!< Print frequency for console output. In units of iterations. (0 = off). */
        inline double   energy_density_target           = 0.5;              /*!< Target energy in [0-1], where 0.5 means middle of spectrum. */
        inline double   energy_density_window           = 0.05;             /*!< Accept states inside of energy_tgt_per_site +- energy_dens_window. */
        inline size_t   max_states                      = 1;                /*!< Max number of random states to find using xDMRG on a single disorder realization */
        inline bool     store_wavefn                    = false;            /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */
        inline bool     finish_if_entanglm_saturated    = true;             /*!< Finish early as soon as entanglement has saturated */
        inline bool     finish_if_variance_saturated    = false;            /*!< Finish early as soon as energy variance has saturated */

    }
}
/* clang-format on */
