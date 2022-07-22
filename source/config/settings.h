#pragma once

#include "debug.h"
#include "enums.h"
#include "tid/enums.h"
#include <string>
#include <vector>

class Loader;
namespace h5pp {
    class File;
}

/* clang-format off */


/*!
 *  \namespace settings
 *  This namespace contains settings such as time-step length, number of iterations and precision parameters for
 *  the different algorithms.
 */
namespace settings {
    extern void load(Loader &dmrg_config);
    extern void load(std::string_view  config_filename);

    extern bool     algorithm_is_on(AlgorithmType algo_type);
    extern size_t   print_freq(AlgorithmType algo_type);
    extern long     get_bond_init(AlgorithmType algo_type);
    extern long     get_bond_max(AlgorithmType algo_type);
    extern bool     store_wave_function(AlgorithmType algo_type);





    /*!  \namespace settings::threading Parameters for multithreading */
    namespace threading{
        inline int omp_threads = 1;                                              /*!< Number of threads for openmp threads used in blas/lapack and Eigen. num_threads <= 0 will try to use as many as possible */
        inline int stl_threads = 1;                                              /*!< Number of threads for c++11 threading. Used in Eigen::Tensor. stl_threads <= 0 will try to use as many as possible */
    }

    /*!  \namespace settings::input Settings for initialization */
    namespace input{
        inline long        seed                                 = 1;                            /*!< Main seed for the random number generator. */
        inline size_t      bitfield                             = -1ul;                         /*!< Number whose bitfield represents the initial product state in the basis given by initial_parity_sector. Disable with -1ul */
        inline std::string config_filename                      = "input/input.cfg";            /*!< Default config filename. Can either be a .cfg file or a .h5 file with a config stored as a string in /common/config_file_contents */
        inline std::string config_file_contents;
    }

    /*!  \namespace settings::storage Settings for output-file generation
     *
     * **NOTE: Storage levels**
     *
     *  enum StorageLevel:
     *       - `NONE`:   no data is saved at all
     *       - `LIGHT`:  Mainly mid-chain data (energy/variance/polarization, schmidt values, entanglement entropy, lambda matrix, truncation error) , simulation status, and timers (if save_timers == true)
     *       - `NORMAL`: Same as `LIGHT` + whole-chain measurements like entanglement entropies, truncation errors and schmidt values (lambda-matrices), and model Hamiltonian parameters
     *       - `FULL`:   Same as `NORMAL` + MPS (Gamma + Lambda matrices) + MPO at each site.
     *
     * **Note: Resume**
     *
     * Simulations can only be resumed from a state or checkpoint saved with `StorageLevel::FULL`.
     *
     * The only exception is `storage_level_model == StorageLevel::NORMAL` which is enough to recreate MPOs, since they can be reconstructed from the Hamiltonian parameter table
     *
     */
    namespace storage {
        inline std::string         output_filepath                 = "output/output.h5";           /*!< Name of the output HDF5 file relative to the execution point  */
        inline bool                save_timers                     = true;                         /*!< Save timer information to file (on storage level FULL|NORMAL) */
        inline bool                savepoint_keep_newest_only      = true;                         /*!< If true, a savepoint will overwrite previous savepoints on file. Otherwise, all iterations are kept (dramaticallay increases file size) */
        inline size_t              savepoint_frequency             = 1;                            /*!< How often, in units of iterations, to make a savepoint. 0 disables regular savepoints but bond-update savepoints can still happen */
        inline bool                checkpoint_keep_newest_only     = true;                         /*!< If true, a checkpoint will overwrite previous checkpoint on file. Otherwise, all iterations are kept (dramaticallay increases file size) */
        inline size_t              checkpoint_frequency            = 1;                            /*!< How often, in units of iterations, to make a checkpoint. 0 disables checkpoints but bond-update checkpoints can still happen */
        inline bool                use_temp_dir                    = true;                         /*!< If true uses a temporary directory for writes in the local drive (usually /tmp) and copies the results afterwards */
        inline size_t              copy_from_temp_freq             = 4;                            /*!< How often, in units of iterations, to copy the hdf5 file in tmp dir to target destination */
        inline std::string         temp_dir                        = "/tmp/DMRG";                  /*!< Local temp directory on the local system. If it does not exist we default to /tmp instead (or whatever is the default) */
        inline unsigned            compression_level               = 1;                            /*!< GZip compression level in HDF5. Choose between [0-9] (0 = off, 9 = max compression) */
        inline FileCollisionPolicy file_collision_policy           = FileCollisionPolicy::RESUME;  /*!< What to do when a prior output file is found. Choose between RESUME, REVIVE, BACKUP, RENAME, REPLACE */
        inline FileResumePolicy    file_resume_policy              = FileResumePolicy::FULL;       /*!< Depends on dataset "common/finished_all=bool" FULL: Ignore bool -> Scan .cfg to add missing items. FAST: exit if true. */
        inline std::string         file_resume_name                = ""  ;                         /*!< On file_collision_policy=RESUME|REVIVE: resume from state candidate matching this string. Empty implies any */
        inline size_t              file_resume_iter                = -1ul;                         /*!< On file_collision_policy=RESUME|REVIVE: which iteration to resume from. -1ul implies resume from last available iteration */

        inline StorageLevel     storage_level_model      = StorageLevel::LIGHT;  /*!< Storage level for the model realization. LIGHT stores nothing. NORMAL stores the Hamiltonian parameter table, and FULL also the MPO's */
        inline StorageLevel     storage_level_savepoint  = StorageLevel::LIGHT;  /*!< Storage level for savepoints, which are snapshots used for resume */
        inline StorageLevel     storage_level_checkpoint = StorageLevel::LIGHT;  /*!< Storage level for checkpoints, which are mid-simulation measurements */
        inline StorageLevel     storage_level_finished   = StorageLevel::NORMAL; /*!< Storage level for final results written when a simulation terminates */
        inline StorageLevel     storage_level_proj_state = StorageLevel::LIGHT;  /*!< Storage level for the parity projected states, a projected version of the state written when a simulation terminates */
        inline StorageLevel     storage_level_init_state = StorageLevel::LIGHT;  /*!< Storage level for the initial states (for instance when launching a simulation or starting a new state) */
        inline StorageLevel     storage_level_emin_state = StorageLevel::LIGHT;  /*!< Storage level for the minimum energy state (ground state) */
        inline StorageLevel     storage_level_emax_state = StorageLevel::LIGHT;  /*!< Storage level for the maximum energy state */
        inline StorageLevel     storage_level_bond_state = StorageLevel::NORMAL; /*!< Storage level for states written on bond limit change */
        inline StorageLevel     storage_level_trnc_state = StorageLevel::NORMAL; /*!< Storage level for states written on truncation error limit change */
        inline StorageLevel     storage_level_fes_state  = StorageLevel::NORMAL; /*!< Storage level for states written during finite entanglement scaling (fes) after the main simulation */

        namespace tmp{
            inline std::string hdf5_temp_path;
            inline std::string hdf5_final_path;
        }
    }


    /*!  \namespace settings::timer Settings for performance profiling */
    namespace timer {
        inline bool         on        = false;                         /*!< If true, turns on timers. These will be shown on console. */
        inline tid::level   level     = tid::normal;                   /*!< How much extra to print on exit [normal | extra | detailed]  */
        inline size_t       precision = 5;                             /*!< Sets precision (number of decimals) of time storage. */
    }

    /*! \namespace settings::console Settings for console output */
    namespace console {
//        inline size_t verbosity     = 2;                            /*!< Level of verbosity desired [0-6]. Level 0 prints everything, 6 nothing. Level 2 or 3 is recommended for normal use */
        inline bool   timestamp     = false;                          /*!< Whether to put a timestamp on console outputs */
        inline size_t loglevel      = 2;                              /*!< Verbosity [0-6]. Level 0 prints everything, 6 nothing. Level 2 or 3 is recommended for normal use */
        inline size_t logh5pp       = 2;                              /*!< Verbosity of h5pp library [0-6] Level 2 or 3 is recommended for normal use */
}

    /*! \namespace settings::strategy Settings affecting the convergence rate of the xDMRG algorithm */
    namespace strategy {
        inline bool          bfgs_fix_rnorm_w_eigs       = true;                                   /*!< Use the eigenvalue solver for (H-E/L)² when BFGS returns with bad gradient */
        inline OptEigs       prefer_eigs_over_bfgs       = OptEigs::WHEN_SATURATED;                /*!< Prefer using the eigenvalue solver for (H-E/L)² over BFGS. Choose {WHEN_SATURATED, ALWAYS} */
        inline bool          expand_envs_when_stuck      = true;                                   /*!< Use environment expansion when stuck in local minima. alpha == lowest_variance */
        inline size_t        project_on_saturation       = 10;                                     /*!< Project to target axis/parity sector every nth iteration when saturated. (0 = turn off) */
        inline size_t        project_on_every_iter       = 5;                                      /*!< Project to target axis/parity sector at the end of every nth iteration. This implies doing it when stuck also. */
        inline bool          project_on_bond_update      = true;                                   /*!< Project to target axis/parity sector before the bond dimension limit is increased (only works if bond_increase_when == true). */
        inline bool          project_initial_state       = false;                                  /*!< Project to target axis/parity sector when initializing a state. */
        inline bool          project_final_state         = false;                                  /*!< Project to target axis/parity sector before writing down the final state */
        inline bool          randomize_on_bond_update    = true;                                   /*!< Randomize MPS by flipping random spins when growing the bond dimension */
        inline bool          randomize_early             = true;                                   /*!< Randomize MPS by flipping random spins before fully converging the first attempt (because the first attempt is biased) */
        inline bool          use_eigenspinors            = false;                                  /*!< Use random pauli-matrix eigenvectors when initializing each mps site along x,y or z  */
        inline size_t        max_resets                  = 1;                                      /*!< Maximum number of resets to product state due to saturation. One must be allowed for initialization */
        inline size_t        max_stuck_iters             = 5;                                      /*!< If stuck for this many iterations -> stop. */
        inline size_t        max_saturation_iters        = 5;                                      /*!< If either variance or entanglement saturated this long -> algorithm saturated = true */
        inline size_t        min_saturation_iters        = 1;                                      /*!< Saturated at least this many iterations before stopping */
        inline size_t        min_converged_iters         = 2;                                      /*!< Converged at least this many iterations before success */
        inline double        max_env_expansion_alpha     = 1e-4;                                   /*!< Maximum value of alpha used in environment expansion */
        inline size_t        multisite_mps_site_def      = 2;                                      /*!< Default number of sites in a multisite mps. More than ~8 is very expensive */
        inline size_t        multisite_mps_site_max      = 4;                                      /*!< Maximum number of sites in a multisite mps (used when stuck). More than ~8 is very expensive */
        inline MultisiteMove multisite_mps_move          = MultisiteMove::ONE;                     /*!< How many sites to move after a multi-site dmrg step, choose between {ONE, MID, MAX} */
        inline MultisiteWhen multisite_mps_when          = MultisiteWhen::OFF;                     /*!< When to increase the number of sites in a DMRG step {OFF, STUCK, SATURATED, ALWAYS} */
        inline std::string   target_sector               = "none";                                 /*!< Find an eigenstate in this parity sector. Choose between [random,randomAxis, none, x,+x,-x, y, +y,-y, z,+z,-z]  */
        inline std::string   initial_sector              = "random";                               /*!< Initialize state in this spin pattern/parity sector. Choose between [random, x,+x,-x, y, +y,-y, z,+z,-z]  */
        inline StateInitType initial_type                = StateInitType::REAL;                    /*!< Initial state can be REAL/CPLX */
        inline StateInit     initial_state               = StateInit::RANDOM_ENTANGLED_STATE;      /*!< Initial configuration for the spin chain (only for finite systems)  */
        inline StateInit     secondary_states            = StateInit::RANDOMIZE_PREVIOUS_STATE;    /*!< Spin configuration for subsequent states (only for finite systems)  */
        inline double        fes_rate                    = 2;                                      /*!< If |fes_rate| > 0, runs a finite entanglement scaling (fes) analysis with this step size in bond dimension, after finishing the main algorithm */
        inline UpdateWhen    bond_increase_when          = UpdateWhen::NEVER;                      /*!< If and when to increase the bond dimension limit {NEVER, TRUNCATED, STUCK, SATURATED, ITERATION}. */
        inline double        bond_increase_rate          = 8;                                      /*!< Bond dimension growth rate. Factor if 1<x<=2, constant shift if x > 2, otherwise invalid. */
        inline UpdateWhen    trnc_decrease_when          = UpdateWhen::NEVER;                      /*!< If and when to decrease SVD truncation error limit {NEVER, TRUNCATED, STUCK, SATURATED, ITERATION} */
        inline double        trnc_decrease_rate          = 1e-2;                                   /*!< Decrease SVD truncation error limit by this factor. Valid if 0 < x < 1 */
}


    /*! \namespace settings::precision Settings for the convergence threshold and precision of MPS, SVD and eigensolvers */
    namespace precision {
        inline size_t   eigs_max_iter                   = 1000  ;                  /*!< Maximum number of iterations for eigenvalue solver. */
        inline double   eigs_tolerance                  = 1e-14 ;                  /*!< Precision tolerance for halting the eigenvalue solver. */
        inline size_t   eigs_default_ncv                = 32    ;                  /*!< Parameter controlling the krylov/column space of the Arnoldi eigenvalue solver */
        inline size_t   bfgs_max_iter                   = 1000  ;                  /*!< Maximum number of iterations for the L-BFGS solver. */
        inline double   svd_truncation_lim              = 5e-32 ;                  /*!< Truncation error limit, i.e. discard singular values while the truncation error is lower than this */
        inline double   svd_truncation_init             = 1e-4 ;                   /*!< If truncation error limit is updated (trnc_decrease_when != NEVER), start from this value */
        inline size_t   svd_switchsize_bdc              = 16    ;                  /*!< Linear size of a matrix, below which SVD will use slower but more precise JacobiSVD instead of BDC (default is 16 , good could be ~64) */
        inline double   max_grad_tolerance              = 1e-8  ;                  /*!< Keep running an opimization step (BFGS/Arnoldi/GD+k) until max(∇log10(Var H)) < max_grad_tolerance */
        inline bool     use_compressed_mpo_squared_all  = false ;                  /*!< Use SVD to compress the bond dimensions of all H² mpos at the end of an iteration */
        inline bool     use_compressed_mpo_squared_otf  = true  ;                  /*!< Use SVD to compress the bond dimensions of the multisite H² mpo on-the-fly, just before an optimization step  */
        inline bool     use_mpo_energy_shift            = true  ;                  /*!< Whether to subtract E/L from ALL mpos to avoid catastrophic cancellation when computing the variance */
        inline bool     use_projection_on_mpo_squared   = true  ;                  /*!< Lift degeneracy by redefining H² --> (H² + Q(σ)) where Q(σ) = 0.5(1 - prod(σ)) = P(-σ) is the (flipped sign) projection operator */
        inline double   variance_convergence_threshold  = 1e-12 ;                  /*!< Desired precision on total energy variance. The MPS state is considered good enough when its energy variance reaches below this value */
        inline double   variance_saturation_sensitivity = 1e-2  ;                  /*!< Energy variance saturates when it stops changing. This sets the sensitivity to change. Good values are 1e-1 to 1e-4   */
        inline double   entropy_saturation_sensitivity  = 1e-6  ;                  /*!< Entanglement entropy saturates when it stops changing. This sets the sensitivity to change. Good values are 1e-3 to 1e-8   */
        inline double   target_subspace_error           = 1e-10 ;                  /*!< The target subspace error 1-Σ|<ϕ_i|ψ>|². Eigenvectors are found until reaching this value. Measures whether the incomplete basis of eigenstates spans the current state. */
        inline size_t   max_subspace_size               = 256   ;                  /*!< Maximum number of candidate eigenstates to keep for a subspace optimization step */
        inline long     max_size_full_diag              = 2048  ;                  /*!< Maximum problem size allowed for full diagonalization of the local (effective) hamiltonian matrix. */
        inline long     max_size_part_diag              = 4096  ;                  /*!< Maximum problem size allowed for partial diagonalization of the local (effective) hamiltonian matrix. */
        inline long     max_size_multisite              = 131072;                  /*!< Maximum problem size for multisite dmrg. If the linear size is larger than this, the algorithm prefers 1-site dmrg. */
        inline double   max_norm_error                  = 1e-10 ;                  /*!< Maximum norm deviation from unity during integrity checks */
    }


    /*! \namespace settings::model Settings for the Hamiltonian spin-model */
    namespace model {
        inline ModelType    model_type = ModelType::ising_tf_rf;   /*!< Choice of model: {ising_tf_rf,ising_sdual, ising_majorana, lbit} */
        inline size_t       model_size = 16;                       /*!< Number of sites on the chain. Only relevant for finite algorithms: fDMRG and xDMRG */

        /*! \namespace settings::model::ising_tf_rf Settings for the Transverse-field Ising model with a random on-site field */
        namespace ising_tf_rf {
            inline double       J1         = 1;                 /*!< Ferromagnetic coupling for nearest neighbors.*/
            inline double       J2         = 1;                 /*!< Ferromagnetic coupling for next-nearest neighbors.*/
            inline double       h_tran     = 1;                 /*!< Transverse field strength */
            inline double       h_mean     = 0;                 /*!< Random field mean of distribution */
            inline double       h_wdth     = 0;                 /*!< Random field width of distribution */
            inline long         spin_dim   = 2;                 /*!< Spin dimension */
            inline std::string  distribution  = "uniform";      /*!< Random distribution for couplings and fields */
        }

        /*! \namespace settings::model::ising_sdual Settings for the Self-dual Ising model */
        namespace ising_sdual {
            inline double       lambda        = 0;              /*!< Lambda parameter related to next nearest neighbor coupling */
            inline double       delta         = 0;              /*!< Delta defined as log(J_mean) - log(h_mean). We get J_mean and h_mean by fixing max(J_mean,h_mean) = 1 */
            inline bool         parity_sep    = false;          /*!< Separation of +-Z parity sectors */
            inline long         spin_dim      = 2;              /*!< Spin dimension */
            inline std::string  distribution  = "uniform";      /*!< Random distribution for couplings and fields */
        }

        /*! \namespace settings::model::ising_majorana Settings for the Ising-Majorana model */
        namespace ising_majorana {
            inline double       g             = 0;              /*!< Interaction parameter for nearest ZZ and next-nearest XX neighbor coupling */
            inline double       delta         = 0;              /*!< Delta defined as log(J_mean) - log(h_mean). We get J_mean and h_mean by fixing delta = 2lnW, W = J_wdth = 1/h_wdth */
            inline bool         parity_sep    = false;          /*!< Separation of +-Z parity sectors */
            inline long         spin_dim      = 2;              /*!< Spin dimension */
            inline std::string  distribution  = "uniform";      /*!< Random distribution for couplings and fields */
        }

        /*! \namespace settings::model::lbit Settings for the l-bit Hamiltonian */
        namespace lbit {
            inline double       J1_mean       = 0;              /*!< Constant offset for on-site */
            inline double       J2_mean       = 0;              /*!< Constant offset for two-body interaction */
            inline double       J3_mean       = 0;              /*!< Constant offset for three-body interaction */
            inline double       J1_wdth       = 0.5;            /*!< Width of the distribution for on-site interactions */
            inline double       J2_wdth       = 0.5;            /*!< Width of the distribution for two-body interaction (st.dev. for normal distribution). */
            inline double       J3_wdth       = 0.5;            /*!< Width of the distribution for three-body interaction */
            inline double       J2_xcls       = 1;              /*!< Exp. decay rate of two-body interactions: exp(-|i-j|/J2_xcls) * J2_rand */
            inline size_t       J2_span       = -1ul;           /*!< Maximum allowed range for pairwise interactions, |i-j| <= J2_span. Use -1 for infinite. Note that J2_span + 1 MPOs are used */
            inline double       f_mixer       = 0.1;            /*!< Mixing factor for unitary transformation to real-space */
            inline size_t       u_layer       = 6;              /*!< Number of unitary 2-site layers which transform lbit <-> real spaces */
            inline long         spin_dim      = 2;              /*!< Spin dimension */
            inline std::string  distribution  = "uniform";      /*!< Random distribution for interaction strengths */
        }
    }

    /*! \namespace settings::idmrg Settings for the infinite DMRG algorithm */
    namespace idmrg {
        inline bool     on                  = false;                               /*!< Turns iDMRG simulation on/off. */
        inline size_t   max_iters           = 5000;                                /*!< Maximum number of iDMRG iterations before forced termination */
        inline long     bond_max            = 32;                                  /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline long     bond_init           = 16;                                  /*!< Initial bond dimension limit. Only used when bond_increase_when == true. */
        inline size_t   print_freq          = 1000;                                /*!< Print frequency for console output. In units of iterations.  (0 = off). */
    }


    /*! \namespace settings::itebd Settings for the imaginary-time infinite TEBD algorithm  */
    namespace itebd {
        inline bool      on                    = false;                            /*!< Turns iTEBD simulation on/off. */
        inline size_t    max_iters             = 100000;                           /*!< Maximum number of iTEBD iterations before forced termination */
        inline double    time_step_init_real   = 0.0;                              /*!< Real part of initial time step delta_t */
        inline double    time_step_init_imag   = 0.1;                              /*!< Imag part of initial time step delta_t */
        inline double    time_step_min         = 0.00001;                          /*!< (Absolute value) Minimum and final time step for iTEBD time evolution. */
        inline size_t    suzuki_order          = 1;                                /*!< Order of the suzuki trotter decomposition (1,2 or 4) */
        inline long      bond_max              = 8;                                /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline long      bond_init             = 4;                                /*!< Initial bond dimension limit. Only used when bond_increase_when == true. */
        inline size_t    print_freq            = 5000;                             /*!< Print frequency for console output. In units of iterations. (0 = off).*/
    }

    /*! \namespace settings::fdmrg Settings for the finite DMRG algorithm */
    namespace fdmrg {
        inline bool      on                  = false;                              /*!< Turns fdmrg simulation on/off. */
        inline size_t    max_iters           = 10;                                 /*!< Max number of iterations. One iterations moves L steps. */
        inline size_t    min_iters           = 4;                                  /*!< Min number of iterations. One iterations moves L steps. */
        inline long      bond_max            = 128;                                /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline long      bond_init           = 8;                                  /*!< Initial bond dimension limit. Only used when bond_increase_when == true. */
        inline size_t    print_freq          = 100;                                /*!< Print frequency for console output. In units of iterations. (0 = off). */
        inline bool      store_wavefn        = false;                              /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */
    }


    /*! \namespace settings::flbit Settings for the finite l-bit algorithm */
    namespace flbit {
        inline bool     on                      = false;                           /*!< Turns flbit simulation on/off. */
        inline size_t   max_iters               = 10000;                           /*!< Max number of iterations. One iterations moves L steps. */
        inline size_t   min_iters               = 4;                               /*!< Min number of iterations. One iterations moves L steps. */
        inline bool     use_swap_gates          = true;                            /*!< Use gate swapping for pairwise long-range interactions rather then building a large multisite operator */
        inline long     bond_max                = 1024;                            /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline long     bond_init               = 8;                               /*!< Initial bond dimension limit. Used during iter <= 1 or when bond_increase_when == true, or starting from an entangled state */
        inline double   time_start_real         = 1e-1;                            /*!< Starting time point (real) */
        inline double   time_start_imag         = 0;                               /*!< Starting time point (imag) */
        inline double   time_final_real         = 1e6;                             /*!< Finishing time point (real) */
        inline double   time_final_imag         = 0;                               /*!< Finishing time point (imag) */
        inline size_t   time_num_steps          = 500;                             /*!< Number of steps from start to finish. Start and final times are included */
        inline double   time_gate_id_threshold  = 1e-8;                            /*!< Skip time evo. gates if exp(-iHt) is ~ 1 within this threshold */
        inline size_t   print_freq              = 1;                               /*!< Print frequency for console output. In units of iterations. (0 = off). */
        inline bool     compute_lbit_length     = false;                           /*!< Calculate the characteristic length-scale of lbits */
        inline bool     compute_lbit_stats      = false;                           /*!< Calculate the statistics of characteristic length-scale for various u and f parameters */
        inline bool     store_wavefn            = false;                           /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */
        inline bool     save_swap_gates         = false;                           /*!< Saves the state and swap gates used later for benchmarks */
}

    /*! \namespace settings::xdmrg Settings for the finite excited-state DMRG algorithm */
    namespace xdmrg {
        inline bool       on                            = false;                   /*!< Turns xDMRG simulation on/off. */
        inline size_t     max_iters                     = 10;                      /*!< Max number of iterations. One iterations moves L steps. */
        inline size_t     min_iters                     = 4;                       /*!< Min number of iterations. One iterations moves L steps. */
        inline long       bond_max                      = 768;                     /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline long       bond_init                     = 8;                       /*!< Initial bond dimension limit. Used during iter <= 1 or when bond_increase_when == true, or starting from an entangled state */
        inline size_t     opt_overlap_iters             = 2;                       /*!< Number of initial iterations selecting the candidate state with best overlap to the current state */
        inline long       opt_overlap_bond_lim          = 16;                      /*!< Bond limit during initial OVERLAP optimization. set to <= 0 for unlimited */
        inline size_t     opt_subspace_iters            = 2;                       /*!< Number of iterations using the subspace optimization of variance, after the overlap iterations */
        inline long       opt_subspace_bond_lim         = 32;                      /*!< Bond limit during initial SUBSPACE optimization. set to <= 0 for unlimited */
        inline size_t     print_freq                    = 1;                       /*!< Print frequency for console output. In units of iterations. (0 = off). */
        inline double     energy_density_target         = 0.5;                     /*!< Target energy in [0-1], where 0.5 means middle of spectrum. */
        inline double     energy_density_window         = 0.05;                    /*!< Accept states inside of energy_tgt +- energy_dens_window. */
        inline size_t     max_states                    = 1;                       /*!< Max number of random states to find using xDMRG on a single disorder realization */
        inline bool       store_wavefn                  = false;                   /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */
        inline bool       finish_if_entanglm_saturated  = true;                    /*!< Finish early as soon as entanglement has saturated */
        inline bool       finish_if_variance_saturated  = false;                   /*!< Finish early as soon as energy variance has saturated */
    }
}
/* clang-format on */
