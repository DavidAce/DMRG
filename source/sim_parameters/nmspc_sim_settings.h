//
// Created by david on 8/7/17.
//

#ifndef DMRG_N_SETTINGS_H
#define DMRG_N_SETTINGS_H
#include <string>
#include <unordered_set>
#include <vector>
/*! \brief General settings like max iterations, time-step, precision, etc.*/

/*!
 *  \namespace settings
 *  This namespace contains settings such as time-step length, number of iterations and precision parameters for
 *  the different algorithms.
 */

class class_settings_reader;
namespace h5pp{
    class File;
}

enum class SimulationType{iDMRG,fDMRG, xDMRG, iTEBD};

namespace settings {
    extern void load_from_file(class_settings_reader &indata);
    extern void load_from_hdf5(h5pp::File &h5ppFile);

    namespace input{
        extern std::string input_filename;
        extern std::string input_file;
    }
    //Parameters for the model Hamiltonian
    namespace model {
        extern std::string  initial_state ;                   /*!< Choose initial state of the MPS: {upup, updown, GHZ(upup+downdown), W(updown+downup), rps (random product state), random_chi (random state with bond dimension chi, only for iDMRG!)} "cat" or "random". Default "rps". */
        extern std::string  model_type    ;                   /*!< Choice of model type: {tf_ising, tf_nn_ising, selfdual_tf_rf_ising} above*/
        extern int          seed          ;                   /*!< Seed for the random number generator if you use random fields in the Hamiltonian. */

        //Parameters for the transverse-field Ising model
        namespace tf_ising {
            extern double       J             ;                 /*!< Ferromagnetic coupling. J < 0  Gives a ferromagnet. J > 0 an antiferromagnet. */
            extern double       g             ;                 /*!< Transverse field strength */
            extern double       w             ;                 /*!< Randomness strength for the random field */
            extern int          d             ;                 /*!< Local dimension */
        }

        //Parameters for the transvese-field next-nearest neighbor Ising model
        namespace tf_nn_ising {
            extern double       J1            ;                 /*!< Ferromagnetic coupling for nearest neighbors.*/
            extern double       J2            ;                 /*!< Ferromagnetic coupling for next-nearest neighbors.*/
            extern double       g             ;                 /*!< Transverse field strength */
            extern double       w             ;                 /*!< Randomness strength for the random field */
            extern int          d             ;                 /*!< Local dimension */
        }

        //Parameters for the selfdual transvese-field random-field next-neighbor Ising model
        namespace selfdual_tf_rf_ising {
            extern double       J_log_mean    ;                 /*!< Average ferromagnetic coupling strength.*/
            extern double       h_log_mean    ;                 /*!< Average transverse magnetic field strength */
            extern double       J_sigma       ;                 /*!< Standard deviation for the lognormal distribution, i.e. = std(log(J)) , for the ferromagnetic coupling */
            extern double       h_sigma       ;                 /*!< Standard deviation for the lognormal distribution, i.e. = std(log(h))   for the transverse magnetic field */
            extern double       lambda        ;                 /*!< Lambda parameter */
            extern int          d             ;                 /*!< Local dimension */
        }
    }

    //Parmaters that control MPS, eigensolver and SVD precision
    namespace precision {
        extern int      eigMaxIter   ;                      /*!< Maximum number of steps for eigenvalue solver. */
        extern double   eigThreshold ;                      /*!< Minimum threshold for halting eigenvalue solver. */
        extern int      eigMaxNcv  ;                        /*!< Parameter controlling the column space? of the Lanczos solver. */
        extern double   SVDThreshold ;                      /*!< Minimum threshold value for keeping singular values. */
        extern double   VarConvergenceThreshold ;           /*!< Variance convergence threshold. The MPS state is considered good enough when its variance reaches below this value */
        extern double   VarSaturationThreshold ;            /*!< Variance saturation  threshold. The variance has saturated when its (absolute) slope reaches below this value */
        extern double   EntEntrSaturationThreshold;         /*!< Entanglement Entropy saturation threshold. The entanglement entropy has saturated when its (absolute) slope reaches below this value*/
        extern int      MaxSizeFullDiag;                    /*!< Maximum linear size allowed for full diagonalization of the local hamiltonian matrix. Use 0 to allow any value */
    }

    //Parameters controlling iDMRG
    namespace idmrg {
        extern bool on         ;                            /*!< Turns iDMRG simulation on/off. */
        extern int  max_steps  ;                             /*!< Final length of 1D quantum chain. */
        extern long chi_max    ;                            /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        extern bool chi_grow   ;                            /*!< Whether to increase chi slowly up to chi_max or go up to chi_max directly. */
        extern int  print_freq ;                            /*!< Print frequency for console output. (0 = off). */
        extern int  store_freq ;                            /*!< Store frequency,for output file buffer. (0 = off). */

    }
    //Parameters controlling fDMRG
    namespace fdmrg {
        extern bool on           ;                          /*!< Turns fDMRG simulation on/off. */
        extern int  num_sites   ;                          /*!< Number sweeps along the 1D quantum chain. */
        extern int  max_sweeps   ;                          /*!< Max number sweeps along the 1D quantum chain. */
        extern int  min_sweeps   ;                          /*!< Min number sweeps along the 1D quantum chain. */
        extern long chi_max      ;                          /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        extern bool chi_grow   ;                            /*!< Whether to increase chi slowly up to chi_max or go up to chi_max directly. */
        extern int  print_freq   ;                          /*!< Print frequency for console output. In units of sweeps. (0 = off). */
        extern int  store_freq ;                            /*!< Store frequency,for output file buffer. In units of sweeps. (0 = off). */
        extern bool store_wavefn ;                         /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */

    }

    //Parameters controlling xDMRG
    namespace xdmrg {
        extern bool    on           ;                       /*!< Turns xDMRG simulation on/off. */
        extern int     num_sites    ;                       /*!< Number sweeps along the 1D quantum chain. */
        extern int     max_sweeps   ;                       /*!< Number sweeps along the 1D quantum chain. */
        extern int     min_sweeps   ;                          /*!< Min number sweeps along the 1D quantum chain. */
        extern long    chi_max      ;                       /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        extern bool    chi_grow     ;                       /*!< Whether to increase chi slowly up to chi_max or go up to chi_max directly. */
        extern int     print_freq   ;                       /*!< Print frequency for console output. In units of sweeps. (0 = off). */
        extern int     store_freq   ;                       /*!< Store frequency,for output file buffer. In units of sweeps. (0 = off). */
        extern bool    store_wavefn ;                       /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */
        extern double  energy_density;                      /*!< Target energy in [0-1], where 0.5 means middle of spectrum. */
        extern double  energy_window;                       /*!< Accept states inside of energy_target +- energy_window. */

    }

    //Parameters controlling iTEBD
    namespace itebd {
        extern bool on              ;                       /*!< Turns iTEBD simulation on/off. */
        extern int      max_steps   ;                       /*!< Number of iTEBD iterations, after which the simulation terminates regardless of convergence. Set high.*/
        extern double   delta_t0    ;                       /*!< Initial time step for iTEBD time evolution.*/
        extern double   delta_tmin  ;                       /*!< Final time step for iTEBD time evolution.*/
        extern int      suzuki_order;                       /*!< Order of the suzuki trotter decomposition (1,2 or 4) */
        extern long     chi_max     ;                       /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        extern bool     chi_grow    ;                       /*!< Whether to increase chi slowly up to chi_max or go up to chi_max directly. */
        extern int      print_freq  ;                       /*!< Print frequency for console output. (0 = off).*/
        extern int      store_freq  ;                       /*!< Store frequency,for output file buffer. (0 = off). */

    }

    namespace hdf5 {
        extern bool         save_to_file            ;        /*!< If true, saves the simulation data to an HDF5 file instead of just outputting to console */
        extern bool         save_progress           ;        /*!< If true, saves the simulation data periodically */
        extern std::string  access_mode             ;        /*!< Choose access mode to the file. Choose between READWRITE, READONLY */
        extern std::string  create_mode             ;        /*!< Choose access mode to the file. Choose between TRUNCATE, OPEN, RENAME */
        extern std::string  output_filename         ;        /*!< Name of the output HDF5 file relative to the execution point  */
        extern bool         full_storage            ;        /*!< If true, saves the full MPS to file. Set to false to reduce output file size. */
        extern bool         store_profiling         ;        /*!< Whether to store profiling information to file. */
    }
    //Profiling
    namespace profiling {
        extern bool     on        ;                  /*!< If true, turns on profiling and timings will be shown on console. */
        extern int      precision ;                  /*!< Sets precision (number of decimals) of time output. */
    }
    //Console settings
    namespace console {
        extern int  verbosity ;                      /*!< Level of verbosity desired [0-2]. Level 0 prints almost nothing, level 2 prints everything */
        extern bool timestamp ;                      /*!< Whether to put a timestamp on console outputs */
    }
}
#endif //DMRG_N_SETTINGS_H
