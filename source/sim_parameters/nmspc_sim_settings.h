//
// Created by david on 8/7/17.
//

#ifndef DMRG_N_SETTINGS_H
#define DMRG_N_SETTINGS_H
#include <string>
#include <unordered_set>
/*! \brief General settings like max iterations, time-step, precision, etc.*/

/*!
 *  \namespace settings
 *  This namespace contains settings such as time-step length, number of iterations and precision parameters for
 *  the different algorithms.
 */

class class_file_reader;

enum class SimulationType{iDMRG,fDMRG, xDMRG, iTEBD};

namespace settings {
    extern void load_from_file(class_file_reader &indata);
    //Parameters for the model Hamiltonian
    namespace model {
        extern std::string  initial_state ;                   /*!< Choose initial state of the MPS: {upup, updown, GHZ(upup+downdown), W(updown+downup), rps (random product state), random_chi (random state with bond dimension chi, only for iDMRG!)} "cat" or "random". Default "rps". */
        extern std::string  model_type    ;                   /*!< Choice of model type: {tf_ising, tf_nn_ising, selfdual_tf_rf_ising} above*/

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
            extern double       J_mu          ;                 /*!< Average ferromagnetic coupling strength.*/
            extern double       h_mu          ;                 /*!< Average transverse magnetic field strength */
            extern double       J_sigma       ;                 /*!< Standard deviation for the lognormal distribution, i.e. = std(log(J)) , for the ferromagnetic coupling */
            extern double       h_sigma       ;                 /*!< Standard deviation for the lognormal distribution, i.e. = std(log(h))   for the transverse magnetic field */
            extern double       lambda        ;                 /*!< Lambda parameter */
            extern int          d             ;                 /*!< Local dimension */
        }
    }

    //Parmaters that control eigensolver and SVD precision
    namespace precision {
        extern int      eigMaxIter   ;                      /*!< Maximum number of steps for eigenvalue solver. */
        extern double   eigThreshold ;                      /*!< Minimum threshold for halting eigenvalue solver. */
        extern int      eigMaxNcv  ;                        /*!< Parameter controlling the column space? of the Lanczos solver. */
        extern double   SVDThreshold ;                      /*!< Minimum threshold value for keeping singular values. */
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
        extern int  max_length   ;                          /*!< Number sweeps along the 1D quantum chain. */
        extern int  max_sweeps   ;                          /*!< Number sweeps along the 1D quantum chain. */
        extern long chi_max      ;                          /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        extern bool chi_grow   ;                            /*!< Whether to increase chi slowly up to chi_max or go up to chi_max directly. */
        extern int  print_freq   ;                          /*!< Print frequency for console output. In units of sweeps. (0 = off). */
        extern int  store_freq ;                            /*!< Store frequency,for output file buffer. In units of sweeps. (0 = off). */

    }

    //Parameters controlling xDMRG
    namespace xdmrg {
        extern bool    on           ;                       /*!< Turns xDMRG simulation on/off. */
        extern int     max_length   ;                       /*!< Number sweeps along the 1D quantum chain. */
        extern int     max_sweeps   ;                       /*!< Number sweeps along the 1D quantum chain. */
        extern long    chi_max      ;                       /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        extern bool    chi_grow     ;                       /*!< Whether to increase chi slowly up to chi_max or go up to chi_max directly. */
        extern int     seed         ;                       /*!< Seed for the random number generator if you use random fields in the Hamiltonian. */
        extern int     print_freq   ;                       /*!< Print frequency for console output. In units of sweeps. (0 = off). */
        extern int     store_freq   ;                       /*!< Store frequency,for output file buffer. In units of sweeps. (0 = off). */

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

    //Save data_struct to hdf5 (NOT FULLY IMPLEMENTED YET)
    namespace hdf5 {
        extern bool         save_to_file            ;        /*!< If true, saves the simulation data to an HDF5 file instead of just outputting to console */
        extern bool         create_dir_if_not_found ;        /*!< If true, an output directory will be created in the project root folder if it isn't found */
        extern bool         overwrite_file_if_found ;        /*!< If true, an hdf5-file with the provided filename will be overwritten if found in output_folder */
        extern std::string  output_filename         ;        /*!< Name of the output HDF5 file */
        extern std::string  output_folder           ;        /*!< Path of the output HDF5 file */
        extern bool         full_storage            ;        /*!< If true, saves more simulation data to file (such as explicit form of MPS). Set to false to reduce output file size. */
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
};
#endif //DMRG_N_SETTINGS_H
