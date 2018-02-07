//
// Created by david on 8/7/17.
//

#ifndef DMRG_N_SETTINGS_H
#define DMRG_N_SETTINGS_H
#include <string>
/*! \brief General settings like max iterations, time-step, precision, etc.*/

/*!
 *  \namespace settings
 *  This namespace contains settings such as time-step length, number of iterations and precision parameters for
 *  the different algorithms.
 *
 */

class class_file_reader;

namespace settings {
    extern void load_from_file(class_file_reader &indata);


    //Parmaters that control eigensolver and SVD precision
    namespace precision {
        extern int      eigSteps     ;                      /*!< Maximum number of steps for eigenvalue solver. */
        extern double   eigThreshold ;                      /*!< Minimum threshold for halting eigenvalue solver. */
        extern int      eig_max_ncv  ;                      /*!< Parameter controlling the column space? of the Lanczos solver. */
        extern double   SVDThreshold ;                      /*!< Minimum threshold value for keeping singular values. */
    }

    //Parameters controlling iDMRG
    namespace idmrg {
        extern bool on         ;                            /*!< Turns iDMRG simulation on/off. */
        extern int  max_length ;                            /*!< Final length of 1D quantum chain. */
        extern long chi_max    ;                            /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        extern int  print_freq ;                            /*!< Print frequency for console output. (0 = off). */

    }
    //Parameters controlling fDMRG
    namespace fdmrg {
        extern bool on           ;                          /*!< Turns fDMRG simulation on/off. */
        extern int  max_length   ;                          /*!< Number sweeps along the 1D quantum chain. */
        extern int  max_sweeps   ;                          /*!< Number sweeps along the 1D quantum chain. */
        extern long chi_max      ;                          /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        extern int  print_freq   ;                          /*!< Print frequency for console output. (0 = off). */
    }

    //Parameters controlling iTEBD
    namespace itebd {
        extern bool on              ;                       /*!< Turns iTEBD simulation on/off. */
        extern int      max_steps   ;                       /*!< Number of iTEBD iterations, after which the simulation terminates regardless of convergence. Set high.*/
        extern double   delta_t0    ;                       /*!< Initial time step for iTEBD time evolution.*/
        extern double   delta_tmin  ;                       /*!< Final time step for iTEBD time evolution.*/
        extern int      suzuki_order;                       /*!< Order of the suzuki trotter decomposition (1,2 or 4) */
        extern long     chi_max     ;                       /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        extern int      print_freq  ;                       /*!< Print frequency for console output. (0 = off).*/
    }
    //Parameters controlling Finite-entanglement scaling (FES) in iTEBD-mode.
    namespace fes_itebd {
        extern bool      on          ;                      /*!< Turns FES-ITEBD simulation on/off. */
        extern int       max_steps   ;                      /*!< Number of iTEBD iterations, after which the simulation terminates regardless of convergence. Set high.*/
        extern double    delta_t0    ;                      /*!< Initial time step for iTEBD time evolution.*/
        extern double    delta_tmin  ;                      /*!< Final time step for iTEBD time evolution.*/
        extern int       suzuki_order;                      /*!< Order of the suzuki trotter decomposition (1,2 or 4) */
        extern long      chi_min     ;                      /*!< Minimum chi-value in range. */
        extern long      chi_max     ;                      /*!< Maximum chi-value in range. */
        extern long      chi_num     ;                      /*!< Number of chi values for in range. */
        extern int       print_freq  ;                      /*!< Print frequency for console output. (0 = off).*/

    }
    //Parameters controlling Finite-entanglement scaling (FES) in iDMRG-mode.
    namespace fes_idmrg {
        extern bool on               ;                      /*!< Turns FES-iDMRG simulation on/off. */
        extern int       max_length  ;                      /*!< Number of FES iterations per chi-value. */
        extern long      chi_min     ;                      /*!< Minimum chi-value in range. */
        extern long      chi_max     ;                      /*!< Maximum chi-value in range. */
        extern long      chi_num     ;                      /*!< Number of chi values for in range. */
        extern int       print_freq  ;                       /*!< Print frequency for console output. (0 = off).*/

    }


    //Save data to hdf5 (NOT FULLY IMPLEMENTED YET)
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
