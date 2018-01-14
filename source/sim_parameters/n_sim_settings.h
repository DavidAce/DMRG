//
// Created by david on 8/7/17.
//

#ifndef DMRG_N_SETTINGS_H
#define DMRG_N_SETTINGS_H
#include <string>
#include <IO/class_file_reader.h>
/*! \brief General settings like max iterations, time-step, precision, etc.*/

/*!
 *  \namespace settings
 *  This namespace contains settings such as time-step length, number of iterations and precision parameters for
 *  the different algorithms.
 *
 */



namespace settings {
    extern void initialize(class_file_reader &indata);

    //Parmaters that control eigensolver and SVD precision
    namespace precision {
        extern int      eigSteps     ;                /*!< Maximum number of steps for eigenvalue solver. */
        extern double   eigThreshold ;                /*!< Minimum threshold for halting eigenvalue solver. */
        extern int      eig_max_ncv  ;                /*!< Parameter controlling the column space? of the Lanczos solver. */
        extern double   SVDThreshold ;                /*!< Minimum threshold value for keeping singular values. */
    }

    //Parameters controlling iDMRG
    namespace idmrg {
        extern bool on         ;
        extern int  max_length ;                     /*!< Final length of 1D quantum chain. */
        extern long chi_max    ;                     /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */

    }
    //Parameters controlling fDMRG
    namespace fdmrg {
        extern bool on           ;
        extern int  max_length   ;                   /*!< Number sweeps along the 1D quantum chain. */
        extern int  max_sweeps   ;                   /*!< Number sweeps along the 1D quantum chain. */
        extern long chi_max      ;
    }

    //Parameters controlling iTEBD
    namespace itebd {
        extern bool on           ;
        extern int      max_steps;                   /*!< Number of iTEBD iterations. */
        extern double   delta_t  ;                   /*!< Time step for iTEBD time evolution.*/
        extern long     chi_max  ;
    }
    //Parameters controlling Finite-entanglement scaling (FES) in iTEBD-mode.
    namespace fes_itebd {
        extern bool on             ;
        extern int       max_steps ;                /*!< Number of FES iterations per chi-value. */
        extern double    delta_t   ;                /*!< Time step for iTEBD time evolution.*/
        extern long      chi_min   ;
        extern long      chi_max   ;
        extern long      chi_num   ;                /*!< Number of chi values for FES. */
    }
    //Parameters controlling Finite-entanglement scaling (FES) in iDMRG-mode.
    namespace fes_idmrg {
        extern bool on               ;
        extern int       max_steps   ;             /*!< Number of FES iterations per chi-value. */
        extern long      chi_min     ;
        extern long      chi_max     ;
        extern long      chi_num     ;             /*!< Number of chi values for FES. */
    }


    //Save data to hdf5 (NOT FULLY IMPLEMENTED YET)
    namespace hdf5 {
        extern bool         save_to_file            ;
        extern bool         create_dir_if_not_found ;
        extern std::string  filename                ;
        extern std::string  path                    ;
        extern bool         full_storage            ;    /*!< If true, saves more simulation data to file (such as explicit form of MPS). Set to false to reduce output file size. */
    }
    //Profiling
    namespace profiling {
        extern bool     on        ;             /*!< Turns profiling options on/off. */
        extern int      precision ;             /*!< Sets precision of time output. */
    }
    //Console settings
    namespace console {
        extern int  verbosity ;                  /*!< Three-level verbosity. */
    }
};
#endif //DMRG_N_SETTINGS_H
