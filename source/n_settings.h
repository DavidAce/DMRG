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
namespace settings {

    //Parmaters that control eigensolver and SVD precision
    namespace precision {
        inline int      eigSteps        = 5000;         /*!< Maximum number of steps for eigenvalue solver. */
        inline double   eigThreshold    = 1e-12;        /*!< Minimum threshold for halting eigenvalue solver. */
        inline int      eig_max_ncv     = 20;
        inline double   SVDThreshold    = 1e-12;        /*!< Minimum threshold value for keeping singular values. */
    }

    //Parameters controlling iDMRG
    namespace idmrg {
        inline int  max_length   = 200;        /*!< Final length of 1D quantum chain. */
        inline long chi_max      = 25;           /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */

    }
    //Parameters controlling fDMRG
    namespace fdmrg {
        inline int  max_length   = 200;            /*!< Number sweeps along the 1D quantum chain. */
        inline int  max_sweeps   = 2;             /*!< Number sweeps along the 1D quantum chain. */
        inline long chi_max      = 8;
    }

    //Parameters controlling iTEBD
    namespace itebd {
        inline int      max_steps   = 10000;            /*!< Number of iTEBD iterations. */
        inline double   delta_t     = 0.01;             /*!< Time step for iTEBD time evolution.*/
        inline long     chi_max     = 25;
    }
    //Parameters controlling Finite-entanglement scaling (FES) in iTEBD-mode.
    namespace fes_itebd {
        inline int       max_steps    = 300000;         /*!< Number of FES iterations per chi-value. */
        inline double    delta_t      = 0.01;          /*!< Time step for iTEBD time evolution.*/
        inline long      chi_min      = 4;
        inline long      chi_max      = 8;
        inline long      chi_num      = 3;           /*!< Number of chi values for FES. */
    }
    //Parameters controlling Finite-entanglement scaling (FES) in iDMRG-mode.
    namespace fes_idmrg {
        inline int       max_steps    = 5000;          /*!< Number of FES iterations per chi-value. */
        inline long      chi_min      = 4;
        inline long      chi_max      = 8;
        inline long      chi_num      = 3;           /*!< Number of chi values for FES. */
    }


    //Save data to hdf5
    namespace hdf5 {
        inline bool         save_to_file               = true;
        inline bool         create_dir_if_not_found    = true;
        inline std::string  filename                   = "data.h5";
        inline std::string  path                       = "../output";
        inline bool         full_storage               = true;           /*!< If true, saves more simulation data to file (such as explicit form of MPS). Set to false to reduce output file size. */
    }
    //Profiling
    namespace profiling {
        inline bool     on        = false;        /*!< Turns profiling options on/off. */
        inline int      precision = 5;            /*!< Sets precision of time output. */
    }
    //Console settings
    namespace console {
        inline int  verbosity = 0;            /*!< Three-level verbosity. */
        inline bool graphics = true;          /*!< Whether to print the chain graphically at each iteration.*/
    }
};
#endif //DMRG_N_SETTINGS_H
