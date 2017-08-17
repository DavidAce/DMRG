//
// Created by david on 8/7/17.
//

#ifndef FINITE_DMRG_EIGEN_N_SETTINGS_H
#define FINITE_DMRG_EIGEN_N_SETTINGS_H
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
        inline int      eigSteps        = 1000;                /*!< Maximum number of steps for eigenvalue solver. */
        inline double   eigThreshold    = 1e-12;        /*!< Minimum threshold for halting eigenvalue solver. */
        inline double   SVDThreshold    = 1e-8;        /*!< Minimum threshold value for keeping singular values. */
    }

    //Parameters controlling iDMRG
    namespace idmrg {
        inline int max_length   = 50;        /*!< Final length of 1D quantum chain. */
        inline int max_chi      = 20;           /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
    }
    //Parameters controlling fDMRG
    namespace fdmrg {
        inline int max_sweeps   = 1;             /*!< Number sweeps along the 1D quantum chain. */
        inline int max_length   = 50;            /*!< Number sweeps along the 1D quantum chain. */
        inline int max_chi      = 20;
    }

    //Parameters controlling iTEBD
    namespace itebd {
        inline int      max_steps   = 10000;            /*!< Number of iTEBD iterations. */
        inline double   delta_t     = 0.01;                /*!< Time step for iTEBD time evolution.*/
        inline int      max_chi     = 25;
    }
    //Parameters controlling Finite-entanglement scaling (FES)
    namespace fes {
        inline int      max_steps    = 10000;         /*!< Number of FES iterations per chi-value. */
        inline double   delta_t      = 0.1;                /*!< Time step for iTEBD time evolution.*/
        inline int      min_chi      = 5;
        inline int      max_chi      = 25;
        inline int      num_chi      = 5;           /*!< Number of chi values for FES. */
    }
    //Save data to hdf5
    namespace hdf5 {
        inline bool         save_to_file               = true;
        inline bool         create_dir_if_not_found    = true;
        inline std::string  filename                   = "data.h5";
        inline std::string  path                       = "../output";
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
#endif //FINITE_DMRG_EIGEN_N_SETTINGS_H
