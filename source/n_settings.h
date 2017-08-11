//
// Created by david on 8/7/17.
//

#ifndef FINITE_DMRG_EIGEN_N_SETTINGS_H
#define FINITE_DMRG_EIGEN_N_SETTINGS_H

namespace settings{
    //Parmaters that control eigensolver and SVD precision
    inline int             chi                 = 10;           /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
    inline int             eigSteps            = 1000;         /*!< Maximum number of steps for eigenvalue solver. */
    inline double          eigThreshold        = 1e-12;        /*!< Minimum threshold for halting eigenvalue solver. */
    inline double          SVDThreshold        = 1e-12;        /*!< Minimum threshold value for keeping singular values. */

    //Parameters controlling iDMRG
    inline int             idmrg_max_length    = 50;           /*!< Final length of 1D quantum chain. */

    //Parameters controlling fDMRG
    inline int             fdmrg_max_sweeps    = 1;             /*!< Number sweeps along the 1D quantum chain. */
    inline int             fdmrg_max_length    = 50;            /*!< Number sweeps along the 1D quantum chain. */


    //Parameters controlling iTEBD
    inline int             itebd_max_steps     = 5000;            /*!< Number of iTEBD iterations. */
    inline double          itebd_delta_t       = 0.01;            /*!< Time step for iTEBD time evolution.*/

    //Parameters controlling Finite-entanglement scaling (FES)
    inline int             fes_max_steps       = 1000;         /*!< Number of FES iterations per chi-value. */
    inline int             fes_min_chi         = 10;
    inline int             fes_max_chi         = 20;
    inline int             fes_num_chi         = 3;           /*!< Number of chi values for FES. */

    //Save data to hdf5
    inline bool            hdf5_save_to_file               = true;
    inline bool            hdf5_create_dir_if_not_found    = true;
    inline std::string     hdf5_filename                   = "data.h5";
    inline std::string     hdf5_path                       = "../output";

    //Profiling
    inline bool            profiling_on                    = false;        /*!< Turns profiling options on/off. */
    inline int             profiling_precision             = 5;            /*!< Sets precision of time output. */

    //Verbosity level
    inline int             console_verbosity           = 0;            /*!< Three-level verbosity. */
    inline bool            console_graphics            = true;         /*!< Whether to print the chain graphically at each iteration.*/

};
#endif //FINITE_DMRG_EIGEN_N_SETTINGS_H
