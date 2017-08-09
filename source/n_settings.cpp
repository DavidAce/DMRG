//
// Created by david on 8/7/17.
//
#include <iostream>
namespace settings{
    //Parmaters that control eigensolver and SVD precision
    int             chi                 = 20;           /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
    int             eigSteps            = 1000;         /*!< Maximum number of steps for eigenvalue solver. */
    double          eigThreshold        = 1e-10;        /*!< Minimum threshold for halting eigenvalue solver. */
    double          SVDThreshold        = 1e-10;        /*!< Minimum threshold value for keeping singular values. */

    //Parameters controlling iDMRG
    int             idmrg_max_length    = 50;           /*!< Final length of 1D quantum chain. */

    //Parameters controlling fDMRG
    int             fdmrg_max_sweeps    = 4;            /*!< Number sweeps along the 1D quantum chain. */
    int             fdmrg_max_length    = 50;            /*!< Number sweeps along the 1D quantum chain. */


    //Parameters controlling iTEBD
    int             itebd_max_steps     = 5000;         /*!< Number of iTEBD iterations. */
    double          itebd_delta_t       = 0.005;        /*!< Time step for iTEBD time evolution.*/

    //Parameters controlling Finite-entanglement scaling (FES)
    int             fes_max_steps       = 2000;         /*!< Number of FES iterations per chi-value. */
    int             fes_min_chi         = 5;
    int             fes_max_chi         = 50;
    int             fes_num_chi         = 5;           /*!< Number of chi values for FES. */

    //Save data to hdf5
    bool            hdf5_save_to_file               = true;
    bool            hdf5_create_dir_if_not_found    = true;
    std::string     hdf5_filename                   = "data.h5";
    std::string     hdf5_path                       = "../output";

    //Profiling
    bool            profiling_on        = false;        /*!< Turns profiling options on/off. */
    int             profiling_precision = 5;            /*!< Sets precision of time output. */

    //Verbosity level
    int             console_verbosity           = 0;            /*!< Three-level verbosity. */
    bool            console_graphics            = true;         /*!< Whether to print the chain graphically at each iteration.*/



};