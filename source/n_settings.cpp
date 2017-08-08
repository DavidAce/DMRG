//
// Created by david on 8/7/17.
//

namespace settings{
    //Parmaters that control eigensolver and SVD precision
    int             chi                 = 20;           /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
    int             eigSteps            = 1000;         /*!< Maximum number of steps for eigenvalue solver. */
    double          eigThreshold        = 1e-10;        /*!< Minimum threshold for halting eigenvalue solver. */
    double          SVDThreshold        = 1e-12;        /*!< Minimum threshold value for keeping singular values. */

    //Parameters controlling iDMRG
    int             max_idmrg_length    = 50;           /*!< Final length of 1D quantum chain. */

    //Parameters controlling fDMRG
    int             max_fdmrg_sweeps    = 4;            /*!< Number sweeps along the 1D quantum chain. */
    int             max_fdmrg_length    = 50;            /*!< Number sweeps along the 1D quantum chain. */


    //Parameters controlling iTEBD
    int             max_itebd_steps     = 5000;         /*!< Number of iTEBD iterations. */
    double          delta_t_itebd       = 0.005;        /*!< Time step for iTEBD time evolution.*/

    //Parameters controlling Finite-entanglement scaling (FES)
    int             max_fes_steps       = 2000;         /*!< Number of FES iterations per chi-value. */
    int             min_fes_chi         = 5;
    int             max_fes_chi         = 50;
    int             num_fes_chi         = 5;           /*!< Number of chi values for FES. */
//    bool            increasing_chi      = false;        /*!< \todo remove */
//    int             chi_max             = 50;           /*!< \todo remove */


    //Profiling
    bool            profiling           = false;        /*!< Turns profiling options on/off. */
    int             time_prec           = 5;            /*!< Sets precision of time output. */

    //Verbosity level
    int             verbosity           = 0;            /*!< Three-level verbosity. */
    bool            graphics            = true;         /*!< Whether to print the chain graphically at each iteration.*/

};