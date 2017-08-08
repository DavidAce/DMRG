//
// Created by david on 8/7/17.
//

#ifndef FINITE_DMRG_EIGEN_N_SETTINGS_H
#define FINITE_DMRG_EIGEN_N_SETTINGS_H
namespace settings{
    //Parmaters that control eigensolver and SVD precision
    extern int             chi          ;       /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
    extern int             eigSteps     ;       /*!< Maximum number of steps for eigenvalue solver. */
    extern double          eigThreshold ;       /*!< Minimum threshold for halting eigenvalue solver. */
    extern double          SVDThreshold ;       /*!< Minimum threshold value for keeping singular values. */

    //Parameters controlling iDMRG
    extern int             max_idmrg_length  ;        /*!< Final length of 1D quantum chain. */

    //Parameters controlling fDMRG
    extern int             max_fdmrg_sweeps  ;        /*!< Number sweeps along the 1D quantum chain. */
    extern int             max_fdmrg_length  ;         /*!< Number sweeps along the 1D quantum chain. */

    //Parameters controlling iTEBD
    extern int             max_itebd_steps   ;        /*!< Number of iTEBD iterations. */
    extern double          delta_t_itebd     ;        /*!< Time step for iTEBD time evolution.*/

    //Parameters controlling Finite-entanglement scaling (FES)
    extern int             max_fes_steps      ;       /*!< Number of FES iterations per chi-value. */
    extern int             min_fes_chi        ;
    extern int             max_fes_chi        ;
    extern int             num_fes_chi        ;       /*!< Number of chi values for FES. */
//  extern   bool            increasing_chi   ;        /*!< \todo remove */
//  extern   int             chi_max          ;         /*!< \todo remove */

    //Profiling
    extern bool            profiling          ;       /*!< Turns profiling options on/off. */
    extern int             time_prec          ;       /*!< Sets precision of time output. */

    //Verbosity level
    extern int             verbosity          ;       /*!< Three-level verbosity. */
    extern bool            graphics           ;       /*!< Whether to print the chain graphically at each iteration.*/

};
#endif //FINITE_DMRG_EIGEN_N_SETTINGS_H
