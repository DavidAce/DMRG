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
    extern int             idmrg_max_length  ;        /*!< Final length of 1D quantum chain. */

    //Parameters controlling fDMRG
    extern int             fdmrg_max_sweeps  ;        /*!< Number sweeps along the 1D quantum chain. */
    extern int             fdmrg_max_length  ;         /*!< Number sweeps along the 1D quantum chain. */

    //Parameters controlling iTEBD
    extern int             itebd_max_steps   ;        /*!< Number of iTEBD iterations. */
    extern double          itebd_delta_t     ;        /*!< Time step for iTEBD time evolution.*/

    //Parameters controlling Finite-entanglement scaling (FES)
    extern int             fes_max_steps      ;       /*!< Number of FES iterations per chi-value. */
    extern int             fes_min_chi        ;
    extern int             fes_max_chi        ;
    extern int             fes_num_chi        ;       /*!< Number of chi values for FES. */

    //Save data to hdf5
    extern bool hdf5_save_to_file                    ;
    extern std::string hdf5_filename                 ;
    extern std::string hdf5_path                     ;

    //Profiling
    extern bool            profiling_on          ;       /*!< Turns profiling options on/off. */
    extern int             profiling_precision          ;       /*!< Sets precision of time output. */

    //Verbosity level
    extern int             console_verbosity          ;       /*!< Three-level verbosity. */
    extern bool            console_graphics           ;       /*!< Whether to print the chain graphically at each iteration.*/

};
#endif //FINITE_DMRG_EIGEN_N_SETTINGS_H
