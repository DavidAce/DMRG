//
// Created by david on 8/7/17.
//

#ifndef DMRG_N_SETTINGS_H
#define DMRG_N_SETTINGS_H
#include <string>
#include <unordered_set>
#include <vector>
/*! \brief General settings like max iterations, time-step, precision, etc.*/

/*!
 *  \namespace settings
 *  This namespace contains settings such as time-step length, number of iterations and precision parameters for
 *  the different algorithms.
 */

class class_settings_reader;
namespace h5pp{
    class File;
}

enum class SimulationType{iDMRG,fDMRG, xDMRG, iTEBD};
enum class StorageLevel:size_t {NONE,LIGHT,NORMAL,FULL};

namespace settings {
    extern void load_from_file(class_settings_reader &indata);
    extern void load_from_hdf5(h5pp::File &h5ppFile);

    namespace input{
        inline std::string input_file = "input/input.cfg";
        inline std::string input_filename = "input.cfg";
    }
    //Parameters for the model Hamiltonian
    namespace model {
        inline std::string  initial_state  = "tf_ising";        /*!< Choose initial state of the MPS: {upup, updown, GHZ(upup+downdown), W(updown+downup), rps (random product state), random_chi (random state with bond dimension chi, only for iDMRG!)} "cat" or "random". Default "rps". */
        inline std::string  model_type     = "rps";             /*!< Choice of model type: {tf_ising, tf_nn_ising, selfdual_tf_rf_ising} above*/
        inline int          seed_init_mpo  = 1;                 /*!< Seed for the random number generator if you use random fields in the Hamiltonian. */
        inline int          seed_init_mps  = 1;                 /*!< Seed for the random number generator when selecting the initial random product state. */
        inline std::string  symmetry       = "sx";              /*!< Initialize in parity symmetry sector: {sx,sy,sz,random,none} */


        //Parameters for the transverse-field Ising model
        namespace tf_ising {
            inline double       J  = 1;                         /*!< Ferromagnetic coupling. J < 0  Gives a ferromagnet. J > 0 an antiferromagnet. */
            inline double       g  = 1;                         /*!< Transverse field strength */
            inline double       w  = 0;                         /*!< Randomness strength for the random field */
            inline int          d  = 2;                         /*!< Local dimension */
        }

        //Parameters for the transvese-field next-nearest neighbor Ising model
        namespace tf_nn_ising {
            inline double       J1  = 1;                         /*!< Ferromagnetic coupling for nearest neighbors.*/
            inline double       J2  = 1;                         /*!< Ferromagnetic coupling for next-nearest neighbors.*/
            inline double       g   = 1;                         /*!< Transverse field strength */
            inline double       w   = 0;                         /*!< Randomness strength for the random field */
            inline int          d   = 2;                         /*!< Local dimension */
        }

        //Parameters for the selfdual transvese-field random-field next-neighbor Ising model
        namespace selfdual_tf_rf_ising {
            inline double       J_log_mean    = 0;               /*!< Average ferromagnetic coupling strength.*/
            inline double       h_log_mean    = 0;               /*!< Average transverse magnetic field strength */
            inline double       J_sigma       = 1;               /*!< Standard deviation for the lognormal distribution, i.e. = std(log(J)) , for the ferromagnetic coupling */
            inline double       h_sigma       = 0;               /*!< Standard deviation for the lognormal distribution, i.e. = std(log(h))   for the transverse magnetic field */
            inline double       lambda        = 0;               /*!< Lambda parameter */
            inline int          d             = 2;               /*!< Local dimension */
        }
    }

    //Parmaters that control MPS, eigensolver and SVD precision
    namespace precision {
        inline int      eigMaxIter                   = 1000  ;   /*!< Maximum number of steps for eigenvalue solver. */
        inline double   eigThreshold                 = 1e-12 ;   /*!< Minimum threshold for halting eigenvalue solver. */
        inline int      eigMaxNcv                    = 16    ;   /*!< Parameter controlling the column space? of the Lanczos solver. */
        inline double   SVDThreshold                 = 1e-8  ;   /*!< Minimum threshold value for keeping singular values. */
        inline double   VarConvergenceThreshold      = 1e-8  ;   /*!< Variance convergence threshold. The MPS state is considered good enough when its variance reaches below this value */
        inline double   VarSaturationThreshold       = 1e-4  ;   /*!< Variance saturation  threshold. The variance has saturated when its (absolute) slope reaches below this value */
        inline double   EntEntrSaturationThreshold   = 1e-4  ;   /*!< Entanglement Entropy saturation threshold. The entanglement entropy has saturated when its (absolute) slope reaches below this value*/
        inline int      MaxSizeFullDiag              = 2048  ;   /*!< Maximum linear size allowed for full diagonalization of the local hamiltonian matrix. Use 0 to allow any value */
    }

    //Parameters controlling iDMRG
    namespace idmrg {
        inline bool on         = true;                           /*!< Turns iDMRG simulation on/off. */
        inline int  max_steps  = 5000;                           /*!< Final length of 1D quantum chain. */
        inline long chi_max    = 8;                              /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool chi_grow   = true;                           /*!< Whether to increase chi slowly up to chi_max or go up to chi_max directly. */
        inline int  print_freq = 1000;                           /*!< Print frequency for console output. (0 = off). */
        inline int  store_freq = 100;                            /*!< Store frequency,for output file buffer. (0 = off). */

    }
    //Parameters controlling fDMRG
    namespace fdmrg {
        inline bool on           = true;                         /*!< Turns fDMRG simulation on/off. */
        inline int  num_sites    = 30;                           /*!< Number sweeps along the 1D quantum chain. */
        inline int  max_sweeps   = 10;                           /*!< Max number sweeps along the 1D quantum chain. */
        inline int  min_sweeps   = 4;                            /*!< Min number sweeps along the 1D quantum chain. */
        inline long chi_max      = 8;                            /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool chi_grow     = true;                         /*!< Whether to increase chi slowly up to chi_max or go up to chi_max directly. */
        inline int  print_freq   = 100;                          /*!< Print frequency for console output. In units of sweeps. (0 = off). */
        inline int  store_freq   = 100;                          /*!< Store frequency,for output file buffer. In units of sweeps. (0 = off). */
        inline bool store_wavefn = false;                        /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */

    }

    //Parameters controlling xDMRG
    namespace xdmrg {
        inline bool    on             = true;                    /*!< Turns xDMRG simulation on/off. */
        inline int     num_sites      = 200;                     /*!< Number sweeps along the 1D quantum chain. */
        inline int     max_sweeps     = 10;                      /*!< Number sweeps along the 1D quantum chain. */
        inline int     min_sweeps     = 4;                       /*!< Min number sweeps along the 1D quantum chain. */
        inline long    chi_max        = 8;                       /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool    chi_grow       = true;                    /*!< Whether to increase chi slowly up to chi_max or go up to chi_max directly. */
        inline int     print_freq     = 100;                     /*!< Print frequency for console output. In units of sweeps. (0 = off). */
        inline int     store_freq     = 100;                     /*!< Store frequency,for output file buffer. In units of sweeps. (0 = off). */
        inline bool    store_wavefn   = false;                   /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */
        inline double  energy_density = 0.5;                     /*!< Target energy in [0-1], where 0.5 means middle of spectrum. */
        inline double  energy_window  = 0.01;                    /*!< Accept states inside of energy_target +- energy_window. */

    }

    //Parameters controlling iTEBD
    namespace itebd {
        inline bool     on           = true;                     /*!< Turns iTEBD simulation on/off. */
        inline int      max_steps    = 100000;                   /*!< Number of iTEBD iterations, after which the simulation terminates regardless of convergence. Set high.*/
        inline double   delta_t0     = 0.1;                      /*!< Initial time step for iTEBD time evolution.*/
        inline double   delta_tmin   = 0.00001;                  /*!< Final time step for iTEBD time evolution.*/
        inline int      suzuki_order = 1;                        /*!< Order of the suzuki trotter decomposition (1,2 or 4) */
        inline long     chi_max      = 8;                        /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
        inline bool     chi_grow     = true;                     /*!< Whether to increase chi slowly up to chi_max or go up to chi_max directly. */
        inline int      print_freq   = 5000;                     /*!< Print frequency for console output. (0 = off).*/
        inline int      store_freq   = 100;                      /*!< Store frequency,for output file buffer. (0 = off). */

    }

    namespace hdf5 {
        inline bool         save_progress        = true;                         /*!< If true, saves the simulation data periodically */
        inline std::string  access_mode          = "output/default.h5";          /*!< Choose access mode to the file. Choose between READWRITE, READONLY */
        inline std::string  create_mode          = "READWRITE";                  /*!< Choose access mode to the file. Choose between TRUNCATE, OPEN, RENAME */
        inline std::string  output_filename      = "RENAME" ;                    /*!< Name of the output HDF5 file relative to the execution point  */
        inline StorageLevel storage_level        = StorageLevel::NORMAL;         /*!< Sets the storage level: choose "0=NONE,1=LIGHT,2=NORMAL,3=FULL */
        inline bool         store_profiling      = true;                         /*!< Whether to store profiling information to file. */
    }
    //Profiling
    namespace profiling {
        inline bool     on        = false;             /*!< If true, turns on profiling and timings will be shown on console. */
        inline int      precision = 5;                 /*!< Sets precision (number of decimals) of time output. */
    }
    //Console settings
    namespace console {
        inline int  verbosity  = 2;                    /*!< Level of verbosity desired [0-2]. Level 0 prints almost nothing, level 2 prints everything */
        inline bool timestamp  = false;                /*!< Whether to put a timestamp on console outputs */
    }
}
#endif //DMRG_N_SETTINGS_H
