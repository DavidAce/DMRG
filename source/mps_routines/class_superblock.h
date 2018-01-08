//
// Created by david on 7/22/17.
//

#ifndef DMRG_CLASS_SUPERBLOCK_H
#define DMRG_CLASS_SUPERBLOCK_H
//#define EIGEN_USE_MKL_ALL

#include "general/n_tensor_extra.h"
#include "mps_routines/class_environment.h"
#include "mps_routines/class_mps.h"
#include "mps_routines/class_mpo.h"
#include "sim_parameters/n_model.h"
#include "class_eigensolver_product.h"

/*!
 # Superblock Class
   This class contains the Hamiltonian MPO, current wave function MPS, left and right environment blocks and routines to contract, diagonalize, truncate
   and update them.
*/


class class_superblock {
public:
    using Scalar = class_mps::Scalar;
private:
    //Bond dimensions and tensor shapes needed by the eigensolver and SVD.

    long d;                                             /*!< Local, or physical, dimension of each particle. */
    long chiL;                                          /*!< Bond dimension of the previous (left) position. */
    long chiR;                                          /*!< Bond dimension of the next (right) position. */

    //Store temporaries for eigensolver and SVD.
    Tensor<Scalar,2> ground_state;                               /*!< Stores the ground state eigenvector from eigenvalue solver */
    Tensor<Scalar,3> U,V;                                        /*!< Stores the left and right unitary matrices \f$U\f$ and \f$V\f$ after an SVD decomposition \f$A = USV^\dagger\f$.*/
public:

    class_superblock();
    class_environment<Side::L>  Lblock;                 /*!< Left  environment block. */
    class_environment<Side::R>  Rblock;                 /*!< Right environment block. */
    class_environment_var<Side::L>  Lblock2;                 /*!< Left  environment block used for variance calculation */
    class_environment_var<Side::R>  Rblock2;                 /*!< Right environment block used for variance calculation */
    class_mps                   MPS;                    /*!< Matrix product states for two sites, A and B, in Vidal Canonical Form \f$\Gamma^A\Lambda^A\Gamma^B\Lambda^B\f$. */
    class_mpo           H;

    array1              shape1;                 /*!< Shape for Tensor1 representation of the MPS. */
    array2              shape2;                 /*!< Shape for Tensor2 representation of the MPS. */
    array4              shape4;                 /*!< Shape for Tensor4 representation of the MPS. */


    long   chi;                                           /*!< Bond dimension of the current (center) position. */
    long   chi_max ;                                      /*!< Maximum allowed bond dimension, which also defines the SVD cutoff */
    double truncation_error;
    int    chain_length = H.mps_sites;

    void find_ground_state(int eigSteps,                /*!< Maximum number of steps for eigenvalue solver. */
                           double eigThreshold          /*!< Minimum threshold for halting eigenvalue solver. */
                            )    __attribute((hot));    /*!< Finds the smallest algebraic eigenvalue and eigenvector (the ground state) using [Spectra](https://github.com/yixuan/spectra). */
    void time_evolve();
    void truncate(long chi,                             /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
                  double SVDThreshold                   /*!< Minimum threshold value for keeping singular values. */
                    )             __attribute((hot));   /*!< Singular value decomposition, SVD, or Schmidt decomposition, of the ground state, where the truncation keeps \f$\chi\f$ (`chi`) singular values. */
    void update_MPS();                                  /*!< Update `MPS.G.A` and `MPS.G.B` i.e., \f$\Gamma^A\f$ and \f$\Gamma^B\f$.*/
    void enlarge_environment(int direction = 0);        /*!< Contract the MPS of current position \f$n\f$ into the left and right environments \f$L\f$ and \f$R\f$, i.e. `Lblock` and `Rblock`.
                                                         * \f[ L \leftarrow L \Lambda^B_{n-1} \Gamma^A_n W \Lambda^B_{n-1} (\Gamma^A_n)^* \f]
                                                         * \f[ R \leftarrow R \Gamma^B_{n+1} \Lambda^B_{n+1} W (\Gamma^B_{n+1})^* \Lambda^B_{n+1} \f] */

    void initialize();
    void canonicalize_iMPS();
    void canonicalize_iMPS_iterative();
    void update_bond_dimensions();                      /*!< Update the MPS by contracting
                                                        * \f[ \Gamma^A \leftarrow (\Lambda^B_{n-1})^{-1} U \f]
                                                        * \f[ \Gamma^B \leftarrow V (\Lambda^B_{n+1})^{-1} \f] */
    void swap_AB();                                     /*!< Swap the roles of A and B. Used in the infinite-DMRG stage.*/

//    void reset();


};


#endif //DMRG_CLASS_SUPERBLOCK_H
