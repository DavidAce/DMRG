//
// Created by david on 7/22/17.
//

#ifndef DMRG_CLASS_SUPERBLOCK_H
#define DMRG_CLASS_SUPERBLOCK_H

#include <general/nmspc_tensor_extra.h>
#include <memory>
#include <general/class_tic_toc.h>
#include <general/nmspc_eigsolver_props.h>

class class_mps_2site;
class class_mpo;
class class_hamiltonian;
class class_environment;
class class_environment_var;
template<typename Scalar> class class_SVD;

/*!
  \class class_superblock
  \brief This class contains the Hamiltonian MPO, current wave function MPS, left and right environment blocks and routines to contract, diagonalize, truncate
   and update them.
*/


class class_superblock {
public:
    using Scalar = std::complex<double>;
public:

    class_superblock();
    std::unique_ptr<class_mps_2site>         MPS;        /*!< Matrix product states for two sites, A and B, in Vidal Canonical Form \f$\Gamma^A\Lambda^A\Gamma^B\Lambda^B\f$. */
    std::unique_ptr<class_mpo>               H;
    std::unique_ptr<class_hamiltonian>       HA;
    std::unique_ptr<class_hamiltonian>       HB;
    std::unique_ptr<class_environment>       Lblock;     /*!< Left  environment block. */
    std::unique_ptr<class_environment>       Rblock;     /*!< Right environment block. */
    std::unique_ptr<class_environment_var>   Lblock2;    /*!< Left  environment block used for variance calculation */
    std::unique_ptr<class_environment_var>   Rblock2;    /*!< Right environment block used for variance calculation */
    std::unique_ptr<class_SVD<Scalar>>       SVD;
    //Bond dimensions and tensor shapes needed by the eigensolver and SVD.
    double E_optimal;                                              /*!< Stores the energy obtained in the eigenvalue solver. This energy corresponds to non-truncated MPS, so it will differ a tiny bit from what you see in final resuls. */

    class_tic_toc t_eig;

    unsigned long    environment_size = 0;

    Eigen::Tensor<Scalar, 4>
    optimize_MPS(Eigen::Tensor<Scalar, 4> &theta, Ritz ritz = Ritz::SR
    )    __attribute((hot));                            /*!< Finds the smallest algebraic eigenvalue and eigenvector (the ground state) using [Spectra](https://github.com/yixuan/spectra). */


    Eigen::Tensor<Scalar, 4>
    evolve_MPS(const Eigen::Tensor<Scalar, 4> &U);
    Eigen::Tensor<Scalar, 4>
    evolve_MPS(const Eigen::Tensor<Scalar, 4> &theta, const Eigen::Tensor<Scalar, 4> &U);
    Eigen::Tensor<Scalar,4> truncate_MPS(
            const Eigen::Tensor<Scalar, 4> &theta,        /*!< The 2-site MPS to truncate */
            long chi,                                      /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
            double SVDThreshold                            /*!< Minimum threshold value for keeping singular values. */
    )             __attribute((hot));                      /*!< Singular value decomposition, SVD, or Schmidt decomposition, of the ground state, where the truncation keeps \f$\chi\f$ (`chi`) singular values. */

    void truncate_MPS(
            const Eigen::Tensor<Scalar, 4> &theta,        /*!< The 2-site MPS to truncate */
            const std::shared_ptr<class_mps_2site> &MPS_out,
            long chi,                                      /*!< Bond dimension of the current position (maximum number of singular values to keep in SVD). */
            double SVDThreshold                            /*!< Minimum threshold value for keeping singular values. */
    )             __attribute((hot));                      /*!< Singular value decomposition, SVD, or Schmidt decomposition, of the ground state, where the truncation keeps \f$\chi\f$ (`chi`) singular values. */



   void enlarge_environment(int direction = 0);          /*!< Contract the MPS of current position \f$n\f$ into the left and right environments \f$L\f$ and \f$R\f$, i.e. `Lblock` and `Rblock`.
                                                         * \f[ L \leftarrow L \Lambda^B_{n-1} \Gamma^A_n W \Lambda^B_{n-1} (\Gamma^A_n)^* \f]
                                                         * \f[ R \leftarrow R \Gamma^B_{n+1} \Lambda^B_{n+1} W (\Gamma^B_{n+1})^* \Lambda^B_{n+1} \f] */

//    void set_current_dimensions()      ;                /*!< Update variables for dimensions */
    void swap_AB();                                     /*!< Swap the roles of A and B. Used in the infinite-DMRG stage.*/
};


#endif //DMRG_CLASS_SUPERBLOCK_H
