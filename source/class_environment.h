//
// Created by david on 7/21/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_ENVIRONMENT_H
#define FINITE_DMRG_EIGEN_CLASS_ENVIRONMENT_H

#include <n_tensor_extra.h>
#include <class_MPS.h>


using namespace Textra;
using namespace std;

/*! \brief Left environment block.
 *
 * # Left environment
 * This class contains the Left environment block as a rank-3 tensor. New sites are contracted into the left environment as
 * \f[
 * L \leftarrow L \Lambda^B_{n-1} \Gamma^A_n W \Lambda^B_{n-1} (\Gamma^A_n)^*
 * \f]
 *
 * where \f$W\f$ is the rank-4 tensor Hamiltonian MPO.
 */
class class_environment_L {
private:

public:
    std::string single_picture = "=";                       /*!< String representation each site in the environment (i.e. left of current MPS). */
    std::string picture;                                    /*!< Pretty string representation of the left environment block. */
    int size;                                               /*!< Number of particles that have been contracted into this left environment. */

    Tensor3 block;                                         /*!< The environment block. */
    class_environment_L();
    void enlarge(const class_MPS &MPS, const Tensor4 &W);  /*!< Contracts a site into the block. */
};




/*! \brief Right environment block.
 *
 * # Right environment
 * This class contains the Right environment block as a rank-3 tensor. New sites are contracted into the Right environment as
 * \f[
 * R \leftarrow R \Gamma^B_{n+1} \Lambda^B_{n+1} W (\Gamma^B_{n+1})^* \Lambda^B_{n+1}
 * \f]
 *
 * where \f$W\f$ is the rank-4 tensor Hamiltonian MPO.
 */

class class_environment_R {
private:

public:
    std::string single_picture = "-";                       /*!< String representation each site in the right environment (i.e. right of current MPS). */
    std::string picture;                                    /*!< Pretty string representation of the right environment block. */
    int size;                                               /*!< Number of particles that have been contracted into this right environment. */

    Tensor3 block;                                          /*!< The environment block. */
    class_environment_R();

    void enlarge(const class_MPS &MPS,const Tensor4 &W);    /*!< Contracts a site into the block. */
};


#endif //FINITE_DMRG_EIGEN_CLASS_ENVIRONMENT_H
