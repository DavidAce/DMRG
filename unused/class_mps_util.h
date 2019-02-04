//
// Created by david on 2018-07-06.
//

#ifndef DMRG_CLASS_MPS_UTIL_H
#define DMRG_CLASS_MPS_UTIL_H
#include "general/nmspc_tensor_extra.h"

class class_mps_2site;


class class_mps_util {
    using Scalar = std::complex<double>;
public:

    class_mps_util()=default;
    class_mps_util(const std::unique_ptr <class_mps_2site> &MPS){
        compute_mps_components(MPS);
    };

    void compute_mps_components(const std::unique_ptr <class_mps_2site> &MPS);


    Eigen::Tensor<Scalar,4> theta, theta_evn_normalized, theta_odd_normalized;
    Eigen::Tensor<Scalar,4> theta_sw ;
    Eigen::Tensor<Scalar,3> LBGA, LAGB;
    Eigen::Tensor<Scalar,2> l_evn, r_evn;
    Eigen::Tensor<Scalar,2> l_odd, r_odd;
    Eigen::Tensor<Scalar,4> transfer_matrix_LBGA;
    Eigen::Tensor<Scalar,4> transfer_matrix_LAGB;
    Eigen::Tensor<Scalar,4> transfer_matrix_evn;
    Eigen::Tensor<Scalar,4> transfer_matrix_odd;

    Eigen::Tensor<Scalar,4> get_theta                       (const std::unique_ptr<class_mps_2site> &MPS, Scalar norm = 1.0)  const;              /*!< Returns rank 4 tensor \f$\Theta\f$.*/
    Eigen::Tensor<Scalar,4> get_theta_swapped               (const std::unique_ptr<class_mps_2site> &MPS, Scalar norm = 1.0)  const;              /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/
    Eigen::Tensor<Scalar,4> get_theta_evn                   (const std::unique_ptr<class_mps_2site> &MPS, Scalar norm = 1.0)  const;              /*!< Returns rank 4 tensor \f$\Theta\f$.*/
    Eigen::Tensor<Scalar,4> get_theta_odd                   (const std::unique_ptr<class_mps_2site> &MPS, Scalar norm = 1.0)  const;              /*!< Returns rank 4 tensor \f$\Theta\f$, with A and B swapped.*/

    Eigen::Tensor<Scalar,4> get_transfer_matrix_zero        (const std::unique_ptr<class_mps_2site> &MPS) const;
    Eigen::Tensor<Scalar,4> get_transfer_matrix_LBGA        (const std::unique_ptr <class_mps_2site> &MPS, Scalar norm = 1.0) const;
    Eigen::Tensor<Scalar,4> get_transfer_matrix_GALC        (const std::unique_ptr <class_mps_2site> &MPS, Scalar norm = 1.0) const;
    Eigen::Tensor<Scalar,4> get_transfer_matrix_GBLB        (const std::unique_ptr <class_mps_2site> &MPS, Scalar norm = 1.0) const;
    Eigen::Tensor<Scalar,4> get_transfer_matrix_LCGB        (const std::unique_ptr <class_mps_2site> &MPS, Scalar norm = 1.0) const;
    Eigen::Tensor<Scalar,4> get_transfer_matrix_theta_evn   (const std::unique_ptr <class_mps_2site> &MPS, Scalar norm = 1.0) const;
    Eigen::Tensor<Scalar,4> get_transfer_matrix_theta_odd   (const std::unique_ptr <class_mps_2site> &MPS, Scalar norm = 1.0) const;
    Eigen::Tensor<Scalar,4> get_transfer_matrix_AB          (const std::unique_ptr <class_mps_2site> &MPS, int p)             const;
};


#endif //DMRG_CLASS_MPS_UTIL_H
