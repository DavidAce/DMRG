//
// Created by david on 2018-05-04.
//

#ifndef DMRG_CLASS_OPTIMIZE_MPS_H
#define DMRG_CLASS_OPTIMIZE_MPS_H
#include <iomanip>
#include <complex>
#include <vector>

template<class T>
class class_contraction {

private:
    using Scalar = std::complex<double>;
    const Scalar *Lblock;
    const Scalar *Rblock;
    const Scalar *HA;
    const Scalar *HB;
    std::array<long,4> shape_theta4;
    std::array<long,2> shape_theta2;
    std::array<long,1> shape_theta1;
    std::array<long,4> shape_mpo4;
    Scalar *testmat;
public:
    int rows()const {return (int)shape_theta1[0];};               /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */
    int cols()const {return (int)shape_theta1[0];};               /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */

    void MultMv(T* theta_in_, T* theta_out_);               /*!< The function that contracts.  */
    int counter = 0;
    class_contraction(
            const Scalar *Lblock_,        /*!< The left block tensor.  */
            const Scalar *Rblock_,        /*!< The right block tensor.  */
            const Scalar *HA_,            /*!< The left Hamiltonian MPO's  */
            const Scalar *HB_,            /*!< The right Hamiltonian MPO's */
            const std::array<long,4> shape_theta4_,         /*!< An array containing the shapes of theta  */
            const std::array<long,4> shape_mpo4_,            /*!< An array containing the shapes of the MPO  */
                  Scalar *testmat_ = nullptr
    ):                                                    /*!< Initializes the custom contraction. */
            Lblock(Lblock_),
            Rblock(Rblock_),
            HA(HA_),
            HB(HB_),
            shape_theta4(shape_theta4_),
            shape_theta2({shape_theta4[0] * shape_theta4[1] , shape_theta4[2] * shape_theta4[3]}),
            shape_theta1({shape_theta4[0] * shape_theta4[1] * shape_theta4[2] * shape_theta4[3]}),
            shape_mpo4(shape_mpo4_),
            testmat(testmat_)
    {
    }
};


class class_optimize_mps{
private:
    using Scalar = std::complex<double>;
    const Scalar *Lblock;
    const Scalar *Rblock;
    const Scalar *HA;
    const Scalar *HB;
    std::array<long,4> shape_theta4;
    std::array<long,2> shape_theta2;
    std::array<long,1> shape_theta1;
    std::array<long,4> shape_mpo4;
    Scalar *resid;
    Scalar *testmat;


public:
    std::vector<Scalar> eigvecs;
    std::vector<Scalar> eigvals;
    int rows = 0;
    int cols = 0;
    int iter = 0;
    int counter = 0;

    class_optimize_mps (            const Scalar *Lblock_,        /*!< The left block tensor.  */
                                    const Scalar *Rblock_,        /*!< The right block tensor.  */
                                    const Scalar *HA_,            /*!< The left Hamiltonian MPO's  */
                                    const Scalar *HB_,            /*!< The right Hamiltonian MPO's */
                                    const std::array<long,4> shape_theta4_,         /*!< An array containing the shapes of theta  */
                                    const std::array<long,4> shape_mpo4_ ,           /*!< An array containing the shapes of the MPO  */
                                          Scalar *resid_,
                                          Scalar *testmat_= nullptr
    ):                                                    /*!< Initializes the custom contraction. */
            Lblock(Lblock_),
            Rblock(Rblock_),
            HA(HA_),
            HB(HB_),
            shape_theta4(shape_theta4_),
            shape_theta2({shape_theta4[0] * shape_theta4[1] , shape_theta4[2] * shape_theta4[3]}),
            shape_theta1({shape_theta4[0] * shape_theta4[1] * shape_theta4[2] * shape_theta4[3]}),
            shape_mpo4(shape_mpo4_),
            resid(resid_),
            testmat(testmat_)
    {

    }

    void optimize_mps(int nev = 1);

};


#endif //DMRG_CLASS_OPTIMIZE_MPS_H
