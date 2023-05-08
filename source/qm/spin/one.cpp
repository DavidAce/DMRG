#include "../spin.h"
#include <Eigen/Core>

namespace qm::spin::one {
    /* clang-format off */

    Eigen::Matrix3cd sx = (Eigen::Matrix3cd() <<  0.0, 1.0, 0.0,
        1.0, 0.0, 1.0,
        0.0, 1.0, 0.0).finished();
    Eigen::Matrix3cd sy = (Eigen::Matrix3cd() <<  0.0 , -1.0i, 0.0,
        1.0i,  0.0 ,-1.0i,
        0.0 ,  1.0i, 0.0).finished();
    Eigen::Matrix3cd sz = (Eigen::Matrix3cd() <<  1.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0,-1.0).finished();
    Eigen::Matrix3cd id = (Eigen::Matrix3cd()  << 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0).finished();

    std::vector<Eigen::MatrixXcd> SX;
    std::vector<Eigen::MatrixXcd> SY;
    std::vector<Eigen::MatrixXcd> SZ;
    std::vector<Eigen::MatrixXcd> II;
    /* clang-format on */

    std::vector<Eigen::MatrixXcd> gen_twobody_spins(const Eigen::Matrix3cd &s, bool swap)
    // Returns a pair of two-body 4x4 spin operators for embedded in a two-site Hilbert space:
    //        (σ ⊗ i, i ⊗ σ)
    // where σ is a 3x3 (pauli) matrix and i is the 3x3 identity matrix.
    // So don't forget to set "swap = true" if you intend to use the result as a tensor.
    {
        return gen_manybody_spins(s, 2, swap);
    }
}