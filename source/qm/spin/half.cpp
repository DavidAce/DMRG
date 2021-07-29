#include "../spin.h"
#include <array>
#include <Eigen/Core>
#include <vector>

namespace qm::spin::half {

    /* clang-format off */
    Eigen::Matrix2cd sx = (Eigen::Matrix2cd() <<
        0.0, 1.0,
        1.0, 0.0).finished();
    Eigen::Matrix2cd sy = (Eigen::Matrix2cd() <<
        0.0, imn,
        imp, 0.0).finished();
    Eigen::Matrix2cd sz = (Eigen::Matrix2cd() <<
        1.0, 0.0,
        0.0, -1.0).finished();
    Eigen::Matrix2cd sp = (Eigen::Matrix2cd() <<
        0.0, 1.0,
        0.0, 0.0).finished();
    Eigen::Matrix2cd sm = (Eigen::Matrix2cd() <<
        0.0, 0.0,
        1.0, 0.0).finished();
    Eigen::Matrix2cd id  = (Eigen::Matrix2cd() << 1.0, 0.0,
        0.0, 1.0).finished();

    std::array<Eigen::Vector2cd,2> sx_spinors{(Eigen::Vector2cd() << 1.0, 1.0).finished()/std::sqrt(2),
                                       (Eigen::Vector2cd() << 1.0,-1.0).finished()/std::sqrt(2)};

    std::array<Eigen::Vector2cd,2> sy_spinors{(Eigen::Vector2cd() << 1.0, imp).finished()/std::sqrt(2),
                                       (Eigen::Vector2cd() << 1.0, imn).finished()/std::sqrt(2)};

    std::array<Eigen::Vector2cd,2> sz_spinors{(Eigen::Vector2cd() << 1.0, 0.0).finished()/std::sqrt(2),
                                       (Eigen::Vector2cd() << 0.0, 1.0).finished()/std::sqrt(2)};

    std::vector<Eigen::MatrixXcd> SX;
    std::vector<Eigen::MatrixXcd> SY;
    std::vector<Eigen::MatrixXcd> SZ;
    std::vector<Eigen::MatrixXcd> II;
    /* clang-format on */

    Eigen::MatrixXcd gen_embedded_spin_half_operator(const Eigen::Matrix2cd &s, size_t at, size_t sites, bool swap) {
        // Don't forget to set "swap = true" if you intend to use the result as a tensor.
        return gen_embedded_spin_operator(s, at, sites, swap);
    }

    std::vector<Eigen::Matrix4cd> gen_twobody_spins(const Eigen::Matrix2cd &s, bool swap)
    // Returns a pair of two-body 4x4 spin operators for embedded in a two-site Hilbert space:
    //        (σ ⊗ i, i ⊗ σ)
    // where σ is a 2x2 (pauli) matrix and i is the 2x2 identity matrix.
    // Don't forget to set "swap = true" if you intend to use the result as a tensor.
    {
        std::vector<Eigen::Matrix4cd> S;
        for(size_t site = 0; site < 2; site++) S.emplace_back(gen_embedded_spin_operator(s, site, 2, swap));
        return S;
    }

}