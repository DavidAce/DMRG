//
// Created by david on 4/25/17.
//

#include "nmspc_model.h"
#include <iomanip>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>
#include <sim_parameters/nmspc_sim_settings.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_math.h>

using namespace std;
using namespace Eigen;
using namespace qm::spinOneHalf;
namespace Model {


    /*! From Zapp, K. Matrix Product States for Lattice Gauge Theories. (2015).
     * and
     * Müller-hermes, B. V. A. Tensor-Network-Methods for Simulating Infinite 1-dimensional Quantum-Many-Body Systems. 1–68 (2010).
     */

    double get_exact_energy() {
        using namespace settings::model::tf_ising;
        return -1.27;
//        return (-1.0 / M_PI / 2.0) *
//               Math::compute_integral([](double x) { return sqrt(1.0 + g * g - 2.0 * g * cos(x)); }, {-M_PI, M_PI});
    }


    MatrixXcd h(int sites, int position) {
        using namespace settings::model::tf_ising;
        int i = Math::mod(position, sites);
        int j = Math::mod(position + 1, sites);
        if (spins_must_be_generated) {
            SX = qm::gen_manybody_spin(sx, sites);
            SY = qm::gen_manybody_spin(sy, sites);
            SZ = qm::gen_manybody_spin(sz, sites);
            spins_must_be_generated = false;
        }
        return (-J * SZ[i] * SZ[j] - 0.5 * g * (SX[i] + SX[j]));
    }

    MatrixXcd H(int sites) {
        MatrixXcd hi = MatrixXcd::Zero((long) pow(2, sites), (long) pow(2, sites));
        for (int position = 0; position < sites; position++) {
            hi += h(sites, position);
        }
        return hi;
    }

//    MatrixXcd Hsq(int sites) {
//        MatrixXcd hi = MatrixXcd::Zero((long) pow(2, sites), (long) pow(2, sites));
////        hi = h(sites, 0) * h(sites, 0) +  h(sites, 1) * h(sites, 1) + h(sites, 0) * h(sites, 1) + h(sites, 1) * h(sites, 0);
//        for (int position1 = 0; position1 < sites; position1++) {
//            for(int position2 = 0; position2 < sites; position2++){
//                hi += h(sites, position1) * h(sites, position2);
//            }
//        }
//        return hi;
//    }

//    MatrixXcd H_MPO(double e) {
//        /*! Returns the MPO as a matrix. Notation following Schollwöck (2010) */
//        using namespace settings::model::tf_ising;
//        MatrixXcd W(6, 6);
//        W.setZero();
//        W.block(0, 0, 2, 2) = I;
//        W.block(2, 0, 2, 2) = sz;
//        W.block(4, 0, 2, 2) = -g * sx - e * I; // Optionally subtract a constant. Default is k = 0.
//        W.block(4, 2, 2, 2) = -J * sz;
//        W.block(4, 4, 2, 2) = I;
//        return W;
//    }
//
//    MatrixXcd H_MPO_random_field(double g, double e) {
//        /*! Returns the MPO as a matrix. Notation following Schollwöck (2010) */
//        using namespace settings::model::tf_ising;
//        MatrixXcd W(6, 6);
//        W.setZero();
//        W.block(0, 0, 2, 2) = I;
//        W.block(2, 0, 2, 2) = sz;
//        W.block(4, 0, 2, 2) = -g * sz - e*I; // Optionally subtract a constant. Default is k = 0.
//        W.block(4, 2, 2, 2) = -J * sz;
//        W.block(4, 4, 2, 2) = I;
//        return W;
//    }
}