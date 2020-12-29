//
// Created by david on 2018-04-17.
//

#pragma once

#include <Eigen/Core>
#include <complex>
#include <general/eigen_tensor_fwd_decl.h>
#include <list>
#include <vector>
#include <physics/class_quantum_gates.h>
namespace qm{
    /* clang-format off */
    using Scalar = std::complex<double>;
    using cplx = std::complex<double>;
    using real = double;

    extern Eigen::MatrixXcd gen_embedded_spin_operator(const Eigen::MatrixXcd &s, size_t at, size_t site, bool swap = false);
    extern std::vector<Eigen::MatrixXcd> gen_manybody_spins(const Eigen::MatrixXcd &s, int sites, bool swap = false);

    constexpr std::complex<double> imp(0.0,1.0);
    constexpr std::complex<double> imn(0.0,-1.0);
    namespace spinHalf {
        extern Eigen::Matrix2cd sx;
        extern Eigen::Matrix2cd sy;
        extern Eigen::Matrix2cd sz;
        extern Eigen::Matrix2cd id;
        extern std::array<Eigen::Vector2cd,2> sx_spinors;
        extern std::array<Eigen::Vector2cd,2> sy_spinors;
        extern std::array<Eigen::Vector2cd,2> sz_spinors;
        extern std::vector<Eigen::MatrixXcd> SX;
        extern std::vector<Eigen::MatrixXcd> SY;
        extern std::vector<Eigen::MatrixXcd> SZ;
        extern std::vector<Eigen::MatrixXcd> II;

        extern Eigen::MatrixXcd gen_embedded_spin_operator(const Eigen::Matrix2cd &s, size_t at, size_t sites, bool swap = false);
        extern std::vector<Eigen::Matrix4cd> gen_twobody_spins(const Eigen::Matrix2cd &s, bool swap = false);
    }

    namespace spinOne{
        extern Eigen::Matrix3cd sx;
        extern Eigen::Matrix3cd sy;
        extern Eigen::Matrix3cd sz;
        extern Eigen::Matrix3cd id;
        extern std::vector<Eigen::MatrixXcd> SX;
        extern std::vector<Eigen::MatrixXcd> SY;
        extern std::vector<Eigen::MatrixXcd> SZ;
        extern std::vector<Eigen::MatrixXcd> II;
        extern std::vector<Eigen::MatrixXcd> gen_twobody_spins(const Eigen::Matrix3cd &s, bool swap = false);
    }

    namespace timeEvolution{
        extern std::vector<Eigen::Tensor<Scalar,2>> Suzuki_Trotter_1st_order(cplx delta_t, const Eigen::Tensor<Scalar,2> &h_evn, const Eigen::Tensor<Scalar,2> &h_odd);
        extern std::vector<Eigen::Tensor<Scalar,2>> Suzuki_Trotter_2nd_order(cplx delta_t, const Eigen::Tensor<Scalar,2> &h_evn, const Eigen::Tensor<Scalar,2> &h_odd);
        extern std::vector<Eigen::Tensor<Scalar,2>> Suzuki_Trotter_4th_order(cplx delta_t, const Eigen::Tensor<Scalar,2> &h_evn, const Eigen::Tensor<Scalar,2> &h_odd);
        extern std::vector<Eigen::Tensor<Scalar,2>> get_twosite_time_evolution_operators(cplx delta_t, size_t susuki_trotter_order, const Eigen::Tensor<Scalar,2> &h_evn, const Eigen::Tensor<Scalar,2> &h_odd);
        extern std::vector<Eigen::Tensor<Scalar,2>> compute_G(cplx a, size_t susuki_trotter_order, const Eigen::Tensor<Scalar,2> &h_evn, const Eigen::Tensor<Scalar,2> &h_odd);
    }

    namespace lbit{
        extern Eigen::Tensor<Scalar,2> get_time_evolution_operator(cplx delta_t, const Eigen::Tensor<Scalar,2> &hamiltonian);
        extern std::vector<qm::Gate> get_time_evolution_gates(cplx delta_t, const std::vector<qm::Gate> &hams_nsite);
        extern std::vector<qm::Gate> get_unitary_2gate_layer(size_t sites, double fmix);
        extern std::vector<Eigen::Tensor<Scalar,2>> get_time_evolution_operators_2site(size_t sites, cplx delta_t, const std::vector<Eigen::Tensor<Scalar,2>> &twosite_hams);
        extern std::vector<Eigen::Tensor<Scalar,2>> get_time_evolution_operators_3site(size_t sites, cplx delta_t, const std::vector<Eigen::Tensor<Scalar,2>> &hams_3site);
        extern std::vector<Eigen::Tensor<Scalar,4>> get_time_evolution_mpos(cplx delta_t, const std::vector<Eigen::Tensor<Scalar,4>> &mpos);

    }

    namespace mpo{
        extern std::tuple<
                Eigen::Tensor<Scalar,4>,
                Eigen::Tensor<Scalar,3>,
                Eigen::Tensor<Scalar,3>>
        pauli_mpo(const Eigen::MatrixXcd &paulimatrix);


        extern std::tuple<
                std::list<Eigen::Tensor<Scalar,4>>,
                Eigen::Tensor<Scalar,3>,
                Eigen::Tensor<Scalar,3>>
        parity_projector_mpos(const Eigen::MatrixXcd &paulimatrix, size_t sites, int sector = 1);

        extern std::tuple<
            std::list<Eigen::Tensor<Scalar,4>>,
            Eigen::Tensor<Scalar,3>,
            Eigen::Tensor<Scalar,3>>
        random_pauli_mpos(const Eigen::MatrixXcd &paulimatrix, size_t sites);

        extern std::tuple<
            std::list<Eigen::Tensor<Scalar,4>>,
            Eigen::Tensor<Scalar,3>,
            Eigen::Tensor<Scalar,3>>
        random_pauli_mpos_x2(const Eigen::MatrixXcd &paulimatrix1, const Eigen::MatrixXcd &paulimatrix2, size_t sites);

        extern std::tuple<
            std::list<Eigen::Tensor<Scalar,4>>,
            Eigen::Tensor<Scalar,3>,
            Eigen::Tensor<Scalar,3>>
        random_pauli_mpos(const std::vector<Eigen::Matrix2cd> &paulimatrices, size_t sites);

        extern std::tuple<
            std::list<Eigen::Tensor<Scalar,4>>,
            Eigen::Tensor<Scalar,3>,
            Eigen::Tensor<Scalar,3>>
        sum_of_pauli_mpo(const std::vector<Eigen::Matrix2cd> &paulimatrices, size_t sites, RandomizerMode mode);

        extern std::tuple<
            std::list<Eigen::Tensor<Scalar,4>>,
            Eigen::Tensor<Scalar,3>,
            Eigen::Tensor<Scalar,3>>
        random_pauli_mpos(const std::vector<Eigen::Matrix2cd> &paulimatrices, const std::vector<double> & uniform_dist_widths, size_t sites);
    }

    /* clang-format on */


}

