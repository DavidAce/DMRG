//
// Created by david on 2018-04-17.
//

#pragma once
#include <general/nmspc_tensor_omp.h>
#include <complex>
#include <vector>
#include <list>

namespace qm{
    using Scalar = std::complex<double>;
    extern std::vector<Eigen::MatrixXcd> gen_manybody_spin(const Eigen::MatrixXcd &s, int sites);
    namespace spinOneHalf {
        extern Eigen::Matrix2cd sx;
        extern Eigen::Matrix2cd sy;
        extern Eigen::Matrix2cd sz;
        extern Eigen::Matrix2cd Id;
        extern std::array<Eigen::Vector2cd,2> sx_eigvecs;
        extern std::array<Eigen::Vector2cd,2> sy_eigvecs;
        extern std::array<Eigen::Vector2cd,2> sz_eigvecs;
        extern std::vector<Eigen::MatrixXcd> SX;
        extern std::vector<Eigen::MatrixXcd> SY;
        extern std::vector<Eigen::MatrixXcd> SZ;
        extern std::vector<Eigen::MatrixXcd> II;

    }

    namespace spinOne{
        extern Eigen::Matrix3cd sx;
        extern Eigen::Matrix3cd sy;
        extern Eigen::Matrix3cd sz;
        extern Eigen::Matrix3cd Id;
        extern std::vector<Eigen::MatrixXcd> SX;
        extern std::vector<Eigen::MatrixXcd> SY;
        extern std::vector<Eigen::MatrixXcd> SZ;
        extern std::vector<Eigen::MatrixXcd> II;
    }

    namespace timeEvolution{
        extern std::vector<Eigen::MatrixXcd> Suzuki_Trotter_1st_order(const std::complex<double> t, const Eigen::MatrixXcd &h_evn, const Eigen::MatrixXcd &h_odd);
        extern std::vector<Eigen::MatrixXcd> Suzuki_Trotter_2nd_order(const std::complex<double> t, const Eigen::MatrixXcd &h_evn, const Eigen::MatrixXcd &h_odd);
        extern std::vector<Eigen::MatrixXcd> Suzuki_Trotter_4th_order(const std::complex<double> t, const Eigen::MatrixXcd &h_evn, const Eigen::MatrixXcd &h_odd);
        extern std::vector<Eigen::Tensor<std::complex<double>,4>> get_2site_evolution_gates(const std::complex<double> t, const int susuki_trotter_order,  const Eigen::MatrixXcd &h_evn, const Eigen::MatrixXcd &h_odd);
        extern std::vector<Eigen::Tensor<std::complex<double>,4>> compute_G(const std::complex<double> a, const int susuki_trotter_order, const Eigen::MatrixXcd &h_evn, const Eigen::MatrixXcd &h_odd);
    }

    namespace mpo{
        extern std::tuple<
                Eigen::Tensor<Scalar,4>,
                Eigen::Tensor<Scalar,3>,
                Eigen::Tensor<Scalar,3>>
        pauli_mpo(const Eigen::MatrixXcd paulimatrix);


        extern std::tuple<
                Eigen::Tensor<Scalar,4>,
                Eigen::Tensor<Scalar,3>,
                Eigen::Tensor<Scalar,3>>
        parity_selector_mpo(const Eigen::MatrixXcd paulimatrix, const int sector = 1);


        extern std::tuple<
                std::list<Eigen::Tensor<Scalar,4>>,
                Eigen::Tensor<Scalar,3>,
                Eigen::Tensor<Scalar,3>>
        parity_projector_mpos(const Eigen::MatrixXcd paulimatrix, const size_t sites, const int sector = 1);




    }



}

