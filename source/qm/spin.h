#pragma once
#include "math/tenx/fwd_decl.h"
#include "qm.h"
#include <array>
#include <complex>
#include <vector>

namespace qm::spin {
    extern Eigen::MatrixXcd              gen_embedded_spin_operator(const Eigen::MatrixXcd &s, size_t at, size_t site, bool reverse = false);
    extern std::vector<Eigen::MatrixXcd> gen_manybody_spins(const Eigen::MatrixXcd &s, size_t sites, bool reverse = false);

    namespace half {
        extern Eigen::Matrix2cd                sx;
        extern Eigen::Matrix2cd                sy;
        extern Eigen::Matrix2cd                sz;
        extern Eigen::Matrix2cd                sp;
        extern Eigen::Matrix2cd                sm;
        extern Eigen::Matrix2cd                nu;
        extern Eigen::Matrix2cd                nd;
        extern Eigen::Matrix2cd                id;
        extern std::array<Eigen::Vector2cd, 2> sx_spinors;
        extern std::array<Eigen::Vector2cd, 2> sy_spinors;
        extern std::array<Eigen::Vector2cd, 2> sz_spinors;
        extern std::vector<Eigen::MatrixXcd>   SX;
        extern std::vector<Eigen::MatrixXcd>   SY;
        extern std::vector<Eigen::MatrixXcd>   SZ;
        extern std::vector<Eigen::MatrixXcd>   II;
        extern Eigen::MatrixXcd                gen_embedded_spin_half_operator(const Eigen::Matrix2cd &s, size_t at, size_t sites, bool swap = false);
        extern std::vector<Eigen::Matrix4cd>   gen_twobody_spins(const Eigen::Matrix2cd &s, bool swap = false);

        inline static constexpr std::array<std::string_view, 9> valid_axis_str = {"x", "+x", "-x", "y", "+y", "-y", "z", "+z", "-z"};
        extern bool                                             is_valid_axis(std::string_view sector);
        extern int                                              get_sign(std::string_view sector);
        extern std::string_view                                 get_axis(std::string_view sector);
        extern Eigen::Vector2cd                                 get_spinor(std::string_view axis, int sign);
        extern Eigen::Vector2cd                                 get_spinor(std::string_view sector);
        extern Eigen::Matrix2cd                                 get_pauli(std::string_view axis);
    }

    namespace qm::spin::one {
        extern Eigen::Matrix3cd              sx;
        extern Eigen::Matrix3cd              sy;
        extern Eigen::Matrix3cd              sz;
        extern Eigen::Matrix3cd              id;
        extern std::vector<Eigen::MatrixXcd> SX;
        extern std::vector<Eigen::MatrixXcd> SY;
        extern std::vector<Eigen::MatrixXcd> SZ;
        extern std::vector<Eigen::MatrixXcd> II;
        extern std::vector<Eigen::MatrixXcd> gen_twobody_spins(const Eigen::Matrix3cd &s, bool swap = false);
    }

}