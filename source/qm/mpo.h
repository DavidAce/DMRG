#pragma once

#include "../config/enums.h"
#include "math/tenx/fwd_decl.h"
#include "qm.h"
#include <vector>

namespace qm::mpo {
    extern std::tuple<Eigen::Tensor<cplx, 4>, Eigen::Tensor<cplx, 3>, Eigen::Tensor<cplx, 3>> pauli_mpo(const Eigen::MatrixXcd &paulimatrix);
    extern std::tuple<Eigen::Tensor<cplx, 4>, Eigen::Tensor<cplx, 3>, Eigen::Tensor<cplx, 3>> pauli_mpo(std::string_view axis);

    extern std::tuple<std::vector<Eigen::Tensor<cplx, 4>>, Eigen::Tensor<cplx, 3>, Eigen::Tensor<cplx, 3>>
        parity_projector_mpos(const Eigen::MatrixXcd &paulimatrix, size_t sites, int sector = 1);

    extern std::tuple<std::vector<Eigen::Tensor<cplx, 4>>, Eigen::Tensor<cplx, 3>, Eigen::Tensor<cplx, 3>>
        random_pauli_mpos(const Eigen::MatrixXcd &paulimatrix, size_t sites);

    extern std::tuple<std::vector<Eigen::Tensor<cplx, 4>>, Eigen::Tensor<cplx, 3>, Eigen::Tensor<cplx, 3>>
        random_pauli_mpos_x2(const Eigen::MatrixXcd &paulimatrix1, const Eigen::MatrixXcd &paulimatrix2, size_t sites);

    extern std::tuple<std::vector<Eigen::Tensor<cplx, 4>>, Eigen::Tensor<cplx, 3>, Eigen::Tensor<cplx, 3>>
        random_pauli_mpos(const std::vector<Eigen::Matrix2cd> &paulimatrices, size_t sites);

    extern std::tuple<std::vector<Eigen::Tensor<cplx, 4>>, Eigen::Tensor<cplx, 3>, Eigen::Tensor<cplx, 3>>
        sum_of_pauli_mpo(const std::vector<Eigen::Matrix2cd> &paulimatrices, size_t sites, RandomizerMode mode);

    extern std::tuple<std::vector<Eigen::Tensor<cplx, 4>>, Eigen::Tensor<cplx, 3>, Eigen::Tensor<cplx, 3>>
        random_pauli_mpos(const std::vector<Eigen::Matrix2cd> &paulimatrices, const std::vector<double> &uniform_dist_widths, size_t sites);
}