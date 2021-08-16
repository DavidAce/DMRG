#include "../spin.h"
#include <Eigen/Core>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <math/linalg/matrix.h>

Eigen::MatrixXcd qm::spin::gen_embedded_spin_operator(const Eigen::MatrixXcd &s, size_t at, size_t sites, bool mirror)
/*
 * Returns a spin operator embedded in a larger Hilbert space. For instance, if at == 1 and sites == 4:
 *
 *   σ¹ = i ⊗ σ ⊗ i ⊗ i
 *
 * where each element is a dxd matrix, resulting in a d^4 * d^4 matrix.

 * Note that if this matrix is converted to a rank-8 tensor, the indexing goes like:
 *
 @verbatim
        3 2 1 0
        | | | |
       [  σ¹  ]
       | | | |
       7 6 5 4
 @endverbatim

 * whereas you would normally want left-to-right indexing in MPS contexts:
 *
 @verbatim
        0 1 2 3
        | | | |
       [  σ¹  ]
       | | | |
       4 5 6 7
 @endverbatim

 * So don't forget to set "reverse = true" if you intend to use the result as a tensor.
 */

{
    if(at >= sites) throw std::logic_error(fmt::format("Expected at < sites. Got: at = {} | sites {}", at, sites));
    Eigen::MatrixXcd id     = Eigen::MatrixXcd::Identity(s.rows(), s.cols());
    Eigen::MatrixXcd result = at == 0 ? s : id;
    for(size_t site = 1; site < sites; site++)
        result = linalg::matrix::kronecker(result, site == at ? s : id, mirror).eval(); // .eval() is required to avoid aliasing!!
    return result;
}

std::vector<Eigen::MatrixXcd> qm::spin::gen_manybody_spins(const Eigen::MatrixXcd &s, size_t sites, bool reverse) {
    std::vector<Eigen::MatrixXcd> S;
    S.reserve(sites);
    for(size_t site = 0; site < sites; site++) S.emplace_back(gen_embedded_spin_operator(s, site, sites, reverse));
    return S;
}