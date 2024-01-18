#include "../opt_mps.h"
#include "general/iter.h"
#include "tools/common/log.h"
#include "tools/finite/opt/opt-internal.h"

using namespace tools::finite::opt;
using namespace tools::finite::opt::internal;

void tools::finite::opt::internal::subspace::filter_subspace(std::vector<opt_mps> &subspace, size_t max_accept) {
    // Sort the eigvec list in order of descending overlaps. If the overlaps are the same, compare instead the distance in energy to the
    // current energy
    std::sort(subspace.begin(), subspace.end(), std::greater<>());
    size_t initial_size           = subspace.size();
    size_t min_accept             = std::min(std::min(max_accept, 32ul), initial_size);
    max_accept                    = std::min(max_accept, initial_size);
    auto   subspace_errors        = subspace::get_subspace_errors(subspace); // Vector with decreasing (cumulative) subspace errors, 1-eps(eigvec_i)
    double subspace_error_initial = subspace_errors.back();
    while(true) {
        if(subspace.size() <= max_accept) break;
        if(subspace.back().is_basis_vector) subspace_errors.pop_back();
        subspace.pop_back();
    }
    double subspace_error_final = subspace_errors.back();
    size_t idx                  = 0;
    for(auto &eigvec : subspace) {
        eigvec.set_name(fmt::format("eigenvector {}", idx++));
        //        tools::log->trace("Filtered {:<16}: overlap {:.16f} | energy {:>20.16f}", eigvec.get_label(), eigvec.get_overlap(),
        //                          eigvec.get_energy());
    }

    tools::log->trace("Filtered from {} down to {} states", initial_size, subspace.size());
    tools::log->trace("Filter changed subspace error = {:8.2e} --> {:8.2e}", subspace_error_initial, subspace_error_final);
    if(subspace.size() < min_accept) throw std::runtime_error("Filtered too many eigvecs");
}

std::optional<size_t> tools::finite::opt::internal::subspace::get_idx_to_eigvec_with_highest_overlap(const std::vector<opt_mps> &eigvecs, double energy_llim,
                                                                                                     double energy_ulim) {
    if(eigvecs.empty()) return std::nullopt;
    auto overlaps = get_overlaps(eigvecs);
    long idx      = 0;
    for(const auto &eigvec : eigvecs) {
        if(not eigvec.is_basis_vector) continue;
        auto energy = eigvec.get_energy();
        if(energy > energy_ulim) overlaps(idx) = 0.0;
        if(energy < energy_llim) overlaps(idx) = 0.0;
        idx++;
    }
    // Now we have a list of overlaps where nonzero elements correspond to eigvecs inside the energy window
    // Get the index to the highest overlapping element
    double max_overlap_val = overlaps.maxCoeff(&idx);
    if(max_overlap_val == 0.0) {
        tools::log->debug("No overlapping eigenstates in given energy window {} to {}.", energy_llim, energy_ulim);
        for(const auto &eigvec : eigvecs) tools::log->info("energy {}: {}", eigvec.get_eigs_idx(), eigvec.get_energy());
        return std::nullopt;
    }
    return idx;
}

std::optional<size_t> tools::finite::opt::internal::subspace::get_idx_to_eigvec_with_lowest_variance(const std::vector<opt_mps> &eigvecs, double energy_llim,
                                                                                                     double energy_ulim) {
    if(eigvecs.empty()) return std::nullopt;
    auto   var = std::numeric_limits<double>::infinity();
    size_t idx = 0;
    for(const auto &[i, eigvec] : iter::enumerate(eigvecs)) {
        if(not eigvec.is_basis_vector) continue;
        if(eigvec.get_variance() < var and eigvec.get_energy() <= energy_ulim and eigvec.get_energy() >= energy_llim) {
            idx = i;
            var = eigvec.get_variance();
        }
    }
    if(std::isinf(var)) return std::nullopt;
    return idx;
}

std::vector<size_t> tools::finite::opt::internal::subspace::get_idx_to_eigvec_with_highest_overlap(const std::vector<opt_mps> &eigvecs, size_t max_eigvecs,
                                                                                                   double energy_llim, double energy_ulim) {
    if(eigvecs.empty()) return std::vector<size_t>();
    auto overlaps = get_overlaps(eigvecs);
    long idx      = 0;
    for(const auto &eigvec : eigvecs) {
        if(not eigvec.is_basis_vector) continue;
        if(eigvec.get_energy() > energy_ulim) overlaps(idx) = 0.0;
        if(eigvec.get_energy() < energy_llim) overlaps(idx) = 0.0;
        idx++;
    }
    // Now we have a list of overlaps where nonzero elements correspond to eigvecs inside the energy window
    // Get the index to the highest overlapping element
    double max_overlap_val = overlaps.maxCoeff();
    if(max_overlap_val == 0.0) {
        tools::log->debug("No overlapping eigenstates in given energy window {} to {}.", energy_llim, energy_ulim);
        return std::vector<size_t>();
    }

    std::vector<size_t> best_idx;
    std::vector<double> best_val;
    // We collect the best eigvecs from this list until the squared sum of their overlaps goes above a certain threshold.
    // Note that the squared sum overlap = 1 if the states in the list form a complete basis for the current state.
    while(true) {
        max_overlap_val = overlaps.maxCoeff(&idx);
        if(max_overlap_val == 0.0) break;
        best_idx.emplace_back(idx);
        const auto &eigvec = *std::next(eigvecs.begin(), static_cast<long>(idx));
        if(eigvec.is_basis_vector) best_val.emplace_back(max_overlap_val);
        if(best_idx.size() > max_eigvecs) break;
        double sq_sum_overlap = overlaps.cwiseAbs2().sum();
        if(sq_sum_overlap > 0.6) break; // Half means cat state.
        overlaps(idx) = 0.0;            // Zero out the current best so the next-best is found in the next iteration
    }
    tools::log->debug("Found {} eigvecs with good overlap in energy window", best_idx.size());
    return best_idx;
}

Eigen::MatrixXcd tools::finite::opt::internal::subspace::get_eigvecs(const std::vector<opt_mps> &eigvecs) {
    long             rows = eigvecs.front().get_tensor().size();
    long             cols = static_cast<long>(eigvecs.size());
    Eigen::MatrixXcd eigvecs_mat(rows, cols);
    long             idx = 0;
    for(const auto &eigvec : eigvecs) {
        if(not eigvec.is_basis_vector) continue;
        eigvecs_mat.col(idx++) = Eigen::Map<const Eigen::VectorXcd>(eigvec.get_tensor().data(), rows);
    }
    eigvecs_mat.conservativeResize(rows, idx);
    return eigvecs_mat;
}

Eigen::VectorXd tools::finite::opt::internal::subspace::get_eigvals(const std::vector<opt_mps> &eigvecs) {
    long            size = static_cast<long>(eigvecs.size());
    Eigen::VectorXd eigvals(size);
    long            idx = 0;
    for(const auto &eigvec : eigvecs) {
        if(not eigvec.is_basis_vector) continue;
        eigvals(idx++) = eigvec.get_eigval();
    }
    eigvals.conservativeResize(idx);
    return eigvals;
}

Eigen::VectorXd tools::finite::opt::internal::subspace::get_energies(const std::vector<opt_mps> &eigvecs) {
    Eigen::VectorXd energies(static_cast<long>(eigvecs.size()));
    long            idx = 0;
    for(const auto &eigvec : eigvecs) {
        if(not eigvec.is_basis_vector) continue;
        energies(idx++) = eigvec.get_energy();
    }
    energies.conservativeResize(idx);
    return energies;
}

double tools::finite::opt::internal::subspace::get_subspace_error(const std::vector<opt_mps> &eigvecs, std::optional<size_t> max_eigvecs) {
    double eps         = 0;
    size_t num_eigvecs = 0;
    if(not max_eigvecs) max_eigvecs = eigvecs.size();
    for(const auto &eigvec : eigvecs) {
        if(not eigvec.is_basis_vector) continue;
        eps += std::pow(eigvec.get_overlap(), 2);
        num_eigvecs += 1;
        if(num_eigvecs >= max_eigvecs.value()) break;
    }
    return 1.0 - eps;
}

std::vector<double> tools::finite::opt::internal::subspace::get_subspace_errors(const std::vector<opt_mps> &eigvecs) {
    double              eps = 0;
    std::vector<double> subspace_errors;
    subspace_errors.reserve(eigvecs.size());
    for(const auto &eigvec : eigvecs) {
        if(not eigvec.is_basis_vector) continue;
        eps += std::pow(eigvec.get_overlap(), 2);
        subspace_errors.emplace_back(1 - eps);
    }
    return subspace_errors;
}

double tools::finite::opt::internal::subspace::get_subspace_error(const std::vector<double> &overlaps) {
    double eps = 0;
    for(const auto &overlap : overlaps) { eps += std::pow(overlap, 2); }
    return 1.0 - eps;
}

Eigen::VectorXd tools::finite::opt::internal::subspace::get_overlaps(const std::vector<opt_mps> &eigvecs) {
    Eigen::VectorXd overlaps(static_cast<long>(eigvecs.size()));
    long            idx = 0;
    for(const auto &eigvec : eigvecs) {
        if(not eigvec.is_basis_vector) continue;
        overlaps(idx++) = eigvec.get_overlap();
    }
    overlaps.conservativeResize(idx);
    return overlaps;
}

std::pair<double, size_t> find_max_overlap(const std::vector<double> &overlaps) {
    auto max_it  = max_element(overlaps.begin(), overlaps.end());
    auto max_val = *max_it;
    auto max_idx = std::distance(overlaps.begin(), max_it);
    return {max_val, static_cast<size_t>(max_idx)};
}

Eigen::VectorXcd tools::finite::opt::internal::subspace::get_vector_in_subspace(const std::vector<opt_mps> &eigvecs, size_t idx) {
    // In this function we project a vector to the subspace spanned by a small set of eigenvectors
    // Essentially this old computation
    //      Eigen::VectorXcd subspace_vector = (eigvecs.adjoint() * fullspace_vector).normalized();

    auto             fullspace_vector = std::next(eigvecs.begin(), static_cast<long>(idx))->get_vector();
    long             subspace_size    = static_cast<long>(eigvecs.size());
    Eigen::VectorXcd subspace_vector(subspace_size);
    long             subspace_idx = 0;
    for(const auto &eigvec : eigvecs) {
        if(not eigvec.is_basis_vector) continue;
        subspace_vector(subspace_idx++) = eigvec.get_vector().dot(fullspace_vector);
    }
    subspace_vector.conservativeResize(subspace_idx);
    return subspace_vector.normalized();
}

Eigen::VectorXcd tools::finite::opt::internal::subspace::get_vector_in_subspace(const std::vector<opt_mps> &eigvecs, const Eigen::VectorXcd &fullspace_vector) {
    // In this function we project a vector to the subspace spanned by a small set of eigenvectors
    // Essentially this old computation
    //      Eigen::VectorXcd subspace_vector = (eigvecs.adjoint() * fullspace_vector).normalized();
    long             subspace_size = static_cast<long>(eigvecs.size());
    Eigen::VectorXcd subspace_vector(subspace_size);
    long             subspace_idx = 0;
    for(const auto &eigvec : eigvecs) {
        if(not eigvec.is_basis_vector) continue;
        subspace_vector(subspace_idx++) = eigvec.get_vector().dot(fullspace_vector);
    }
    subspace_vector.conservativeResize(subspace_idx);
    return subspace_vector.normalized();
}

Eigen::VectorXcd tools::finite::opt::internal::subspace::get_vector_in_fullspace(const std::vector<opt_mps> &eigvecs, const Eigen::VectorXcd &subspace_vector) {
    // In this function we project a subspace vector back to the full space
    // Essentially this old computation
    //     Eigen::VectorXcd fullspace_vector = (eigvecs * subspace_vector.asDiagonal()).rowwise().sum().normalized();
    //
    // The subspace_vector is a list of weights corresponding to each eigenvector that we use to form a linear combination.
    // Therefore all we need to do to obtain the fullspace vector is to actually add up the linear combination with those weights.

    long             fullspace_size = eigvecs.front().get_vector().size();
    Eigen::VectorXcd fullspace_vector(fullspace_size);
    fullspace_vector.setZero();
    long subspace_idx = 0;
    for(const auto &eigvec : eigvecs) {
        if(not eigvec.is_basis_vector) continue;
        fullspace_vector += eigvec.get_vector() * subspace_vector(subspace_idx++);
    }
    return fullspace_vector.normalized();
}

Eigen::Tensor<cplx, 3> tools::finite::opt::internal::subspace::get_tensor_in_fullspace(const std::vector<opt_mps>        &eigvecs,
                                                                                       const Eigen::VectorXcd            &subspace_vector,
                                                                                       const std::array<Eigen::Index, 3> &dims) {
    return Eigen::TensorMap<Eigen::Tensor<cplx, 3>>(subspace::get_vector_in_fullspace(eigvecs, subspace_vector).data(), dims);
}