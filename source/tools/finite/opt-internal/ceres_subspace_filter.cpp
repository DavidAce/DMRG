#include <general/nmspc_tensor_extra.h>
// -- (textra first)
#include "ceres_subspace_functor.h"
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <iostream>
#include <math/rnd.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt_tensor.h>
using namespace tools::finite::opt;
using namespace tools::finite::opt::internal;



void tools::finite::opt::internal::subspace::filter_candidates(std::vector<opt_tensor> &candidate_list, double maximum_subspace_error, size_t max_accept) {
    // Sort the candidate list in order of descending overlaps. If the overlaps are the same, compare instead the distance in energy to the
    // current energy
    std::sort(candidate_list.begin(),candidate_list.end(), std::greater<>());
    size_t initial_size = candidate_list.size();
    size_t min_accept   = std::min(std::min(max_accept,32ul), candidate_list.size());
    max_accept          = std::min(max_accept, candidate_list.size());
    double subspace_error = subspace::get_subspace_error(candidate_list); // This number is as low as it gets now, and will start growing
    while(true){
        if(candidate_list.size() <= max_accept){
            if(candidate_list.size() <= min_accept) break;
            subspace_error = subspace::get_subspace_error(candidate_list);
            if(subspace_error >= maximum_subspace_error) break;
        }
       candidate_list.pop_back();
    }
    size_t idx = 0;
    for(auto & candidate : candidate_list){
        candidate.set_name(fmt::format("eigenvector {}", idx++));
//        tools::log->trace("Filtered {:<16}: overlap {:.16f} | energy {:>20.16f}", candidate.get_name(), candidate.get_overlap(),
//                          candidate.get_energy_per_site());
    }

    tools::log->trace("Filtered from {} down to {} states", initial_size, candidate_list.size());
    tools::log->trace("Subspace error after filter = {}", std::log10(subspace_error));
    if(candidate_list.size() < min_accept) throw std::runtime_error("Filtered too many candidates");
}




std::optional<size_t> tools::finite::opt::internal::subspace::get_idx_to_candidate_with_highest_overlap(const std::vector<opt_tensor> &candidate_list,
                                                                                                        double energy_llim_per_site, double energy_ulim_per_site) {
    if(candidate_list.empty()) return std::nullopt;
    auto   overlaps = get_overlaps(candidate_list);
    long idx      = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        if(candidate.get_energy_per_site() > energy_ulim_per_site) overlaps(idx) = 0.0;
        if(candidate.get_energy_per_site() < energy_llim_per_site) overlaps(idx) = 0.0;
        idx++;
    }
    // Now we have a list of overlaps where nonzero elements correspond do candidates inside the energy window
    // Get the index to the highest overlapping element
    double max_overlap_val = overlaps.maxCoeff(&idx);
    if(max_overlap_val == 0.0) {
        tools::log->debug("No overlapping eigenstates in given energy window {} to {}.", energy_llim_per_site, energy_ulim_per_site);
        return std::nullopt;
    }
    return idx;
}

std::vector<size_t> tools::finite::opt::internal::subspace::get_idx_to_candidates_with_highest_overlap(const std::vector<opt_tensor> &candidate_list, size_t max_candidates, double energy_llim_per_site,
                                                               double energy_ulim_per_site) {
    if(candidate_list.empty()) return std::vector<size_t>();
    auto   overlaps = get_overlaps(candidate_list);
    long idx      = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        if(candidate.get_energy_per_site() > energy_ulim_per_site) overlaps(idx) = 0.0;
        if(candidate.get_energy_per_site() < energy_llim_per_site) overlaps(idx) = 0.0;
        idx++;
    }
    // Now we have a list of overlaps where nonzero elements correspond do candidates inside the energy window
    // Get the index to the highest overlapping element
    double max_overlap_val = overlaps.maxCoeff();
    if(max_overlap_val == 0.0) {
        tools::log->debug("No overlapping eigenstates in given energy window {} to {}.", energy_llim_per_site, energy_ulim_per_site);
        return std::vector<size_t>();
    }

    std::vector<size_t> best_idx;
    std::vector<double> best_val;
    // We collect the best candidates from this list until the squared sum of their overlaps goes above a certain threshold.
    // Note that the squared sum overlap = 1 if the states in the list form a complete basis for the current state.
    while(true) {
        max_overlap_val = overlaps.maxCoeff(&idx);
        if(max_overlap_val == 0.0) break;
        best_idx.emplace_back(idx);
        const auto &candidate = *std::next(candidate_list.begin(), static_cast<long>(idx));
        if(candidate.is_basis_vector) best_val.emplace_back(max_overlap_val);
        if(best_idx.size() > max_candidates) break;
        double sq_sum_overlap = overlaps.cwiseAbs2().sum();
        if(sq_sum_overlap > 0.6) break;  // Half means cat state.
        overlaps(idx) = 0.0; // Zero out the current best so the next-best is found in the next iteration
    }
    tools::log->debug("Found {} candidates with good overlap in energy window", best_idx.size());
    return best_idx;
}


Eigen::MatrixXcd tools::finite::opt::internal::subspace::get_eigvecs(const std::vector<opt_tensor> &candidate_list) {
    long             rows = candidate_list.front().get_tensor().size();
    long             cols = static_cast<long>(candidate_list.size());
    Eigen::MatrixXcd eigvecs(rows, cols);
    long             idx = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        eigvecs.col(idx++) = Eigen::Map<const Eigen::VectorXcd>(candidate.get_tensor().data(), rows);
    }
    eigvecs.conservativeResize(rows, idx);
    return eigvecs;
}

Eigen::VectorXd tools::finite::opt::internal::subspace::get_eigvals(const std::vector<opt_tensor> &candidate_list) {
    long            size = static_cast<long>(candidate_list.size());
    Eigen::VectorXd eigvals(size);
    long            idx = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        eigvals(idx++) = candidate.get_eigval();
    }
    eigvals.conservativeResize(idx);
    return eigvals;
}

Eigen::VectorXd tools::finite::opt::internal::subspace::get_energies(const std::vector<opt_tensor> &candidate_list) {
    Eigen::VectorXd energies(static_cast<long>(candidate_list.size()));
    long            idx = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        energies(idx++) = candidate.get_energy();
    }
    energies.conservativeResize(idx);
    return energies;
}

Eigen::VectorXd tools::finite::opt::internal::subspace::get_energies_per_site(const std::vector<opt_tensor> &candidate_list) {
    Eigen::VectorXd energies(static_cast<long>(candidate_list.size()));
    long            idx = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        energies(idx++) = candidate.get_energy_per_site();
    }
    energies.conservativeResize(idx);
    return energies;
}


double tools::finite::opt::internal::subspace::get_subspace_error(const std::vector<opt_tensor> &candidate_list) {
    double eps = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        eps += std::pow(candidate.get_overlap(), 2);
    }
    return 1.0 - eps;
}


double tools::finite::opt::internal::subspace::get_subspace_error(const std::vector<double> &overlaps) {
    double eps = 0;
    for(const auto &overlap : overlaps) {
        eps += std::pow(overlap, 2);
    }
    return 1.0 - eps;
}


Eigen::VectorXd tools::finite::opt::internal::subspace::get_overlaps(const std::vector<opt_tensor> &candidate_list) {
    Eigen::VectorXd overlaps(static_cast<long>(candidate_list.size()));
    long            idx = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        overlaps(idx++) = candidate.get_overlap();
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

Eigen::VectorXcd tools::finite::opt::internal::subspace::get_vector_in_subspace(const std::vector<opt_tensor> &candidate_list, size_t idx) {
    // In this function we project a vector to the subspace spanned by a small set of eigenvectors
    // Essentially this old computation
    //      Eigen::VectorXcd subspace_vector = (eigvecs.adjoint() * fullspace_vector).normalized();

    auto             fullspace_vector = std::next(candidate_list.begin(), static_cast<long>(idx))->get_vector();
    long             subspace_size    = static_cast<long>(candidate_list.size());
    Eigen::VectorXcd subspace_vector(subspace_size);
    long             subspace_idx = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        subspace_vector(subspace_idx++) = candidate.get_vector().dot(fullspace_vector);
    }
    subspace_vector.conservativeResize(subspace_idx);
    return subspace_vector.normalized();
}

Eigen::VectorXcd tools::finite::opt::internal::subspace::get_vector_in_subspace(const std::vector<opt_tensor> &candidate_list, const Eigen::VectorXcd & fullspace_vector) {
    // In this function we project a vector to the subspace spanned by a small set of eigenvectors
    // Essentially this old computation
    //      Eigen::VectorXcd subspace_vector = (eigvecs.adjoint() * fullspace_vector).normalized();
    long             subspace_size    = static_cast<long>(candidate_list.size());
    Eigen::VectorXcd subspace_vector(subspace_size);
    long             subspace_idx = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        subspace_vector(subspace_idx++) = candidate.get_vector().dot(fullspace_vector);
    }
    subspace_vector.conservativeResize(subspace_idx);
    return subspace_vector.normalized();
}


Eigen::VectorXcd tools::finite::opt::internal::subspace::get_vector_in_fullspace(const std::vector<opt_tensor> &candidate_list,
                                                                                 const Eigen::VectorXcd & subspace_vector) {
    // In this function we project a subspace vector back to the full space
    // Essentially this old computation
    //     Eigen::VectorXcd fullspace_vector = (eigvecs * subspace_vector.asDiagonal()).rowwise().sum().normalized();
    //
    // The subspace_vector is a list of weights corresponding to each eigenvector that we use to form a linear combination.
    // Therefore all we need to do to obtain the fullspace vector is to actually add up the linear combination with those weights.

    long             fullspace_size = candidate_list.front().get_vector().size();
    Eigen::VectorXcd fullspace_vector(fullspace_size);
    fullspace_vector.setZero();
    long subspace_idx = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        fullspace_vector += candidate.get_vector() * subspace_vector(subspace_idx++);
    }
    return fullspace_vector.normalized();
}
