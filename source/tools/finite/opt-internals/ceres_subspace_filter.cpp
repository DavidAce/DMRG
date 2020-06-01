#include <general/nmspc_tensor_extra.h>
// -- (textra first)
#include "ceres_subspace_functor.h"
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <iostream>
#include <math/arpack_extra/matrix_product_stl.h>
#include <math/class_eigsolver.h>
#include <math/nmspc_random.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
using namespace tools::finite::opt;
using namespace tools::finite::opt::internal;



void tools::finite::opt::internal::subspace::filter_candidates(std::vector<candidate_tensor> &candidate_list, double maximum_subspace_error, size_t max_accept) {
    // Sort the candidate list in order of descending overlaps. If the overlaps are the same, compare instead the distance in energy to the
    // current energy
    std::sort(candidate_list.begin(),candidate_list.end(), std::greater<>());
    size_t initial_size = candidate_list.size();
    size_t min_accept   = std::min(8ul, candidate_list.size());
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


    tools::log->trace("Filtered from {} down to {} states", initial_size, candidate_list.size());
    tools::log->trace("Subspace error after filter = {}", std::log10(subspace_error));
    if(candidate_list.size() < min_accept) throw std::runtime_error("Filtered too many candidates");
}


//
//void tools::finite::opt::internal::subspace::filter_candidates(std::vector<candidate_tensor> &candidate_list, double maximum_subspace_error, size_t max_accept,
//                                                               double energy_lbound, double energy_ubound) {
//    /* This function filters candidate states down to an acceptable number "max_accept".
//     * We do not need to assume that overlaps or eigvecs/eigvals are sorted in any way, but
//     * we do assume that the elements in overlap correspond to elements in eigvecs/eigvals.
//     * in the same order.
//     * This is a two-step process
//     * Step 1) Collect the highest overlaps in order
//     *          - Copy the highest overlap from overlps_filtered to overlap_accepted (as well as the index, to keep track of positions)
//     *          - Zero out that element in overlap_filtered.
//     *          - Compute the subspace error of overlaps_accepted_val
//     *          - Keep collecting while the number of candidates is still min_accept. Stop collecting if
//     *              - Collected too many ( more than max accepted
//     *              - The subspace error is low enough
//     * Step 2) Now you have a list of indices in overlap_accepted_idx, pointing to
//     *         the candidates in decreasing order of overlap. The indices also correspond
//     *         to columns in the matrix eigvecs (each eigvec is a candidate!)
//     *         Copy these these candidates and return them.
//     *
//     */
//    size_t initial_size = candidate_list.size();
//    size_t min_accept   = std::min(8ul, candidate_list.size());
//    max_accept          = std::min(max_accept, (candidate_list.size()));
//    if(min_accept == max_accept) return;
//    Eigen::VectorXd     energies          = subspace::get_energies(candidate_list);
//    Eigen::VectorXd     overlaps_filtered = subspace::get_overlaps(candidate_list);
//    std::vector<size_t> overlaps_accepted_idx;
//    std::vector<double> overlaps_accepted_val;
//    double              epsilon        = std::numeric_limits<double>::epsilon();
//    double              subspace_error = subspace::get_subspace_error(candidate_list);
//    maximum_subspace_error = epsilon + std::min(subspace_error, maximum_subspace_error); // Make sure you don't actually increase the allowed subspace error
//
//    // Step 1)
//    while(true) {
//        long idx;
//        double overlap = overlaps_filtered.maxCoeff(&idx);
//        if(energies(idx) > energy_lbound and energies(idx) < energy_ubound) {
//            overlaps_accepted_idx.push_back(static_cast<size_t>(idx));
//            overlaps_accepted_val.push_back(overlap);
//            subspace_error = subspace::get_subspace_error(overlaps_accepted_val);
//        }
//
//        overlaps_filtered(idx) = 0;
//        if(overlaps_accepted_val.size() >= min_accept) {
//            if(subspace_error < maximum_subspace_error) break;
//            if(overlaps_accepted_val.size() >= max_accept) break;
//        }
//    }
//
//    // Step 2) Move items from the original list to the final filtered list
//    std::vector<candidate_tensor> candidate_filtered_list;
//    for(const auto &idx : overlaps_accepted_idx) {
//        auto& candidate = *std::next(candidate_list.begin(), static_cast<long>(idx));
//        candidate_filtered_list.emplace_back(candidate);
//        candidate.clear(); // Zero out this element to save memory. This way we emulate a move.
//    }
//    // Now all the important data is in candidate_filtered_list. Just swap the two lists!
//    candidate_list.swap(candidate_filtered_list);
//    if(candidate_list.size() != overlaps_accepted_idx.size())
//        throw std::runtime_error(fmt::format("Size mismatch: {} != {}", candidate_list.size() , overlaps_accepted_idx.size()));
//    for(const auto & candidate : candidate_list)
//        if(candidate.get_tensor().size() == 0) throw std::runtime_error("No tensor in candidate!");
//
//}

std::tuple<Eigen::MatrixXcd, Eigen::VectorXd, Eigen::VectorXd, double>
    tools::finite::opt::internal::subspace::filter_states(const Eigen::MatrixXcd &eigvecs, const Eigen::VectorXd &eigvals, Eigen::VectorXd &overlaps,
                                                          double maximum_subspace_error, size_t max_accept) {
    /* This function filters candidate states down to an acceptable number "max_accept".
     * We do not need to assume that overlaps or eigvecs/eigvals are sorted in any way, but
     * we do assume that the elements in overlap correspond to elements in eigvecs/eigvals.
     * in the same order.
     * This is a two-step process
     * Step 1) Collect the highest overlaps in order
     *          - Copy the highest overlap from overlps_filtered to overlap_accepted (as well as the index, to keep track of positions)
     *          - Zero out that element in overlap_filtered.
     *          - Compute the subspace error of overlaps_accepted_val
     *          - Keep collecting while the number of candidates is still min_accept. Stop collecting if
     *              - Collected too many ( more than max accepted
     *              - The subspace error is low enough
     * Step 2) Now you have a list of indices in overlap_accepted_idx, pointing to
     *         the candidates in decreasing order of overlap. The indices also correspond
     *         to columns in the matrix eigvecs (each eigvec is a candidate!)
     *         Copy these these candidates and return them.
     *
     */

    size_t min_accept = std::min(8ul, static_cast<size_t>(eigvals.size()));
    max_accept        = std::min(max_accept, static_cast<size_t>(eigvals.size()));
    if(min_accept == max_accept) return std::make_tuple(eigvecs, eigvals, overlaps, 1.0 - overlaps.cwiseAbs2().sum());
    Eigen::VectorXd     overlaps_filtered = overlaps;
    std::vector<int>    overlaps_accepted_idx;
    std::vector<double> overlaps_accepted_val;
    double              epsilon        = std::numeric_limits<double>::epsilon();
    double              subspace_error = 1.0 - overlaps.cwiseAbs2().sum();
    maximum_subspace_error = epsilon + std::min(subspace_error, maximum_subspace_error); // Make sure you don't actually increase the allowed subspace error

    // Step 1)
    while(true) {
        int    idx;
        double overlap = overlaps_filtered.maxCoeff(&idx);
        overlaps_accepted_idx.push_back(idx);
        overlaps_accepted_val.push_back(overlap);
        Eigen::Map<Eigen::VectorXd> overlaps_accepted_map(overlaps_accepted_val.data(), static_cast<long>(overlaps_accepted_val.size()));
        subspace_error = 1.0 - overlaps_accepted_map.cwiseAbs2().sum();
        if(overlaps_accepted_val.size() >= min_accept) {
            if(subspace_error < maximum_subspace_error) break;
            if(overlaps_accepted_val.size() >= max_accept) break;
        }
        overlaps_filtered(idx) = 0;
        if(overlaps_filtered.sum() == 0) break;
    }

    // Step 2)
    Eigen::MatrixXcd eigvecs_filtered(eigvecs.rows(), overlaps_accepted_val.size());
    Eigen::VectorXd  eigvals_filtered(overlaps_accepted_val.size());
    overlaps_filtered.resize(static_cast<long>(overlaps_accepted_val.size()));
    int col_num = 0;
    for(auto &idx : overlaps_accepted_idx) {
        eigvecs_filtered.col(col_num) = eigvecs.col(idx);
        eigvals_filtered(col_num)     = eigvals(idx);
        overlaps_filtered(col_num)    = overlaps(idx);
        col_num++;
    }
    tools::log->trace("Filtered from {} down to {} states", eigvals.size(), eigvals_filtered.size());
    tools::log->trace("Subspace error after filter log10(1-eps) = {}", std::log10(epsilon + subspace_error));
    return std::make_tuple(eigvecs_filtered, eigvals_filtered, overlaps_filtered, subspace_error);
}

std::pair<double, int> tools::finite::opt::internal::subspace::get_best_variance_in_window(const class_tensors_finite &tensors, const Eigen::MatrixXcd &eigvecs,
                                                                                           const Eigen::VectorXd &energies_per_site, double lbound,
                                                                                           double ubound) {
    Eigen::VectorXd variances(eigvecs.cols());
    for(long idx = 0; idx < eigvecs.cols(); idx++) {
        if(energies_per_site(idx) <= ubound and energies_per_site(idx) >= lbound) {
            auto multisite_mps = Textra::MatrixTensorMap(eigvecs.col(idx), tensors.state->active_dimensions());
            variances(idx)     = tools::finite::measure::energy_variance_per_site(multisite_mps, tensors);
        } else {
            variances(idx) = std::numeric_limits<double>::infinity();
        }
    }

    if(variances.minCoeff() == std::numeric_limits<double>::infinity()) {
        tools::log->debug("No eigenstates with good variance in given energy window {} to {}.", lbound, ubound);
        tools::log->debug("Subspace energy range is {} to {}.", energies_per_site.minCoeff(), energies_per_site.maxCoeff());
        return std::make_pair(std::numeric_limits<double>::quiet_NaN(), -1);
    }
    int    min_variance_idx;
    double min_variance_val = variances.minCoeff(&min_variance_idx);
    return std::make_pair(min_variance_val, min_variance_idx);
}

std::optional<size_t> tools::finite::opt::internal::subspace::get_idx_to_candidate_with_highest_overlap(const std::vector<candidate_tensor> &candidate_list,
                                                                                                        double energy_lbound, double energy_ubound) {
    if(candidate_list.empty()) return std::nullopt;
    auto   overlaps = get_overlaps(candidate_list);
    long idx      = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        if(candidate.get_energy() > energy_ubound) overlaps(idx) = 0.0;
        if(candidate.get_energy() < energy_lbound) overlaps(idx) = 0.0;
        idx++;
    }
    // Now we have a list of overlaps where nonzero elements correspond do candidates inside the energy window
    // Get the index to the highest overlapping element
    double max_overlap_val = overlaps.maxCoeff(&idx);
    if(max_overlap_val == 0.0) {
        tools::log->debug("No overlapping eigenstates in given energy window {} to {}.", energy_lbound, energy_ubound);
        return std::nullopt;
    }
    return idx;
}

std::pair<double, int> tools::finite::opt::internal::subspace::get_best_overlap_in_window(const Eigen::VectorXd &overlaps,
                                                                                          const Eigen::VectorXd &energies_per_site, double lbound,
                                                                                          double ubound) {
    assert(overlaps.size() == energies_per_site.size() and "get_best_overlap_in_window: Mismatch in overlaps and energies_per_site sizes");
    Eigen::VectorXd overlaps_in_window = overlaps;
    for(long i = 0; i < overlaps.size(); i++) {
        if(energies_per_site(i) > ubound) overlaps_in_window(i) = 0.0;
        if(energies_per_site(i) < lbound) overlaps_in_window(i) = 0.0;
    }
    if(overlaps_in_window.maxCoeff() == 0.0) {
        tools::log->debug("No overlapping eigenstates in given energy window {} to {}.", lbound, ubound);
        tools::log->debug("Subspace energy range is {} to {}.", energies_per_site.minCoeff(), energies_per_site.maxCoeff());
        return std::make_pair(std::numeric_limits<double>::quiet_NaN(), -1);
    }

    int    max_overlap_idx;
    double max_overlap = overlaps_in_window.maxCoeff(&max_overlap_idx);
    return std::make_pair(max_overlap, max_overlap_idx);
}

std::vector<size_t> tools::finite::opt::internal::subspace::get_idx_to_candidates_with_highest_overlap(const std::vector<candidate_tensor> &candidate_list, size_t max_candidates, double energy_lbound,
                                                               double energy_ubound) {
    if(candidate_list.empty()) return std::vector<size_t>();
    auto   overlaps = get_overlaps(candidate_list);
    long idx      = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        if(candidate.get_energy() > energy_ubound) overlaps(idx) = 0.0;
        if(candidate.get_energy() < energy_lbound) overlaps(idx) = 0.0;
        idx++;
    }
    // Now we have a list of overlaps where nonzero elements correspond do candidates inside the energy window
    // Get the index to the highest overlapping element
    double max_overlap_val = overlaps.maxCoeff();
    if(max_overlap_val == 0.0) {
        tools::log->debug("No overlapping eigenstates in given energy window {} to {}.", energy_lbound, energy_ubound);
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

std::vector<std::pair<double, int>> tools::finite::opt::internal::subspace::get_best_candidates_in_window(const Eigen::VectorXd &overlaps,
                                                                                                          const Eigen::VectorXd &energies_per_site,
                                                                                                          double lbound, double ubound) {
    assert(overlaps.size() == energies_per_site.size() and "get_best_overlap_in_window: Mismatch in overlaps and energies_per_site sizes");
    std::vector<std::pair<double, int>> overlaps_in_window;
    for(long i = 0; i < overlaps.size(); i++) {
        if(energies_per_site(i) < ubound and energies_per_site(i) > lbound) overlaps_in_window.emplace_back(std::make_pair(overlaps(i), i));
    }
    if(overlaps_in_window.empty()) {
        tools::log->debug("No candidate eigenstates in given energy window {} to {}.", lbound, ubound);
        tools::log->debug("Subspace energy range is {} to {}.", energies_per_site.minCoeff(), energies_per_site.maxCoeff());
        return overlaps_in_window;
    }

    std::sort(overlaps_in_window.begin(), overlaps_in_window.end()); // Sort in ascending order of overlap
    std::vector<std::pair<double, int>> candidates;

    // We have a list of states in the energy window. We collect "candidates"
    // from this list until the squared sum of their overlaps goes above a certain threshold.
    // Note that the squared sum overlap = 1 if the states in the list form a complete basis for the current state.
    tools::log->debug("Collecting initial guesses");
    auto lambda_sq_sum = [&](double acc, std::pair<double, int> &p) { return acc + p.first * p.first; };
    while(true) {
        if(overlaps_in_window.empty()) break;
        double sq_sum_overlap = std::accumulate(candidates.begin(), candidates.end(), 0.0, lambda_sq_sum);
        tools::log->debug("Sq_sum_overlap:  {:.16f}", sq_sum_overlap);
        if(sq_sum_overlap > 0.6) break; // Half means cat state.
        else {
            candidates.emplace_back(overlaps_in_window.back());
            overlaps_in_window.pop_back();
        }
    }
    tools::log->debug("Found {} candidates with good overlap in energy window", candidates.size());
    return candidates;
}

Eigen::MatrixXcd tools::finite::opt::internal::subspace::get_eigvecs(const std::vector<candidate_tensor> &candidate_list) {
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

Eigen::VectorXd tools::finite::opt::internal::subspace::get_eigvals(const std::vector<candidate_tensor> &candidate_list) {
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

Eigen::VectorXd tools::finite::opt::internal::subspace::get_energies(const std::vector<candidate_tensor> &candidate_list) {
    Eigen::VectorXd energies(static_cast<long>(candidate_list.size()));
    long            idx = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        energies(idx++) = candidate.get_energy();
    }
    energies.conservativeResize(idx);
    return energies;
}

double tools::finite::opt::internal::subspace::get_subspace_error(const std::vector<candidate_tensor> &candidate_list) {
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


Eigen::VectorXd tools::finite::opt::internal::subspace::get_overlaps(const std::vector<candidate_tensor> &candidate_list) {
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

Eigen::VectorXcd tools::finite::opt::internal::subspace::get_vector_in_subspace(const std::vector<candidate_tensor> &candidate_list, size_t idx) {
    // In this function we project a vector to the subspace spanned by a small set of eigenvectors
    // Essentially this old computation
    //      Eigen::VectorXcd subspace_vector = (eigvecs.adjoint() * fullspace_vector).normalized();

    auto             fullspace_vector = std::next(candidate_list.begin(), static_cast<long>(idx))->get_vector();
    long             subspace_size    = static_cast<long>(candidate_list.size());
    Eigen::VectorXcd subspace_vector(subspace_size);
    long             subspace_idx = 0;
    for(const auto &candidate : candidate_list) {
        if(not candidate.is_basis_vector) continue;
        subspace_vector(subspace_idx++) = candidate.get_vector().adjoint() * fullspace_vector;
    }
    subspace_vector.conservativeResize(subspace_idx);
    return subspace_vector.normalized();
}

Eigen::VectorXcd tools::finite::opt::internal::subspace::get_vector_in_fullspace(const std::vector<candidate_tensor> &candidate_list,
                                                                                 const Eigen::VectorXcd &             subspace_vector) {
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
