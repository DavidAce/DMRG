//
// Created by david on 2020-05-31.
//

#include "candidate_tensor.h"
using namespace tools::finite::opt::internal;

candidate_tensor::candidate_tensor(const Eigen::Tensor<cplx, 3> &tensor_, double eigval_, double energy_reduced_ ,std::optional<double> variance_, double overlap_)
    : tensor(tensor_), eigval(eigval_), energy_r(energy_reduced_), energy(eigval_ + energy_reduced_), variance(variance_), overlap(overlap_) {}

const Eigen::Tensor<candidate_tensor::cplx, 3> & candidate_tensor::get_tensor() const {
    if(tensor) return tensor.value();
    else
        throw std::runtime_error("Candidate tensor not set");
}

double candidate_tensor::get_energy() const {
    if(energy) return energy.value();
    else
        throw std::runtime_error("Candidate energy not set");
}

double candidate_tensor::get_eigval() const {
    if(eigval) return eigval.value();
    else
        throw std::runtime_error("Candidate eigval not set");
}

double candidate_tensor::get_variance() const {
    if(variance) return variance.value();
    else
        throw std::runtime_error("Candidate variance not set");
}

void candidate_tensor::set_variance(double var) {
    variance = var;
}

double candidate_tensor::get_overlap() const {
    if(overlap) return overlap.value();
    else
        throw std::runtime_error("Candidate overlap not set");
}

Eigen::Map<const Eigen::VectorXcd> candidate_tensor::get_vector() const {
    return Eigen::Map<const Eigen::VectorXcd>(get_tensor().data(), get_tensor().size());
}

void candidate_tensor::clear() { tensor = std::nullopt; }

bool candidate_tensor::operator<(const candidate_tensor & rhs) const {
    bool overlap_eqvl = this->get_overlap() == rhs.get_overlap();
    if(this->get_overlap() < rhs.get_overlap()) return true;
    else if(this->get_overlap() > rhs.get_overlap()) return false;

    // Now we know that the overlaps are equal (probably 0)
    // To disambiguate we can use the eigenvalue, which should
    // be spread around 0 as well (because we use reduced energies)
    // Eigenvalues nearest 0 are near the target energy, and while
    // not strictly relevant, it's probably the best we can do here.
    return std::abs(this->get_eigval()) < std::abs(rhs.get_overlap());
}

bool candidate_tensor::operator>(const candidate_tensor & rhs) const {
    return not (*this < rhs);
}
