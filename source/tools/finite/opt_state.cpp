//
// Created by david on 2020-05-31.
//

#include "opt_state.h"
#include <tools/common/fmt.h>

using namespace tools::finite::opt;

opt_state::opt_state(const std::string &name_, const Eigen::Tensor<cplx, 3> &tensor_, const std::vector<size_t> &sites_, double eigval_, double energy_reduced_,
                     std::optional<double> variance_, double overlap_, size_t length_)
    : name(name_), tensor(tensor_), sites(sites_), eigval(eigval_), energy_r(energy_reduced_), energy(eigval_ + energy_reduced_), variance(variance_),
      overlap(overlap_), length(length_) {
    norm    = get_vector().norm();
    iter    = 0;
    counter = 0;
    time    = 0;
}

opt_state::opt_state(const std::string &name_, const Eigen::Tensor<cplx, 3> &tensor_, const std::vector<size_t> &sites_, double energy_, double variance_,
                     double overlap_, size_t length_, size_t iter_, size_t counter_, size_t time_)
    : name(name_), tensor(tensor_), sites(sites_), energy(energy_), variance(variance_), overlap(overlap_), length(length_), iter(iter_), counter(counter_),
      time(time_) {
    norm = get_vector().norm();
}

const std::string &opt_state::get_name() const {
    if(name) return name.value();
    else
        throw std::runtime_error("opt_state: name not set");
}

const Eigen::Tensor<opt_state::cplx, 3> &opt_state::get_tensor() const {
    if(tensor) return tensor.value();
    else
        throw std::runtime_error("opt_state: tensor not set");
}

Eigen::Map<const Eigen::VectorXcd> opt_state::get_vector() const { return Eigen::Map<const Eigen::VectorXcd>(get_tensor().data(), get_tensor().size()); }

Eigen::Map<Eigen::VectorXd> opt_state::get_vector_cplx_as_2xreal() {
    if(not tensor) throw std::runtime_error("opt_state: tensor not set");
    return Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double *>(tensor.value().data()), 2 * tensor.value().size());
}

Eigen::VectorXd opt_state::get_vector_cplx_as_1xreal() const {
    if(not tensor) throw std::runtime_error("opt_state: tensor not set");
    return Eigen::Map<const Eigen::VectorXcd>(tensor.value().data(), tensor.value().size()).real();
}

const std::vector<size_t> &opt_state::get_sites() const {
    if(not sites) throw std::runtime_error("opt_state: sites not set");
    return sites.value();
}

double opt_state::get_energy() const {
    if(energy) return energy.value();
    else
        throw std::runtime_error("opt_state: energy not set");
}

double opt_state::get_energy_reduced() const {
    if(energy_r) return energy_r.value();
    else
        throw std::runtime_error("opt_state: energy_r not set");
}

double opt_state::get_energy_per_site() const { return get_energy() / static_cast<double>(get_length()); }

double opt_state::get_eigval() const {
    if(eigval) return eigval.value();
    else
        throw std::runtime_error("opt_state: eigval not set");
}

double opt_state::get_variance() const {
    if(variance) return variance.value();
    else
        throw std::runtime_error("opt_state: variance not set");
}

double opt_state::get_variance_per_site() const { return get_variance() / static_cast<double>(get_length()); }

size_t opt_state::get_iter() const {
    if(iter) return iter.value();
    else
        return 0.0;
}
size_t opt_state::get_counter() const {
    if(counter) return counter.value();
    else
        return 0.0;
}

double opt_state::get_time() const {
    if(time) return time.value();
    else
        return 0.0;
}

double opt_state::get_overlap() const {
    if(overlap) return overlap.value();
    else
        throw std::runtime_error("opt_state: overlap not set");
}

double opt_state::get_norm() const {
    if(norm) return norm.value();
    else
        throw std::runtime_error("opt_state: norm not set");
}

size_t opt_state::get_length() const {
    if(length) return length.value();
    else
        throw std::runtime_error("opt_state: length not set");
}

OptSpace opt_state::get_optspace() const {
    if(optSpace) return optSpace.value();
    else
        throw std::runtime_error("opt_state: optSpace not set");
}
OptMode opt_state::get_optmode() const {
    if(optMode) return optMode.value();
    else
        throw std::runtime_error("opt_state: optMode not set");
}



void opt_state::clear() { tensor = std::nullopt; }
void opt_state::normalize() {
    if(not tensor) throw std::runtime_error("opt_state: tensor not set");
    Eigen::Map<Eigen::VectorXcd> vector(tensor.value().data(), tensor.value().size());
    vector.normalize();
    norm = vector.norm();
}

void opt_state::set_name(const std::string &name_) { name = name_; }
void opt_state::set_tensor(const Eigen::Tensor<cplx, 3> &tensor_) {
    tensor = tensor_;
    norm   = get_vector().norm();
}
void opt_state::set_tensor(const Eigen::VectorXcd &vector, const Eigen::DSizes<long, 3> &dims) {
    tensor = Eigen::TensorMap<const Eigen::Tensor<const cplx, 3>>(vector.data(), dims);
    norm   = get_vector().norm();
}
void opt_state::set_sites(const std::vector<size_t> &sites_) { sites = sites_; }
void opt_state::set_eigval(double eigval_) {
    eigval = eigval_;
    if(energy and not energy_r) energy_r = energy.value() - eigval.value();
    if(energy_r and not energy) energy = eigval.value() + energy_r.value();
}
void opt_state::set_energy_reduced(double energy_reduced_) {
    energy_r = energy_reduced_;
    if(energy and not eigval) eigval = energy.value() - energy_r.value();
    if(eigval and not energy) energy = eigval.value() + energy_r.value();
}

void opt_state::set_energy(double energy_) {
    energy = energy_;
    if(energy_r and not eigval) eigval = energy.value() - energy_r.value();
    if(eigval and not energy_r) energy_r = energy.value() - eigval.value();
}
void opt_state::set_energy_per_site(double energy_per_site_) { set_energy(energy_per_site_ * static_cast<double>(get_length())); }
void opt_state::set_variance(double variance_) { variance = variance_; }
void opt_state::set_variance_per_site(double variance_per_site_) { set_variance(variance_per_site_ * static_cast<double>(get_length())); }
void opt_state::set_overlap(double overlap_) { overlap = overlap_; }
void opt_state::set_length(size_t length_) { length = length_; }
void opt_state::set_iter(size_t iter_) { iter = iter_; }
void opt_state::set_counter(size_t counter_) { counter = counter_; }
void opt_state::set_time(double time_) { time = time_; }

void opt_state::set_tensor_cplx(const double *data, const Eigen::DSizes<long, 3> &dims) {
    // Here we set a complex tensor from double data.
    // In the real representation, a pair of doubles corresponds to a single complex<double>
    // We expect the size of "dims" to corresponds to the number of complex numbers in the data buffer
    tensor = Eigen::TensorMap<const Eigen::Tensor<const cplx, 3>>(reinterpret_cast<const cplx *>(data), dims);
    norm   = get_vector().norm();
}
void opt_state::set_tensor_real(const double *data, const Eigen::DSizes<long, 3> &dims) {
    // Here we set a complex tensor from double data, but this time
    // only the real part is present in the data buffer, so 1 double corresponds to 1 std::complex<double>.
    // Simplest solution is to cast each element to complex before copy.
    tensor = Eigen::TensorMap<const Eigen::Tensor<const real, 3>>(data, dims).cast<cplx>();
    norm   = get_vector().norm();
}

void opt_state::set_optspace(OptSpace optSpace_){
    optSpace = optSpace_;
}
void opt_state::set_optmode(OptMode optMode_){
    optMode = optMode_;

}


void opt_state::validate_candidate() const {
    std::string error_msg;
    /* clang-format off */
    if(not name)     error_msg.append("\t name    \n");
    if(not tensor)   error_msg.append("\t tensor  \n");
    if(not sites)    error_msg.append("\t sites   \n");
    if(not eigval)   error_msg.append("\t eigval  \n");
    if(not energy_r) error_msg.append("\t energy_r\n");
    if(not energy)   error_msg.append("\t energy  \n");
    if(not overlap)  error_msg.append("\t overlap \n");
    if(not norm)     error_msg.append("\t norm    \n");
    if(not length)   error_msg.append("\t length  \n");
    if(not iter)     error_msg.append("\t iter    \n");
    if(not counter)  error_msg.append("\t counter \n");
    if(not time)     error_msg.append("\t time    \n");
    /* clang-format on */
    if(not error_msg.empty()) { throw std::runtime_error(fmt::format("opt_state error: Missing fields:\n{}", error_msg)); }
}

void opt_state::validate_result() const {
    std::string error_msg;
    /* clang-format off */
    if(not name)     error_msg.append("\t name    \n");
    if(not tensor)   error_msg.append("\t tensor  \n");
    if(not sites)    error_msg.append("\t sites   \n");
    if(not eigval)   error_msg.append("\t eigval  \n");
    if(not energy_r) error_msg.append("\t energy_r\n");
    if(not energy)   error_msg.append("\t energy  \n");
    if(not variance) error_msg.append("\t variance\n");
    if(not overlap)  error_msg.append("\t overlap \n");
    if(not norm)     error_msg.append("\t norm    \n");
    if(not length)   error_msg.append("\t length  \n");
    if(not iter)     error_msg.append("\t iter    \n");
    if(not counter)  error_msg.append("\t counter \n");
    if(not time)     error_msg.append("\t time    \n");
    /* clang-format on */
    if(not error_msg.empty()) { throw std::runtime_error(fmt::format("opt_state error: Missing fields:\n{}", error_msg)); }
}

bool opt_state::operator<(const opt_state &rhs) const {
    if(overlap and rhs.overlap){
        if(this->get_overlap() < rhs.get_overlap()) return true;
        else if(this->get_overlap() > rhs.get_overlap())
            return false;
    }
    // If we reached this point then the overlaps are equal (probably 0)
    // To disambiguate we can use the variance
    if (this->variance and rhs.variance)
        return std::abs(this->get_variance()) < std::abs(rhs.get_variance());


    // If we reached this point then the overlaps are equal (probably 0) and there are no variances defined
    // To disambiguate we can use the eigenvalue, which should
    // be spread around 0 as well (because we use reduced energies)
    // Eigenvalues nearest 0 are near the target energy, and while
    // not strictly relevant, it's probably the best we can do here.
    // Of course this is only valid if the eigenvalues are defined
    if (eigval and rhs.eigval)
        return std::abs(this->get_eigval()) < std::abs(rhs.get_eigval());


    // If we have reched this point the overlaps, variances and eigenvalues are not defined.
    // There is probably a logical bug somewhere
    fmt::print("Checking that this opt_state is valid\n");
    this->validate_candidate();
    fmt::print("Checking that rhs opt_state is valid\n");
    rhs.validate_candidate();
    return true;
}

bool opt_state::operator>(const opt_state &rhs) const { return not(*this < rhs); }
