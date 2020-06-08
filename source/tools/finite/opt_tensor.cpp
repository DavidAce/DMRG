//
// Created by david on 2020-05-31.
//

#include "opt_tensor.h"
using namespace tools::finite::opt;

opt_tensor::opt_tensor(const std::string &name_, const Eigen::Tensor<cplx, 3> &tensor_, const std::vector<size_t> &sites_, double eigval_,
                       double energy_reduced_, std::optional<double> variance_, double overlap_,size_t length_)
    : name(name_), tensor(tensor_), sites(sites_), eigval(eigval_), energy_r(energy_reduced_), energy(eigval_ + energy_reduced_), variance(variance_),
      overlap(overlap_),length(length_) {
    norm = get_vector().norm();
}

opt_tensor::opt_tensor(const std::string &name_, const Eigen::Tensor<cplx, 3> &tensor_, const std::vector<size_t> &sites_, double energy_, double variance_,
                       double overlap_, size_t length_, size_t iter_, size_t counter_, size_t time_)
    : name(name_), tensor(tensor_), sites(sites_), energy(energy_), variance(variance_), overlap(overlap_), length(length_), iter(iter_), counter(counter_), time(time_) {
    norm = get_vector().norm();
}

const std::string &opt_tensor::get_name() const {
    if(name) return name.value();
    else
        throw std::runtime_error("opt_tensor: name not set");
}

const Eigen::Tensor<opt_tensor::cplx, 3> &opt_tensor::get_tensor() const {
    if(tensor) return tensor.value();
    else
        throw std::runtime_error("opt_tensor: tensor not set");
}

Eigen::Map<const Eigen::VectorXcd> opt_tensor::get_vector() const { return Eigen::Map<const Eigen::VectorXcd>(get_tensor().data(), get_tensor().size()); }

Eigen::Map<Eigen::VectorXd> opt_tensor::get_vector_cplx_as_2xreal() {
    if(not tensor) throw std::runtime_error("opt_tensor: tensor not set");
    return Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double *>(tensor.value().data()), 2 * tensor.value().size());
}

Eigen::VectorXd opt_tensor::get_vector_cplx_as_1xreal() const {
    if(not tensor) throw std::runtime_error("opt_tensor: tensor not set");
    return Eigen::Map<const Eigen::VectorXcd>(tensor.value().data(), tensor.value().size()).real();
}

const std::vector<size_t> &opt_tensor::get_sites() const {
    if(not sites) throw std::runtime_error("opt_tensor: sites not set");
    return sites.value();
}

double opt_tensor::get_energy() const {
    if(energy) return energy.value();
    else
        throw std::runtime_error("opt_tensor: energy not set");
}

double opt_tensor::get_energy_per_site() const {
    return get_energy() /static_cast<double>(get_length());
}


double opt_tensor::get_eigval() const {
    if(eigval) return eigval.value();
    else
        throw std::runtime_error("opt_tensor: eigval not set");
}

double opt_tensor::get_variance() const {
    if(variance) return variance.value();
    else
        throw std::runtime_error("opt_tensor: variance not set");
}

double opt_tensor::get_variance_per_site() const {
    return get_variance()/static_cast<double>(get_length());
}

size_t opt_tensor::get_iter() const {
    if(iter) return iter.value();
    else
        return 0.0;
}
size_t opt_tensor::get_counter() const {
    if(counter) return counter.value();
    else
        return 0.0;
}

double opt_tensor::get_time() const {
    if(time) return time.value();
    else
        return 0.0;
}

double opt_tensor::get_overlap() const {
    if(overlap) return overlap.value();
    else
        throw std::runtime_error("opt_tensor: overlap not set");
}

double opt_tensor::get_norm() const {
    if(norm) return norm.value();
    else
        throw std::runtime_error("opt_tensor: norm not set");
}

size_t opt_tensor::get_length() const {
    if(length) return length.value();
    else         throw std::runtime_error("opt_tensor: length not set");
}

void opt_tensor::clear() { tensor = std::nullopt; }
void opt_tensor::normalize(){
    if(not tensor) throw std::runtime_error("opt_tensor: tensor not set");
    Eigen::Map<Eigen::VectorXcd> vector (tensor.value().data(), tensor.value().size());
    vector.normalize();
    norm = vector.norm();
}


void opt_tensor::set_name(const std::string &name_) { name = name_; }
void opt_tensor::set_tensor(const Eigen::Tensor<cplx, 3> &tensor_) {
    tensor = tensor_;
    norm   = get_vector().norm();
}
void opt_tensor::set_tensor(const Eigen::VectorXcd &vector, const Eigen::DSizes<long, 3> &dims) {
    tensor = Eigen::TensorMap<const Eigen::Tensor<const cplx, 3>>(vector.data(), dims);
    norm   = get_vector().norm();
}
void opt_tensor::set_sites(const std::vector<size_t> &sites_) { sites = sites_; }
void opt_tensor::set_eigval(double eigval_) { eigval = eigval_; }
void opt_tensor::set_energy_reduced(double energy_reduced_) { energy_r = energy_reduced_; }
void opt_tensor::set_energy(double energy_) { energy = energy_; }
void opt_tensor::set_energy_per_site(double energy_per_site_) { set_energy(energy_per_site_ * static_cast<double>(get_length())); }
void opt_tensor::set_variance(double variance_) { variance = variance_; }
void opt_tensor::set_variance_per_site(double variance_per_site_) { set_variance(variance_per_site_ * static_cast<double>(get_length())); }
void opt_tensor::set_overlap(double overlap_) { overlap = overlap_; }
void opt_tensor::set_length(size_t length_) { length = length_; }
void opt_tensor::set_iter(size_t iter_) { iter = iter_; }
void opt_tensor::set_counter(size_t counter_) { counter = counter_; }
void opt_tensor::set_time(double time_) { time = time_; }

void opt_tensor::set_tensor_cplx(const double *data, const Eigen::DSizes<long, 3> &dims) {
    // Here we set a complex tensor from double data.
    // In the real representation, a pair of doubles corresponds to a single complex<double>
    // We expect the size of "dims" to corresponds to the number of complex numbers in the data buffer
    tensor = Eigen::TensorMap<const Eigen::Tensor<const cplx, 3>>(reinterpret_cast<const cplx *>(data), dims);
    norm   = get_vector().norm();
}
void opt_tensor::set_tensor_real(const double *data, const Eigen::DSizes<long, 3> &dims) {
    // Here we set a complex tensor from double data, but this time
    // only the real part is present in the data buffer, so 1 double corresponds to 1 std::complex<double>.
    // Simplest solution is to cast each element to complex before copy.
    tensor = Eigen::TensorMap<const Eigen::Tensor<const real, 3>>(data, dims).cast<cplx>();
    norm   = get_vector().norm();
}

bool opt_tensor::operator<(const opt_tensor &rhs) const {
    bool overlap_eqvl = this->get_overlap() == rhs.get_overlap();
    if(this->get_overlap() < rhs.get_overlap()) return true;
    else if(this->get_overlap() > rhs.get_overlap())
        return false;

    // Now we know that the overlaps are equal (probably 0)
    // To disambiguate we can use the eigenvalue, which should
    // be spread around 0 as well (because we use reduced energies)
    // Eigenvalues nearest 0 are near the target energy, and while
    // not strictly relevant, it's probably the best we can do here.
    return std::abs(this->get_eigval()) < std::abs(rhs.get_overlap());
}

bool opt_tensor::operator>(const opt_tensor &rhs) const { return not(*this < rhs); }
