//
// Created by david on 2018-07-04.
//

#include "class_model_base.h"
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <iomanip>
#include <math/nmspc_random.h>

using namespace qm;
using Scalar = std::complex<double>;

class_model_base::class_model_base(size_t position_) { position = position_; }

const Eigen::Tensor<Scalar, 4> &class_model_base::MPO() const {
    if(all_mpo_parameters_have_been_set) {
        return mpo_internal;
    } else {
        throw std::runtime_error("All MPO parameters haven't been set yet.");
    }
}

bool class_model_base::isReal() const { return Textra::isReal(MPO(), "MPO"); }

bool class_model_base::hasNaN() const {
    for(auto &param : get_parameters()) {
        if(param.second.type() == typeid(double))
            if(std::isnan(std::any_cast<double>(param.second))) {
                return true;
            }
    }
    return (Textra::hasNaN(mpo_internal, "MPO"));
}

void class_model_base::assertValidity() const {
    for(auto &param : get_parameters())
        if(param.second.type() == typeid(double))
            if(std::isnan(std::any_cast<double>(param.second))) {
                print_parameter_names();
                print_parameter_values();
                throw std::runtime_error("Param [" + param.first + "] = " + std::to_string(std::any_cast<double>(param.second)) + "");
            }
    if(Textra::hasNaN(mpo_internal, "MPO")) throw std::runtime_error("MPO has NAN on position " + std::to_string(get_position()));
}

bool class_model_base::isReduced() const { return e_reduced != 0.0; }

double class_model_base::get_reduced_energy() const { return e_reduced; }

void class_model_base::set_reduced_energy(double site_energy) {
    e_reduced    = site_energy;
    mpo_internal = MPO_reduced_view();
}

void class_model_base::set_position(size_t position_) { position = position_; }

std::vector<std::string> class_model_base::get_parameter_names() const {
    std::vector<std::string> parameter_names;
    for(auto &item : get_parameters()) parameter_names.push_back(item.first);
    return parameter_names;
}
std::vector<std::any> class_model_base::get_parameter_values() const {
    std::vector<std::any> parameter_values;
    for(auto &item : get_parameters()) parameter_values.push_back(item.second);
    return parameter_values;
}

size_t class_model_base::get_position() const {
    if(position) {
        return position.value();
    } else {
        throw std::runtime_error("Position of MPO has not been set");
    }
}

void class_model_base::print_parameter_names() const {
    std::cout << std::setprecision(10);
    for(auto &item : get_parameters()) std::cout << std::setw(16) << std::left << item.first;
    std::cout << std::endl;
}

void class_model_base::print_parameter_values() const {
    std::cout << std::setprecision(10);
    for(auto &item : get_parameters()) {
        if(item.second.type() == typeid(int)) std::cout << std::setw(16) << std::left << std::any_cast<int>(item.second);
        if(item.second.type() == typeid(bool)) std::cout << std::setw(16) << std::left << std::boolalpha << std::any_cast<bool>(item.second);
        if(item.second.type() == typeid(double)) std::cout << std::setw(16) << std::left << std::any_cast<double>(item.second);
        if(item.second.type() == typeid(size_t)) std::cout << std::setw(16) << std::left << std::any_cast<size_t>(item.second);
        if(item.second.type() == typeid(std::string)) std::cout << std::setw(16) << std::left << std::any_cast<std::string>(item.second);
    }
    std::cout << std::endl;
}
//
const std::any &class_model_base::find_val(const Parameters &parameters, std::string_view key) const {
    for(auto &param : parameters) {
        if(key == param.first) return param.second;
    }
    throw std::runtime_error("No parameter named [" + std::string(key) + "]");
}

std::any &class_model_base::find_val(Parameters &parameters, std::string_view key) const {
    for(auto &param : parameters) {
        if(key == param.first) return param.second;
    }
    throw std::runtime_error("No parameter named [" + std::string(key) + "]");
}