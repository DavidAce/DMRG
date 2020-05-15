//
// Created by david on 2018-07-04.
//

#include "class_mpo_base.h"
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <h5pp/h5pp.h>
#include <iomanip>
#include <math/nmspc_random.h>

using namespace qm;
using Scalar = std::complex<double>;

class_mpo_base::class_mpo_base(ModelType model_type_, size_t position_)
    : model_type(model_type_), position(position_){}

const Eigen::Tensor<Scalar, 4> &class_mpo_base::MPO() const {
    if(all_mpo_parameters_have_been_set) {
        return mpo_internal;
    } else {
        throw std::runtime_error("All MPO parameters haven't been set yet.");
    }
}

bool class_mpo_base:: is_real() const { return Textra::isReal(MPO(), "MPO"); }

bool class_mpo_base:: has_nan() const {
    for(auto &param : get_parameters()) {
        if(param.second.type() == typeid(double))
            if(std::isnan(std::any_cast<double>(param.second))) {
                return true;
            }
    }
    return (Textra::hasNaN(mpo_internal, "MPO"));
}

void class_mpo_base::assertValidity() const {
    for(auto &param : get_parameters())
        if(param.second.type() == typeid(double))
            if(std::isnan(std::any_cast<double>(param.second))) {
                print_parameter_names();
                print_parameter_values();
                throw std::runtime_error("Param [" + param.first + "] = " + std::to_string(std::any_cast<double>(param.second)) + "");
            }
    if(Textra::hasNaN(mpo_internal, "MPO")) throw std::runtime_error("MPO has NAN on position " + std::to_string(get_position()));
}


void class_mpo_base::set_position(size_t position_) { position = position_; }

std::vector<std::string> class_mpo_base::get_parameter_names() const {
    std::vector<std::string> parameter_names;
    for(auto &item : get_parameters()) parameter_names.push_back(item.first);
    return parameter_names;
}
std::vector<std::any> class_mpo_base::get_parameter_values() const {
    std::vector<std::any> parameter_values;
    for(auto &item : get_parameters()) parameter_values.push_back(item.second);
    return parameter_values;
}

size_t class_mpo_base::get_position() const {
    if(position) {
        return position.value();
    } else {
        throw std::runtime_error("Position of MPO has not been set");
    }
}

bool class_mpo_base::is_damped() const { return alpha != 0.0 or beta != 0.0; }

bool class_mpo_base::is_reduced() const { return e_reduced != 0.0; }

double class_mpo_base::get_reduced_energy() const { return e_reduced; }

void class_mpo_base::set_reduced_energy(double site_energy) {
    e_reduced    = site_energy;
    mpo_internal = MPO_reduced_view();
}


void class_mpo_base::print_parameter_names() const {
    std::cout << std::setprecision(10);
    for(auto &item : get_parameters()) std::cout << std::setw(16) << std::left << item.first;
    std::cout << std::endl;
}

void class_mpo_base::print_parameter_values() const {
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
//const std::any &class_model_base::find_val(const Parameters &parameters, std::string_view key) const {
//    for(auto &param : parameters) {
//        if(key == param.first) return param.second;
//    }
//    throw std::runtime_error("No parameter named [" + std::string(key) + "]");
//}
//
//std::any &class_model_base::find_val(Parameters &parameters, std::string_view key) const {
//    for(auto &param : parameters) {
//        if(key == param.first) return param.second;
//    }
//    throw std::runtime_error("No parameter named [" + std::string(key) + "]");
//}
//
//


void class_mpo_base::write_mpo(h5pp::File &file, const std::string &model_prefix) const {
    std::string mpo_path     = model_prefix + "/mpo";
    std::string dataset_name = mpo_path + "/H_" + std::to_string(get_position());
    file.writeDataset(MPO(), dataset_name, H5D_layout_t::H5D_COMPACT);
    file.writeAttribute(get_position(), "position", dataset_name);
    for(auto &params : get_parameters()) {
        if(params.second.type() == typeid(double)) file.writeAttribute(std::any_cast<double>(params.second), params.first, dataset_name);
        if(params.second.type() == typeid(size_t)) file.writeAttribute(std::any_cast<size_t>(params.second), params.first, dataset_name);
        if(params.second.type() == typeid(uint64_t)) file.writeAttribute(std::any_cast<uint64_t>(params.second), params.first, dataset_name);
        if(params.second.type() == typeid(int)) file.writeAttribute(std::any_cast<int>(params.second), params.first, dataset_name);
        if(params.second.type() == typeid(bool)) file.writeAttribute(std::any_cast<bool>(params.second), params.first, dataset_name);
        if(params.second.type() == typeid(std::string)) file.writeAttribute(std::any_cast<std::string>(params.second), params.first, dataset_name);
    }
}

void class_mpo_base::read_mpo(const h5pp::File &file, const std::string &model_prefix) {
    std::string mpo_dset = model_prefix + "/mpo"+ "H_" + std::to_string(get_position());
    TableMap    map;
    if(file.linkExists(mpo_dset)) {
        auto param_names = file.getAttributeNames(mpo_dset);
        for(auto &param_name : param_names) {
            auto param_type = file.getAttributeTypeInfo(mpo_dset, param_name);
            if(param_type.cpp_type_index) {
                if(param_type.cpp_type_index.value() == typeid(double)) map[param_name] = file.readAttribute<double>(mpo_dset, param_name);
                if(param_type.cpp_type_index.value() == typeid(size_t)) map[param_name] = file.readAttribute<size_t>(mpo_dset, param_name);
                if(param_type.cpp_type_index.value() == typeid(uint64_t)) map[param_name] = file.readAttribute<uint64_t>(mpo_dset, param_name);
                if(param_type.cpp_type_index.value() == typeid(int)) map[param_name] = file.readAttribute<int>(mpo_dset, param_name);
                if(param_type.cpp_type_index.value() == typeid(bool)) map[param_name] = file.readAttribute<bool>(mpo_dset, param_name);
                if(param_type.cpp_type_index.value() == typeid(std::string)) map[param_name] = file.readAttribute<std::string>(mpo_dset, param_name);
            }
        }
        set_parameters(map);
        build_mpo();
        if(Textra::Tensor_to_Vector(MPO()) != Textra::Tensor_to_Vector(file.readDataset<Eigen::Tensor<Scalar,4>>(mpo_dset)))
            throw std::runtime_error("Built MPO does not match the MPO on file");
    }else{
        throw std::runtime_error(fmt::format("Could not load MPO. Dataset [{}] does not exist", mpo_dset));
    }
}
