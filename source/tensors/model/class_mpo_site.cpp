//
// Created by david on 2018-07-04.
//

#include "class_mpo_site.h"
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_tensor_extra.h>
#include <h5pp/h5pp.h>
#include <math/rnd.h>

using namespace qm;
using Scalar = std::complex<double>;

class_mpo_site::class_mpo_site(ModelType model_type_, size_t position_) : model_type(model_type_), position(position_) {}

Eigen::Tensor<Scalar, 4>class_mpo_site::get_uncompressed_mpo() const {
    const auto &mpo = MPO();
    auto        d0  = mpo.dimension(0) * mpo.dimension(0);
    auto        d1  = mpo.dimension(1) * mpo.dimension(1);
    auto        d2  = mpo.dimension(2);
    auto        d3  = mpo.dimension(3);
    return mpo.contract(mpo, Textra::idx({3}, {2})).shuffle(Textra::array6{0, 3, 1, 4, 2, 5}).reshape(Textra::array4{d0, d1, d2, d3});
}

void class_mpo_site::build_mpo_squared() {
    mpo_squared     = get_uncompressed_mpo();
    if(Textra::hasNaN(mpo_squared.value())) {
        print_parameter_names();
        print_parameter_values();
        throw std::runtime_error(fmt::format("MPO squared at position {} has NAN's", get_position()));
    }
}

void class_mpo_site::set_mpo_squared(const Eigen::Tensor<Scalar, 4> &mpo_sq) { mpo_squared = mpo_sq; }

const Eigen::Tensor<Scalar, 4> &class_mpo_site::MPO() const {
    if(all_mpo_parameters_have_been_set) {
        return mpo_internal;
    } else {
        throw std::runtime_error("All MPO parameters haven't been set yet.");
    }
}

const Eigen::Tensor<Scalar, 4> &class_mpo_site::MPO2() const {
    if(mpo_squared and all_mpo_parameters_have_been_set) return mpo_squared.value();
    else
        throw std::runtime_error("MPO squared has not been set.");
}

Eigen::Tensor<Scalar, 4> &class_mpo_site::MPO2() {
    if(mpo_squared and all_mpo_parameters_have_been_set) return mpo_squared.value();
    else {
        build_mpo_squared();
        return mpo_squared.value();
    }
}

bool class_mpo_site::is_real() const { return Textra::isReal(MPO(), "MPO"); }

bool class_mpo_site::has_nan() const {
    for(auto &param : get_parameters()) {
        if(param.second.type() == typeid(double))
            if(std::isnan(std::any_cast<double>(param.second))) { return true; }
    }
    return (Textra::hasNaN(mpo_internal, "MPO"));
}

void class_mpo_site::assert_validity() const {
    for(auto &param : get_parameters())
        if(param.second.type() == typeid(double))
            if(std::isnan(std::any_cast<double>(param.second))) {
                print_parameter_names();
                print_parameter_values();
                throw std::runtime_error(fmt::format("Param [{}] = {}", param.first, std::any_cast<double>(param.second)));
            }
    if(Textra::hasNaN(mpo_internal, "MPO")) throw std::runtime_error(fmt::format("MPO has NAN on position {}", get_position()));
    if(not Textra::isReal(mpo_internal, "MPO")) throw std::runtime_error(fmt::format("MPO has IMAG on position {}", get_position()));
    if(mpo_squared) {
        if(Textra::hasNaN(mpo_squared.value(), "MPO2")) throw std::runtime_error(fmt::format("MPO2 squared has NAN on position {}", get_position()));
        if(not Textra::isReal(mpo_squared.value(), "MPO2")) throw std::runtime_error(fmt::format("MPO2 squared has IMAG on position {}", get_position()));
    }
}

void class_mpo_site::set_position(size_t position_) { position = position_; }

std::vector<std::string> class_mpo_site::get_parameter_names() const {
    std::vector<std::string> parameter_names;
    for(auto &item : get_parameters()) parameter_names.push_back(item.first);
    return parameter_names;
}
std::vector<std::any> class_mpo_site::get_parameter_values() const {
    std::vector<std::any> parameter_values;
    for(auto &item : get_parameters()) parameter_values.push_back(item.second);
    return parameter_values;
}

size_t class_mpo_site::get_position() const {
    if(position) {
        return position.value();
    } else {
        throw std::runtime_error("Position of MPO has not been set");
    }
}

bool class_mpo_site::is_damped() const { return alpha != 0.0 or beta != 0.0; }

bool class_mpo_site::is_reduced() const { return e_reduced != 0.0; }

double class_mpo_site::get_reduced_energy() const { return e_reduced; }

void class_mpo_site::set_reduced_energy(double site_energy) {
    if(e_reduced != site_energy) {
        e_reduced    = site_energy;
        mpo_internal = MPO_reduced_view();
        mpo_squared  = std::nullopt;
    }
}

void class_mpo_site::print_parameter_names() const {
    std::cout << std::setprecision(10);
    for(auto &item : get_parameters()) fmt::print("{:<16}", item.first);
    fmt::print("\n");
}

void class_mpo_site::print_parameter_values() const {
    for(auto &item : get_parameters()) {
        if(item.second.type() == typeid(int)) fmt::print("{:<16}", std::any_cast<int>(item.second));
        if(item.second.type() == typeid(bool)) fmt::print("{:<16}", std::any_cast<bool>(item.second));
        if(item.second.type() == typeid(size_t)) fmt::print("{:<16}", std::any_cast<size_t>(item.second));
        if(item.second.type() == typeid(std::string)) fmt::print("{:<16}", std::any_cast<std::string>(item.second));
        if(item.second.type() == typeid(double)) fmt::print("{:<16.12f}", std::any_cast<double>(item.second));
    }
    fmt::print("\n");
}

//
// const std::any &class_model_base::find_val(const Parameters &parameters, std::string_view key) const {
//    for(auto &param : parameters) {
//        if(key == param.first) return param.second;
//    }
//    throw std::runtime_error("No parameter named [" + std::string(key) + "]");
//}
//
// std::any &class_model_base::find_val(Parameters &parameters, std::string_view key) const {
//    for(auto &param : parameters) {
//        if(key == param.first) return param.second;
//    }
//    throw std::runtime_error("No parameter named [" + std::string(key) + "]");
//}
//
//

void class_mpo_site::save_mpo(h5pp::File &file, const std::string &mpo_prefix) const {
    std::string dataset_name = fmt::format("{}/H_{}", mpo_prefix, get_position());
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

void class_mpo_site::load_mpo(const h5pp::File &file, const std::string &mpo_prefix) {
    std::string mpo_dset = fmt::format("{}/H_{}", mpo_prefix, get_position());
    TableMap    map;
    if(file.linkExists(mpo_dset)) {
        auto param_names = file.getAttributeNames(mpo_dset);
        for(auto &param_name : param_names) {
            auto param_type = file.getTypeInfoAttribute(mpo_dset, param_name);
            if(param_type.cppTypeIndex) {
                if(param_type.cppTypeIndex.value() == typeid(double)) map[param_name] = file.readAttribute<double>(mpo_dset, param_name);
                if(param_type.cppTypeIndex.value() == typeid(size_t)) map[param_name] = file.readAttribute<size_t>(mpo_dset, param_name);
                if(param_type.cppTypeIndex.value() == typeid(uint64_t)) map[param_name] = file.readAttribute<uint64_t>(mpo_dset, param_name);
                if(param_type.cppTypeIndex.value() == typeid(int)) map[param_name] = file.readAttribute<int>(mpo_dset, param_name);
                if(param_type.cppTypeIndex.value() == typeid(bool)) map[param_name] = file.readAttribute<bool>(mpo_dset, param_name);
                if(param_type.cppTypeIndex.value() == typeid(std::string)) map[param_name] = file.readAttribute<std::string>(mpo_dset, param_name);
            }
        }
        set_parameters(map);
        build_mpo();
        if(Textra::Tensor_to_Vector(MPO()) != Textra::Tensor_to_Vector(file.readDataset<Eigen::Tensor<Scalar, 4>>(mpo_dset)))
            throw std::runtime_error("Built MPO does not match the MPO on file");
    } else {
        throw std::runtime_error(fmt::format("Could not load MPO. Dataset [{}] does not exist", mpo_dset));
    }
}
