#include "MpoSite.h"
#include <debug/exceptions.h>
#include <h5pp/h5pp.h>
#include <math/hash.h>
#include <math/rnd.h>
#include <math/tenx.h>
#include <qm/qm.h>
#include <utility>

MpoSite::MpoSite(ModelType model_type_, size_t position_) : model_type(model_type_), position(position_) {}

Eigen::Tensor<MpoSite::cplx, 4> MpoSite::get_non_compressed_mpo_squared() const {
    tools::log->debug("mpo({}): building mpo²", get_position());
    const auto &mpo = MPO();
    auto        d0  = mpo.dimension(0) * mpo.dimension(0);
    auto        d1  = mpo.dimension(1) * mpo.dimension(1);
    auto        d2  = mpo.dimension(2);
    auto        d3  = mpo.dimension(3);
    return mpo.contract(mpo, tenx::idx({3}, {2})).shuffle(tenx::array6{0, 3, 1, 4, 2, 5}).reshape(tenx::array4{d0, d1, d2, d3});
}

void MpoSite::build_mpo_squared() {
    mpo_squared  = get_non_compressed_mpo_squared();
    unique_id_sq = std::nullopt;
    if(tenx::hasNaN(mpo_squared.value())) {
        print_parameter_names();
        print_parameter_values();
        throw std::runtime_error(fmt::format("MPO squared at position {} has NAN's", get_position()));
    }
}

void MpoSite::set_mpo_squared(const Eigen::Tensor<cplx, 4> &mpo_sq) {
    mpo_squared  = mpo_sq;
    unique_id_sq = std::nullopt;
}

void MpoSite::clear_mpo_squared() {
    mpo_squared  = std::nullopt;
    unique_id_sq = std::nullopt;
}

bool MpoSite::has_mpo_squared() const { return mpo_squared.has_value(); }

const Eigen::Tensor<MpoSite::cplx, 4> &MpoSite::MPO() const {
    if(all_mpo_parameters_have_been_set) {
        return mpo_internal;
    } else {
        throw std::runtime_error("All MPO parameters haven't been set yet.");
    }
}

const Eigen::Tensor<MpoSite::cplx, 4> &MpoSite::MPO2() const {
    if(mpo_squared and all_mpo_parameters_have_been_set)
        return mpo_squared.value();
    else
        throw std::runtime_error("MPO squared has not been set.");
}

Eigen::Tensor<MpoSite::cplx, 4> &MpoSite::MPO2() {
    if(mpo_squared and all_mpo_parameters_have_been_set)
        return mpo_squared.value();
    else {
        build_mpo_squared();
        return mpo_squared.value();
    }
}

Eigen::Tensor<MpoSite::cplx, 4> MpoSite::MPO2_nbody_view(std::optional<std::vector<size_t>> nbody, std::optional<std::vector<size_t>> skip) const {
    if(not nbody) return MPO2();
    auto mpo1 = MPO_nbody_view(nbody, std::move(skip));
    auto dim0 = mpo1.dimension(0) * mpo1.dimension(0);
    auto dim1 = mpo1.dimension(1) * mpo1.dimension(1);
    auto dim2 = mpo1.dimension(2);
    auto dim3 = mpo1.dimension(3);
    return mpo1.contract(mpo1, tenx::idx({3}, {2})).shuffle(tenx::array6{0, 3, 1, 4, 2, 5}).reshape(tenx::array4{dim0, dim1, dim2, dim3});
}

bool MpoSite::is_real() const { return tenx::isReal(MPO()); }

bool MpoSite::has_nan() const {
    for(auto &param : get_parameters()) {
        if(param.second.type() == typeid(double))
            if(std::isnan(std::any_cast<double>(param.second))) { return true; }
    }
    return (tenx::hasNaN(mpo_internal));
}

void MpoSite::assert_validity() const {
    for(auto &param : get_parameters())
        if(param.second.type() == typeid(double))
            if(std::isnan(std::any_cast<double>(param.second))) {
                print_parameter_names();
                print_parameter_values();
                throw std::runtime_error(fmt::format("Param [{}] = {}", param.first, std::any_cast<double>(param.second)));
            }
    if(tenx::hasNaN(mpo_internal)) throw std::runtime_error(fmt::format("MPO has NAN on position {}", get_position()));
    if(not tenx::isReal(mpo_internal)) throw std::runtime_error(fmt::format("MPO has IMAG on position {}", get_position()));
    if(mpo_squared) {
        if(tenx::hasNaN(mpo_squared.value())) throw std::runtime_error(fmt::format("MPO2 squared has NAN on position {}", get_position()));
        if(not tenx::isReal(mpo_squared.value())) throw std::runtime_error(fmt::format("MPO2 squared has IMAG on position {}", get_position()));
    }
}

void MpoSite::set_position(size_t position_) { position = position_; }

std::vector<std::string> MpoSite::get_parameter_names() const {
    std::vector<std::string> parameter_names;
    for(auto &item : get_parameters()) parameter_names.push_back(item.first);
    return parameter_names;
}
std::vector<std::any> MpoSite::get_parameter_values() const {
    std::vector<std::any> parameter_values;
    for(auto &item : get_parameters()) parameter_values.push_back(item.second);
    return parameter_values;
}

size_t MpoSite::get_position() const {
    if(position) {
        return position.value();
    } else {
        throw std::runtime_error("Position of MPO has not been set");
    }
}

bool MpoSite::is_reduced() const { return e_reduced != 0.0; }

bool MpoSite::is_compressed_mpo_squared() const {
    // When H² = mpo*mpo is compressed, we typically find that the virtual bonds
    // have become smaller than they would otherwise. We can simply check that if they are smaller.
    /*           2
     *           |
     *      0---H²---1
     *           |
     *           3
     */

    const auto &mpo         = MPO();
    const auto &mpo_sq      = MPO2();
    const auto  bond_mpo    = std::min(mpo.dimension(0), mpo.dimension(1));
    const auto  bond_mpo_sq = std::min(mpo_sq.dimension(0), mpo_sq.dimension(1));
    return bond_mpo_sq < bond_mpo * bond_mpo;
}

double MpoSite::get_reduced_energy() const { return e_reduced; }

void MpoSite::set_reduced_energy(double site_energy) {
    if(e_reduced != site_energy) {
        e_reduced    = site_energy;
        mpo_internal = MPO_reduced_view();
        mpo_squared  = std::nullopt;
        unique_id    = std::nullopt;
        unique_id_sq = std::nullopt;
    }
}

Eigen::Tensor<MpoSite::cplx, 1> MpoSite::get_MPO_edge_left() const {
    if(mpo_internal.size() == 0) throw except::runtime_error("mpo({}): can't build left edge: mpo has not been built yet", get_position());
    auto                   ldim = mpo_internal.dimension(0);
    Eigen::Tensor<cplx, 1> ledge(ldim);
    ledge.setZero();
    ledge(ldim - 1) = 1;
    return ledge;
}

Eigen::Tensor<MpoSite::cplx, 1> MpoSite::get_MPO_edge_right() const {
    if(mpo_internal.size() == 0) throw except::runtime_error("mpo({}): can't build right edge: mpo has not been built yet", get_position());
    auto                   rdim = mpo_internal.dimension(1);
    Eigen::Tensor<cplx, 1> redge(rdim);
    redge.setZero();
    redge(0) = 1;
    return redge;
}

Eigen::Tensor<MpoSite::cplx, 1> MpoSite::get_MPO2_edge_left() const {
    auto edge = get_MPO_edge_left();
    auto dim  = edge.dimension(0);
    return edge.contract(edge, tenx::idx()).reshape(tenx::array1{dim * dim});
}

Eigen::Tensor<MpoSite::cplx, 1> MpoSite::get_MPO2_edge_right() const {
    auto edge = get_MPO_edge_right();
    auto dim  = edge.dimension(0);
    return edge.contract(edge, tenx::idx()).reshape(tenx::array1{dim * dim});
}

void MpoSite::print_parameter_names() const {
    for(auto &item : get_parameters()) fmt::print("{:<16}", item.first);
    fmt::print("\n");
}

void MpoSite::print_parameter_values() const {
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

void MpoSite::save_mpo(h5pp::File &file, std::string_view mpo_prefix) const {
    std::string dataset_name = fmt::format("{}/H_{}", mpo_prefix, get_position());
    file.writeDataset(MPO(), dataset_name, H5D_layout_t::H5D_CONTIGUOUS);
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

void MpoSite::load_mpo(const h5pp::File &file, std::string_view mpo_prefix) {
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
        if(tenx::VectorMap(MPO()) != tenx::VectorCast(file.readDataset<Eigen::Tensor<cplx, 4>>(mpo_dset)))
            throw std::runtime_error("Built MPO does not match the MPO on file");
    } else {
        throw std::runtime_error(fmt::format("Could not load MPO. Dataset [{}] does not exist", mpo_dset));
    }
}

std::size_t MpoSite::get_unique_id() const {
    if(unique_id) return unique_id.value();
    unique_id = hash::hash_buffer(MPO().data(), static_cast<size_t>(MPO().size()));
    return unique_id.value();
}

std::size_t MpoSite::get_unique_id_sq() const {
    if(unique_id_sq) return unique_id_sq.value();
    unique_id_sq = hash::hash_buffer(MPO2().data(), static_cast<size_t>(MPO2().size()));
    return unique_id_sq.value();
}