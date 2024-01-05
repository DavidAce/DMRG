#include "MpoSite.h"
#include "config/debug.h"
#include "debug/exceptions.h"
#include "math/float.h"
#include "math/hash.h"
#include "math/rnd.h"
#include "math/tenx.h"
#include "qm/qm.h"
#include "qm/spin.h"
#include "tools/common/log.h"
#include <h5pp/h5pp.h>
#include <utility>

MpoSite::MpoSite(ModelType model_type_, size_t position_) : model_type(model_type_), position(position_) {}

Eigen::Tensor<cplx, 4> MpoSite::get_non_compressed_mpo_squared() const {
    tools::log->trace("mpo({}): building mpo²", get_position());
    Eigen::Tensor<cplx, 4> mpo = MPO();
    Eigen::Tensor<cplx, 4> mpo2;

    {
        auto d0 = mpo.dimension(0) * mpo.dimension(0);
        auto d1 = mpo.dimension(1) * mpo.dimension(1);
        auto d2 = mpo.dimension(2);
        auto d3 = mpo.dimension(3);
        mpo2    = mpo.contract(mpo.conjugate(), tenx::idx({3}, {2})).shuffle(tenx::array6{0, 3, 1, 4, 2, 5}).reshape(tenx::array4{d0, d1, d2, d3});
    }

    if(parity_shift_sign != 0 and not parity_shift_axus.empty()) {
        // This redefines H² --> H² + Q(σ), where
        //      * Q(σ) = 0.5 * ( I - prod(σ) ) = Proj(-σ), i.e. the "conjugate" projection operator (sign flipped).
        //      * σ is a pauli matrix (usually σ^z)
        //
        // Example: Let H² be the "energy-shifted" Hamiltonian squared, i.e.
        //                  H² = (H-E)²,
        //          which is the one we use as objective function in xDMRG. This is because we can exploit
        //                  H²|ψ⟩ = Var(H) |ψ⟩,
        //          i.e. the smallest eigenvalue of H² is actually the eigenstate |ψ⟩ with energy E with the smallest variance.
        //          Now, let |ψ+⟩ and |ψ-⟩ be degenerate eigenstates with energy E in each sector of a Z2 symmetry of
        //          the Hamiltonian H, along some axis X,Y or Z. It turns out that our eigenvalue solvers for H² have a
        //          very hard time resolving the degeneracy and will usually converge to some superposition
        //                  |ψ⟩ = a|ψ+⟩ + b|ψ-⟩
        //          with indefinite parity (i.e. the global spin component in the relevant axis will likely be near 0).
        //          By adding Q(σ) we lift this degeneracy, since:
        //
        //              (H² + Q(σ)) |ψ+⟩ = (σ² + 0.5(1-1)) |ψ+⟩ = (Var(H) + 0) |ψ+⟩
        //              (H² + Q(σ)) |ψ-⟩ = (σ² + 0.5(1+1)) |ψ-⟩ = (Var(H) + 1) |ψ-⟩
        //
        //          Note:
        //          1) Var(H) is typically a number close to 0, so 1 adds a very large gap.
        //          2) We could in principle add the projecton on the mpo for H instead by defining
        //                      H  -->  (H + iQ(σ))
        //                      H² --> H^† H  = H² + Q(σ)² = H² + Q(σ)
        //             but this introduces imaginaries in all the MPOs which gives us a performance
        //             penalty by forcing us to use complex versions of all the expensive operations.

        auto d0 = mpo2.dimension(0);
        auto d1 = mpo2.dimension(1);
        auto d2 = mpo2.dimension(2);
        auto d3 = mpo2.dimension(3);
        auto pl = qm::spin::half::get_pauli(parity_shift_axus);
        auto id = qm::spin::half::id;

        Eigen::Tensor<cplx, 4> mpo2_with_parity_shift_op(d0 + 2, d1 + 2, d2, d3);
        mpo2_with_parity_shift_op.setZero();
        mpo2_with_parity_shift_op.slice(tenx::array4{0, 0, 0, 0}, mpo2.dimensions())                  = mpo2;
        mpo2_with_parity_shift_op.slice(tenx::array4{d0, d1, 0, 0}, extent4).reshape(extent2)         = tenx::TensorMap(id);
        mpo2_with_parity_shift_op.slice(tenx::array4{d0 + 1, d1 + 1, 0, 0}, extent4).reshape(extent2) = tenx::TensorMap(pl);
        return mpo2_with_parity_shift_op;
    }
    return mpo2;
}

void MpoSite::build_mpo_squared() {
    mpo_squared  = get_non_compressed_mpo_squared();
    unique_id_sq = std::nullopt;
    if constexpr(settings::debug) {
        if(tenx::hasNaN(mpo_squared.value())) {
            print_parameter_names();
            print_parameter_values();
            throw except::runtime_error("MPO squared at position {} has NAN's", get_position());
        }
    }
}

void MpoSite::set_mpo(const Eigen::Tensor<cplx, 4> &mpo) {
    mpo_internal = mpo;
    unique_id    = std::nullopt;
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

const Eigen::Tensor<cplx, 4> &MpoSite::MPO() const {
    if(all_mpo_parameters_have_been_set) {
        return mpo_internal;
    } else {
        throw std::runtime_error("All MPO parameters haven't been set yet.");
    }
}

const Eigen::Tensor<cplx_t, 4> &MpoSite::MPO_t() const {
    if(all_mpo_parameters_have_been_set) {
        return mpo_internal_t;
    } else {
        throw std::runtime_error("All MPO parameters haven't been set yet.");
    }
}

const Eigen::Tensor<cplx, 4> &MpoSite::MPO2() const {
    if(mpo_squared and all_mpo_parameters_have_been_set)
        return mpo_squared.value();
    else
        throw std::runtime_error("MPO squared has not been set.");
}

Eigen::Tensor<cplx, 4> MpoSite::MPO2_nbody_view(std::optional<std::vector<size_t>> nbody, std::optional<std::vector<size_t>> skip) const {
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
        if(param.second.type() == typeid(long double))
            if(std::isnan(std::any_cast<long double>(param.second))) { return true; }
#if defined(USE_QUADMATH)
        if(param.second.type() == typeid(__float128))
            if(isnanq(std::any_cast<__float128>(param.second))) { return true; }
#endif
    }
    return (tenx::hasNaN(mpo_internal));
}

void MpoSite::assert_validity() const {
    for(auto &param : get_parameters()) {
        if(param.second.type() == typeid(double)) {
            if(std::isnan(std::any_cast<double>(param.second))) {
                print_parameter_names();
                print_parameter_values();
                throw except::runtime_error("Param [{}] = {}", param.first, std::any_cast<double>(param.second));
            }
        }
        if(param.second.type() == typeid(long double)) {
            if(std::isnan(std::any_cast<long double>(param.second))) {
                print_parameter_names();
                print_parameter_values();
                throw except::runtime_error("Param [{}] = {}", param.first, std::any_cast<long double>(param.second));
            }
        }
#if defined(USE_QUADMATH)
        if(param.second.type() == typeid(__float128)) {
            if(isnanq(std::any_cast<__float128>(param.second))) {
                print_parameter_names();
                print_parameter_values();
                throw except::runtime_error("Param [{}] = {}", param.first, f128_t(std::any_cast<__float128>(param.second)));
            }
        }
#endif
    }
    if(tenx::hasNaN(mpo_internal)) throw except::runtime_error("MPO has NAN on position {}", get_position());
    if(not tenx::isReal(mpo_internal)) throw except::runtime_error("MPO has IMAG on position {}", get_position());
    if(mpo_squared) {
        if(tenx::hasNaN(mpo_squared.value())) throw except::runtime_error("MPO2 squared has NAN on position {}", get_position());
        if(not tenx::isReal(mpo_squared.value())) throw except::runtime_error("MPO2 squared has IMAG on position {}", get_position());
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

bool MpoSite::is_shifted() const { return e_shift != 0.0; }

bool MpoSite::is_compressed_mpo_squared() const {
    // When H² = mpo*mpo is compressed, we typically find that the virtual bonds
    // have become smaller than they would otherwise. We can simply check that if they are smaller.
    /*           2
     *           |
     *      0---H²---1
     *           |
     *           3
     */
    if(not has_mpo_squared()) return false;
    const auto &mpo    = MPO();
    const auto &mpo_sq = MPO2();
    auto        d0     = mpo.dimension(0);
    auto        d1     = mpo.dimension(1);
    auto        d2     = parity_shift_sign == 0 ? 0 : 2;
    return mpo_sq.dimension(0) != d0 * d0 + d2 or mpo_sq.dimension(1) != d1 * d1 + d2;
}

double MpoSite::get_energy_shift() const { return e_shift; }

void MpoSite::set_energy_shift(double site_energy) {
    if(e_shift != site_energy) {
        e_shift      = site_energy;
        mpo_internal = MPO_shifted_view();
        mpo_squared  = std::nullopt;
        unique_id    = std::nullopt;
        unique_id_sq = std::nullopt;
    }
}

void MpoSite::set_parity_shift_mpo_squared(int sign, std::string_view axis) {
    if(sign == 0) return;
    if(not qm::spin::half::is_valid_axis(axis)) return;
    if(std::abs(sign) != 1) throw except::logic_error("MpoSite::set_parity_shift_mpo_squared: wrong sign value [{}] | expected -1, 0 or 1", sign);
    auto axus = qm::spin::half::get_axis_unsigned(axis);
    if(sign != parity_shift_sign or axus != parity_shift_axus) {
        tools::log->trace("MpoSite[{}]: Setting MPO2 proj: {:+}{}", get_position(), sign, axus);
        parity_shift_sign = sign;
        parity_shift_axus = axus;
        build_mpo_squared();
    }
}

std::pair<int, std::string_view> MpoSite::get_parity_shift_mpo_squared() const { return {parity_shift_sign, parity_shift_axus}; }

Eigen::Tensor<cplx, 1> MpoSite::get_MPO_edge_left() const {
    if(mpo_internal.size() == 0) throw except::runtime_error("mpo({}): can't build left edge: mpo has not been built yet", get_position());
    auto                   ldim = mpo_internal.dimension(0);
    Eigen::Tensor<cplx, 1> ledge(ldim);
    ledge.setZero();
    ledge(ldim - 1) = 1;
    return ledge;
}

Eigen::Tensor<cplx, 1> MpoSite::get_MPO_edge_right() const {
    if(mpo_internal.size() == 0) throw except::runtime_error("mpo({}): can't build right edge: mpo has not been built yet", get_position());
    auto                   rdim = mpo_internal.dimension(1);
    Eigen::Tensor<cplx, 1> redge(rdim);
    redge.setZero();
    redge(0) = 1;
    return redge;
}

Eigen::Tensor<cplx, 1> MpoSite::get_MPO2_edge_left() const {
    auto edge  = get_MPO_edge_left();
    auto d0    = edge.dimension(0);
    auto edge2 = edge.contract(edge, tenx::idx()).reshape(tenx::array1{d0 * d0});
    if(parity_shift_sign != 0) {
        Eigen::Tensor<cplx, 1> edge2_projz(d0 * d0 + 2);
        edge2_projz.setZero();
        edge2_projz.slice(tenx::array1{0}, edge2.dimensions()) = edge2;
        edge2_projz(d0 * d0 + 0)                               = 1.0;
        edge2_projz(d0 * d0 + 1)                               = 1.0;
        return edge2_projz;
    }
    return edge2;
}

Eigen::Tensor<cplx, 1> MpoSite::get_MPO2_edge_right() const {
    auto edge  = get_MPO_edge_right();
    auto d0    = edge.dimension(0);
    auto edge2 = edge.contract(edge.conjugate(), tenx::idx()).reshape(tenx::array1{d0 * d0});
    if(parity_shift_sign != 0) {
        auto                   q = position == 0 ? -parity_shift_sign : 1.0; // Selects the opposite sector sign (only needed on one MPO)
        Eigen::Tensor<cplx, 1> edge2_projz(d0 * d0 + 2);
        edge2_projz.setZero();
        edge2_projz.slice(tenx::array1{0}, edge2.dimensions()) = edge2;
        edge2_projz(d0 * d0 + 0)                               = 0.5;
        edge2_projz(d0 * d0 + 1)                               = 0.5 * q;
        return edge2_projz;
    }
    return edge2;
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
    file.writeAttribute(get_position(), dataset_name, "position");
    for(auto &params : get_parameters()) {
        if(params.second.type() == typeid(double)) file.writeAttribute(std::any_cast<double>(params.second), dataset_name, params.first);
        if(params.second.type() == typeid(size_t)) file.writeAttribute(std::any_cast<size_t>(params.second), dataset_name, params.first);
        if(params.second.type() == typeid(uint64_t)) file.writeAttribute(std::any_cast<uint64_t>(params.second), dataset_name, params.first);
        if(params.second.type() == typeid(int)) file.writeAttribute(std::any_cast<int>(params.second), dataset_name, params.first);
        if(params.second.type() == typeid(bool)) file.writeAttribute(std::any_cast<bool>(params.second), dataset_name, params.first);
        if(params.second.type() == typeid(std::string)) file.writeAttribute(std::any_cast<std::string>(params.second), dataset_name, params.first);
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
        throw except::runtime_error("Could not load MPO. Dataset [{}] does not exist", mpo_dset);
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