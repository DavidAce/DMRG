#include "MpoSite.h"
#include "config/debug.h"
#include "debug/exceptions.h"
#include "io/fmt_f128_t.h"
#include "math/cast.h"
#include "math/float.h"
#include "math/hash.h"
#include "math/rnd.h"
#include "math/tenx.h"
#include "qm/qm.h"
#include "qm/spin.h"
#include "tools/common/log.h"
#include <config/settings.h>
#include <general/sfinae.h>
#include <h5pp/h5pp.h>
#include <utility>

MpoSite::MpoSite(ModelType model_type_, size_t position_) : model_type(model_type_), position(position_) {}

void MpoSite::build_mpo() {
    mpo_internal = get_mpo(energy_shift_mpo);
    mpo_internal = apply_edge_left(mpo_internal, get_MPO_edge_left(mpo_internal));
    mpo_internal = apply_edge_right(mpo_internal, get_MPO_edge_right(mpo_internal));
    unique_id    = std::nullopt;
}

void MpoSite::build_mpo_t() {
    mpo_internal_t = get_mpo_t(energy_shift_mpo);
    mpo_internal_t = apply_edge_left(mpo_internal_t, get_MPO_edge_left(mpo_internal_t));
    mpo_internal_t = apply_edge_right(mpo_internal_t, get_MPO_edge_right(mpo_internal_t));
}

void MpoSite::build_mpo_squared() {
    mpo_squared  = get_non_compressed_mpo_squared();
    unique_id_sq = std::nullopt;
}

Eigen::Tensor<cplx, 4> MpoSite::get_non_compressed_mpo_squared() const {
    tools::log->trace("mpo({}): building mpo²", get_position());
    Eigen::Tensor<cplx, 4> mpo = MPO_energy_shifted_view(energy_shift_mpo2);
    Eigen::Tensor<cplx, 4> mpo2;

    {
        auto d0 = mpo.dimension(0) * mpo.dimension(0);
        auto d1 = mpo.dimension(1) * mpo.dimension(1);
        auto d2 = mpo.dimension(2);
        auto d3 = mpo.dimension(3);
        mpo2    = mpo.contract(mpo.conjugate(), tenx::idx({3}, {2})).shuffle(tenx::array6{0, 3, 1, 4, 2, 5}).reshape(tenx::array4{d0, d1, d2, d3});
    }

    if(parity_shift_sign_mpo2 != 0 and not parity_shift_axus_mpo2.empty()) {
        // This redefines H² --> H² + Q(σ), where
        //      * Q(σ) = 0.5 * ( I - prod(σ) ) = Proj(-σ), i.e. the "conjugate" projection operator (sign flipped).
        //      * σ is a pauli matrix (usually σ^z)
        // Observe that Q(σ)|ψ+-⟩ = (1 -+ 1) |ψ+-⟩
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
        //  Note:
        //  1) Var(H) is typically a number close to 0, so 1 adds a very large gap.
        //  2) We could in principle add the projection on the mpo for H instead by defining
        //              H  -->  (H + iQ(σ))
        //              H² --> H^† H  = H² + Q(σ)² = H² + Q(σ)
        //     but this introduces imaginaries in all the MPOs which gives us a performance
        //     penalty by forcing us to use complex versions of all the expensive operations.
        //  3) For ground state DMRG (fDMRG) we can add the projection on H directly, as
        //              H --> (H + Q(σ))
        //     such that
        //              (H + Q(σ)) |ψ+⟩ = (σ² + 0.5(1-1)) |ψ+⟩ = (E + 0) |ψ+⟩
        //              (H + Q(σ)) |ψ-⟩ = (σ² + 0.5(1+1)) |ψ-⟩ = (E + 1) |ψ-⟩

        auto d0 = mpo2.dimension(0);
        auto d1 = mpo2.dimension(1);
        auto d2 = mpo2.dimension(2);
        auto d3 = mpo2.dimension(3);
        auto pl = qm::spin::half::get_pauli(parity_shift_axus_mpo2);
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

Eigen::Tensor<cplx_t, 4> MpoSite::get_mpo_t(real_t energy_shift_per_site, std::optional<std::vector<size_t>> nbody,
                                            std::optional<std::vector<size_t>> skip) const {
    tools::log->warn("MpoSite::get_mpo_t(): Pointless upcast {} -> {}", sfinae::type_name<cplx>(), sfinae::type_name<cplx_t>());
    auto mpo = get_mpo(static_cast<double>(energy_shift_per_site), nbody, skip);
    return mpo.unaryExpr([](auto z) { return std::complex<real_t>(static_cast<real_t>(z.real()), static_cast<real_t>(z.imag())); });
}

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

Eigen::Tensor<cplx, 4> MpoSite::MPO_energy_shifted_view(double energy_shift_per_site) const {
    auto mpo_shifted = get_mpo(energy_shift_per_site);
    mpo_shifted      = apply_edge_left(mpo_shifted, get_MPO_edge_left(mpo_shifted));
    mpo_shifted      = apply_edge_right(mpo_shifted, get_MPO_edge_right(mpo_shifted));
    return mpo_shifted;
}

Eigen::Tensor<cplx, 4> MpoSite::MPO_nbody_view(std::optional<std::vector<size_t>> nbody, std::optional<std::vector<size_t>> skip) const {
    return get_mpo(energy_shift_mpo, nbody, skip);
}

Eigen::Tensor<cplx_t, 4> MpoSite::MPO_nbody_view_t(std::optional<std::vector<size_t>> nbody, std::optional<std::vector<size_t>> skip) const {
    return get_mpo_t(energy_shift_mpo, nbody, skip);
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

bool MpoSite::has_energy_shifted_mpo() const { return energy_shift_mpo != 0.0; }
bool MpoSite::has_energy_shifted_mpo2() const { return energy_shift_mpo2 != 0.0; }

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
    auto        d2     = parity_shift_sign_mpo2 == 0 ? 0 : 2;
    if(get_position() == 0)
        return mpo_sq.dimension(1) < d1 * d1 + d2;
    else if(get_position() + 1 == settings::model::model_size)
        return mpo_sq.dimension(0) < d0 * d0 + d2;
    else
        return mpo_sq.dimension(0) < d0 * d0 + d2 or mpo_sq.dimension(1) < d1 * d1 + d2;
}

double MpoSite::get_energy_shift_mpo() const { return energy_shift_mpo; }
double MpoSite::get_energy_shift_mpo2() const { return energy_shift_mpo2; }

void MpoSite::set_energy_shift(double site_energy) {
    if(energy_shift_mpo != site_energy and settings::precision::use_energy_shifted_mpo) {
        energy_shift_mpo = site_energy;
        build_mpo();
    }
    if(energy_shift_mpo2 != site_energy and settings::precision::use_energy_shifted_mpo_squared) {
        energy_shift_mpo2 = site_energy;
        build_mpo_squared();
    }
}

void MpoSite::set_parity_shift_mpo(int sign, std::string_view axis) {
    if(sign == 0) return;
    if(not qm::spin::half::is_valid_axis(axis)) return;
    if(std::abs(sign) != 1) throw except::logic_error("MpoSite::set_parity_shift_mpo: wrong sign value [{}] | expected -1 or 1", sign);
    auto axus = qm::spin::half::get_axis_unsigned(axis);
    if(sign != parity_shift_sign_mpo or axus != parity_shift_axus_mpo) {
        tools::log->trace("MpoSite[{}]: Setting MPO parity shift: {:+}{}", get_position(), sign, axus);
        parity_shift_sign_mpo = sign;
        parity_shift_axus_mpo = axus;
        build_mpo();
    }
}

std::pair<int, std::string_view> MpoSite::get_parity_shift_mpo() const { return {parity_shift_sign_mpo, parity_shift_axus_mpo}; }

void MpoSite::set_parity_shift_mpo_squared(int sign, std::string_view axis) {
    if(sign == 0) return;
    if(not qm::spin::half::is_valid_axis(axis)) return;
    if(std::abs(sign) != 1) throw except::logic_error("MpoSite::set_parity_shift_mpo_squared: wrong sign value [{}] | expected -1 or 1", sign);
    auto axus = qm::spin::half::get_axis_unsigned(axis);
    if(sign != parity_shift_sign_mpo2 or axus != parity_shift_axus_mpo2) {
        tools::log->trace("MpoSite[{}]: Setting MPO² parity shift: {:+}{}", get_position(), sign, axus);
        parity_shift_sign_mpo2 = sign;
        parity_shift_axus_mpo2 = axus;
        build_mpo_squared();
    }
}

std::pair<int, std::string_view> MpoSite::get_parity_shift_mpo_squared() const { return {parity_shift_sign_mpo2, parity_shift_axus_mpo2}; }

template<typename Scalar>
Eigen::Tensor<Scalar, 1> MpoSite::get_MPO_edge_left(const Eigen::Tensor<Scalar, 4> &mpo) const {
    if(mpo.size() == 0) throw except::runtime_error("mpo({}): can't build the left edge: mpo has not been built yet", get_position());
    auto                     ldim = mpo.dimension(0);
    Eigen::Tensor<Scalar, 1> ledge(ldim);
    if(ldim == 1) {
        // Thin edge (it was probably already applied to the left-most MPO
        ledge.setConstant(cplx(1.0, 0.0));
    } else {
        ledge.setZero();
        if(parity_shift_sign_mpo == 0) {
            /*
             *  MPO = |1 0|
             *        |h 1|
             *  So the left edge picks out the row for h
             */
            ledge(ldim - 1) = 1;
        } else {
            /*
             *  MPO = |1 0 0 0|
             *        |h 1 0 0|
             *        |0 0 1 0|
             *        |0 0 0 σ|
             *  So the left edge picks out the row for h, as well as 1 and σ along the diagonal
             */
            ledge(ldim - 3) = 1; // The bottom left corner
            ledge(ldim - 2) = 1;
            ledge(ldim - 1) = 1;
        }
    }
    return ledge;
}
template Eigen::Tensor<cplx, 1>   MpoSite::get_MPO_edge_left(const Eigen::Tensor<cplx, 4> &mpo) const;
template Eigen::Tensor<cplx_t, 1> MpoSite::get_MPO_edge_left(const Eigen::Tensor<cplx_t, 4> &mpo) const;

template<typename Scalar>
Eigen::Tensor<Scalar, 1> MpoSite::get_MPO_edge_right(const Eigen::Tensor<Scalar, 4> &mpo) const {
    if(mpo.size() == 0) throw except::runtime_error("mpo({}): can't build the right edge: mpo has not been built yet", get_position());
    auto                     rdim = mpo.dimension(1);
    Eigen::Tensor<Scalar, 1> redge(rdim);
    if(rdim == 1) {
        // Thin edge (it was probably already applied to the right-most MPO
        redge.setConstant(Scalar(1.0, 0.0));
    } else {
        redge.setZero();

        if(parity_shift_sign_mpo == 0) {
            /*
             *  MPO = |1 0|
             *        |h 1|
             *  So the right edge picks out the column for h
             */
            redge(0) = 1; // The bottom left corner
        } else {
            /*
             *  MPO = |1 0 0 0|
             *        |h 1 0 0|
             *        |0 0 1 0|
             *        |0 0 0 σ|
             *  So the right edge picks out the column for h, as well as 1 and σ along the diagonal
             *  We also put the overall - sign on prod(σ) if this is the first site.
             */
            auto q          = position == 0 ? -parity_shift_sign_mpo : 1.0; // Selects the opposite sector sign (only needed on one MPO)
            redge(0)        = 1;                                            // The bottom left corner of the original non-parity-shifted mpo
            redge(rdim - 2) = 0.5;
            redge(rdim - 1) = 0.5 * q;
        }
    }
    return redge;
}

template Eigen::Tensor<cplx, 1>   MpoSite::get_MPO_edge_right(const Eigen::Tensor<cplx, 4> &mpo) const;
template Eigen::Tensor<cplx_t, 1> MpoSite::get_MPO_edge_right(const Eigen::Tensor<cplx_t, 4> &mpo) const;

template<typename Scalar>
Eigen::Tensor<Scalar, 4> MpoSite::apply_edge_left(const Eigen::Tensor<Scalar, 4> &mpo, const Eigen::Tensor<Scalar, 1> &edgeL) const {
    if(mpo.dimension(0) == 1 or get_position() != 0) return mpo;
    if(mpo.dimension(0) != edgeL.dimension(0))
        throw except::logic_error("apply_edge_left: dimension mismatch: mpo {} | edgeL {}", mpo.dimensions(), edgeL.dimensions());
    auto  tmp     = mpo;
    auto  dim     = tmp.dimensions();
    auto &threads = tenx::threads::get();
    tmp.resize(tenx::array4{1, dim[1], dim[2], dim[3]});
    tmp.device(*threads->dev) = edgeL.reshape(tenx::array2{1, edgeL.size()}).contract(mpo, tenx::idx({1}, {0}));
    return tmp;
}

template Eigen::Tensor<cplx, 4>   MpoSite::apply_edge_left(const Eigen::Tensor<cplx, 4> &mpo, const Eigen::Tensor<cplx, 1> &edgeL) const;
template Eigen::Tensor<cplx_t, 4> MpoSite::apply_edge_left(const Eigen::Tensor<cplx_t, 4> &mpo, const Eigen::Tensor<cplx_t, 1> &edgeL) const;

template<typename Scalar>
Eigen::Tensor<Scalar, 4> MpoSite::apply_edge_right(const Eigen::Tensor<Scalar, 4> &mpo, const Eigen::Tensor<Scalar, 1> &edgeR) const {
    if(mpo.dimension(1) == 1 or get_position() + 1 != settings::model::model_size) return mpo;
    if(mpo.dimension(1) != edgeR.dimension(0))
        throw except::logic_error("apply_edge_right: dimension mismatch: mpo {} | edgeR {}", mpo.dimensions(), edgeR.dimensions());
    auto  tmp     = mpo;
    auto  dim     = tmp.dimensions();
    auto &threads = tenx::threads::get();
    tmp.resize(tenx::array4{dim[0], 1, dim[2], dim[3]});
    tmp.device(*threads->dev) = mpo.contract(edgeR.reshape(tenx::array2{edgeR.size(), 1}), tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2});
    return tmp;
}
template Eigen::Tensor<cplx, 4>   MpoSite::apply_edge_right(const Eigen::Tensor<cplx, 4> &mpo, const Eigen::Tensor<cplx, 1> &edgeR) const;
template Eigen::Tensor<cplx_t, 4> MpoSite::apply_edge_right(const Eigen::Tensor<cplx_t, 4> &mpo, const Eigen::Tensor<cplx_t, 1> &edgeR) const;

template<typename Scalar>
Eigen::Tensor<Scalar, 1> MpoSite::get_MPO_edge_left() const {
    if constexpr(std::is_same_v<Scalar, cplx>)
        return get_MPO_edge_left(mpo_internal);
    else if constexpr(std::is_same_v<Scalar, cplx_t>)
        return get_MPO_edge_left(mpo_internal_t);
    else {
        static_assert(sfinae::invalid_type_v<Scalar>);
        throw std::logic_error("Invalid type");
    }
}
template Eigen::Tensor<cplx, 1>   MpoSite::get_MPO_edge_left() const;
template Eigen::Tensor<cplx_t, 1> MpoSite::get_MPO_edge_left() const;

template<typename Scalar>
Eigen::Tensor<Scalar, 1> MpoSite::get_MPO_edge_right() const {
    if constexpr(std::is_same_v<Scalar, cplx>)
        return get_MPO_edge_right(mpo_internal);
    else if constexpr(std::is_same_v<Scalar, cplx_t>)
        return get_MPO_edge_right(mpo_internal_t);
    else {
        static_assert(sfinae::invalid_type_v<Scalar>);
        throw std::logic_error("Invalid type");
    }
}
template Eigen::Tensor<cplx, 1>   MpoSite::get_MPO_edge_right() const;
template Eigen::Tensor<cplx_t, 1> MpoSite::get_MPO_edge_right() const;

template<typename Scalar>
Eigen::Tensor<Scalar, 1> MpoSite::get_MPO2_edge_left() const {
    auto ledge  = get_MPO_edge_left<Scalar>();
    auto d0     = ledge.dimension(0);
    auto ledge2 = ledge.contract(ledge, tenx::idx()).reshape(tenx::array1{d0 * d0});
    if(parity_shift_sign_mpo2 != 0) {
        auto ledge2_with_shift = Eigen::Tensor<Scalar, 1>(d0 * d0 + 2);
        ledge2_with_shift.setZero();
        ledge2_with_shift.slice(tenx::array1{0}, ledge2.dimensions()) = ledge2;
        ledge2_with_shift(d0 * d0 + 0)                                = 1.0;
        ledge2_with_shift(d0 * d0 + 1)                                = 1.0;
        return ledge2_with_shift;
    }
    return ledge2;
}
template Eigen::Tensor<cplx, 1>   MpoSite::get_MPO2_edge_left() const;
template Eigen::Tensor<cplx_t, 1> MpoSite::get_MPO2_edge_left() const;

template<typename Scalar>
Eigen::Tensor<Scalar, 1> MpoSite::get_MPO2_edge_right() const {
    auto redge = get_MPO_edge_right<Scalar>();
    auto d0    = redge.dimension(0);
    auto edge2 = redge.contract(redge.conjugate(), tenx::idx()).reshape(tenx::array1{d0 * d0});
    if(parity_shift_sign_mpo2 != 0) {
        auto q                 = position == 0 ? -parity_shift_sign_mpo2 : 1.0; // Selects the opposite sector sign (only needed on one MPO)
        auto redge2_with_shift = Eigen::Tensor<Scalar, 1>(d0 * d0 + 2);
        redge2_with_shift.setZero();
        redge2_with_shift.slice(tenx::array1{0}, edge2.dimensions()) = edge2;
        redge2_with_shift(d0 * d0 + 0)                               = 0.5;
        redge2_with_shift(d0 * d0 + 1)                               = 0.5 * q;
        return redge2_with_shift;
    }
    return edge2;
}
template Eigen::Tensor<cplx, 1>   MpoSite::get_MPO2_edge_right() const;
template Eigen::Tensor<cplx_t, 1> MpoSite::get_MPO2_edge_right() const;

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
    file.writeDataset(MPO_energy_shifted_view(0.0), dataset_name, H5D_layout_t::H5D_CONTIGUOUS);
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
    unique_id = hash::hash_buffer(MPO().data(), safe_cast<size_t>(MPO().size()));
    return unique_id.value();
}

std::size_t MpoSite::get_unique_id_sq() const {
    if(unique_id_sq) return unique_id_sq.value();
    unique_id_sq = hash::hash_buffer(MPO2().data(), safe_cast<size_t>(MPO2().size()));
    return unique_id_sq.value();
}