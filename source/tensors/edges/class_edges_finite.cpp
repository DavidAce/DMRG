//
// Created by david on 2020-05-12.
//

#include "class_edges_finite.h"
#include "class_env_ene.h"
#include "class_env_var.h"
#include <math/nmspc_math.h>
#include <tools/common/log.h>

template<typename T, typename = std::void_t<>>
struct has_validity : public std::false_type {};
template<typename T>
struct has_validity<T, std::void_t<decltype(std::declval<T>().assertValidity())>> : public std::true_type {};
template<typename T>
inline constexpr bool has_validity_v = has_validity<T>::value;

template<typename env_type>
void class_edges_finite::env_pair<env_type>::assert_validity() const {
    if constexpr(has_validity_v<env_type>) {
        L.assert_validity();
        R.assert_validity();
    }
}
// Explicit instantiations
//template struct class_edges_finite::env_pair<class_env_ene>;
//template struct class_edges_finite::env_pair<class_env_var>;
//template struct class_edges_finite::env_pair<Eigen::Tensor<class_edges_finite::Scalar,3>>;
//template struct class_edges_finite::env_pair<Eigen::Tensor<class_edges_finite::Scalar,4>>;
//template struct class_edges_finite::env_pair<const class_env_ene>;
//template struct class_edges_finite::env_pair<const class_env_var>;
//template struct class_edges_finite::env_pair<const Eigen::Tensor<class_edges_finite::Scalar,3>>;
//template struct class_edges_finite::env_pair<const Eigen::Tensor<class_edges_finite::Scalar,4>>;


size_t class_edges_finite::get_length() const {
    if(not math::all_equal(eneL.size(), eneR.size(), varL.size(), varR.size()))
        throw std::runtime_error(
            fmt::format("Size mismatch in environments: eneL {} | eneR {} | varL {} | varR {}", eneL.size(), eneR.size(), varL.size(), varR.size()));
    return eneL.size();
}


/* clang-format off */
bool class_edges_finite::is_real() const {
    for(const auto &env : eneL) if(not env->is_real()) return false;
    for(const auto &env : eneR) if(not env->is_real()) return false;
    for(const auto &env : varL) if(not env->is_real()) return false;
    for(const auto &env : varR) if(not env->is_real()) return false;
    return true;
}

bool class_edges_finite::has_nan() const {
    for(const auto &env : eneL) if(not env->has_nan()) return true;
    for(const auto &env : eneR) if(not env->has_nan()) return true;
    for(const auto &env : varL) if(not env->has_nan()) return true;
    for(const auto &env : varR) if(not env->has_nan()) return true;
    return false;
}

void class_edges_finite::assert_validity() const {
    for(const auto &env : eneL) env->assert_validity();
    for(const auto &env : eneR) env->assert_validity();
    for(const auto &env : varL) env->assert_validity();
    for(const auto &env : varR) env->assert_validity();
}
/* clang-format on */



const class_env_ene &class_edges_finite::get_eneL(size_t pos) const {
    if(pos >= get_length()) throw std::range_error(fmt::format("get_eneL(pos): pos is out of range {} | system size {}", pos, get_length()));
    return **std::next(eneL.begin(), static_cast<long>(pos));
}
const class_env_ene &class_edges_finite::get_eneR(size_t pos) const {
    if(pos >= get_length()) throw std::range_error(fmt::format("get_eneR(pos): pos is out of range {} | system size {}", pos, get_length()));
    return **std::next(eneR.begin(), static_cast<long>(pos));
}
const class_env_var &class_edges_finite::get_varL(size_t pos) const {
    if(pos >= get_length()) throw std::range_error(fmt::format("get_varL(pos): pos is out of range {} | system size {}", pos, get_length()));
    return **std::next(varL.begin(), static_cast<long>(pos));
}
const class_env_var &class_edges_finite::get_varR(size_t pos) const {
    if(pos >= get_length()) throw std::range_error(fmt::format("get_varR(pos): pos is out of range {} | system size {}", pos, get_length()));
    return **std::next(varR.begin(), static_cast<long>(pos));
}

class_env_ene &class_edges_finite::get_eneL(size_t pos) { return const_cast<class_env_ene &>(std::as_const(*this).get_eneL(pos)); }
class_env_ene &class_edges_finite::get_eneR(size_t pos) { return const_cast<class_env_ene &>(std::as_const(*this).get_eneR(pos)); }
class_env_var &class_edges_finite::get_varL(size_t pos) { return const_cast<class_env_var &>(std::as_const(*this).get_varL(pos)); }
class_env_var &class_edges_finite::get_varR(size_t pos) { return const_cast<class_env_var &>(std::as_const(*this).get_varR(pos)); }

class_edges_finite::env_pair<const class_env_ene> class_edges_finite::get_ene(size_t posL, size_t posR) const {
    if(posL > posR) throw std::range_error(fmt::format("get_ene(posL,posR): posL is out of range posL {} > posR {}", posL, posR));
    return {get_eneL(posL), get_eneR(posR)};
}
class_edges_finite::env_pair<const class_env_var> class_edges_finite::get_var(size_t posL, size_t posR) const {
    if(posL > posR) throw std::range_error(fmt::format("get_var(posL,posR): posL is out of range posL {} > posR {}", posL, posR));
    return {get_varL(posL), get_varR(posR)};
}

class_edges_finite::env_pair<class_env_ene> class_edges_finite::get_ene(size_t posL, size_t posR) {
    if(posL > posR) throw std::range_error(fmt::format("get_ene(posL,posR): posL is out of range posL {} > posR {}", posL, posR));
    return {get_eneL(posL), get_eneR(posR)};
}
class_edges_finite::env_pair<class_env_var> class_edges_finite::get_var(size_t posL, size_t posR) {
    if(posL > posR) throw std::range_error(fmt::format("get_var(posL,posR): posL is out of range posL {} > posR {}", posL, posR));
    return {get_varL(posL), get_varR(posR)};
}

class_edges_finite::env_pair<const class_env_ene> class_edges_finite::get_ene(size_t pos) const { return {get_eneL(pos), get_eneR(pos)}; }
class_edges_finite::env_pair<const class_env_var> class_edges_finite::get_var(size_t pos) const { return {get_varL(pos), get_varR(pos)}; }
class_edges_finite::env_pair<class_env_ene>       class_edges_finite::get_ene(size_t pos) { return {get_eneL(pos), get_eneR(pos)}; }
class_edges_finite::env_pair<class_env_var>       class_edges_finite::get_var(size_t pos) { return {get_varL(pos), get_varR(pos)}; }

class_edges_finite::env_pair<const Eigen::Tensor<class_edges_finite::Scalar, 3>> class_edges_finite::get_ene_blk(size_t posL, size_t posR) const {
    return {get_ene(posL).L.block, get_ene(posR).R.block};
}

class_edges_finite::env_pair<const Eigen::Tensor<class_edges_finite::Scalar, 4>> class_edges_finite::get_var_blk(size_t posL, size_t posR) const {
    return {get_var(posL).L.block, get_var(posR).R.block};
}

class_edges_finite::env_pair<Eigen::Tensor<class_edges_finite::Scalar, 3>> class_edges_finite::get_ene_blk(size_t posL, size_t posR) {
    return {get_ene(posL).L.block, get_ene(posR).R.block};
}

class_edges_finite::env_pair<Eigen::Tensor<class_edges_finite::Scalar, 4>> class_edges_finite::get_var_blk(size_t posL, size_t posR) {
    return {get_var(posL).L.block, get_var(posR).R.block};
}


class_edges_finite::env_pair<const class_env_ene> class_edges_finite::get_multisite_ene(std::optional<std::list<size_t>> sites) const {
    if(not sites) sites = active_sites;
    if(sites.value().empty()) throw std::runtime_error("Could not get edges: active site list is empty");
    return get_ene(sites.value().front(), sites.value().back());
}

class_edges_finite::env_pair<const class_env_var> class_edges_finite::get_multisite_var(std::optional<std::list<size_t>> sites) const {
    if(not sites) sites = active_sites;
    if(sites.value().empty()) throw std::runtime_error("Could not get edges: active site list is empty");
    return get_var(sites.value().front(), sites.value().back());
}

class_edges_finite::env_pair<class_env_ene> class_edges_finite::get_multisite_ene(std::optional<std::list<size_t>> sites) {
    if(not sites) sites = active_sites;
    if(sites.value().empty()) throw std::runtime_error("Could not get edges: active site list is empty");
    return get_ene(sites.value().front(), sites.value().back());
}

class_edges_finite::env_pair<class_env_var> class_edges_finite::get_multisite_var(std::optional<std::list<size_t>> sites) {
    if(not sites) sites = active_sites;
    if(sites.value().empty()) throw std::runtime_error("Could not get edges: active site list is empty");
    return get_var(sites.value().front(), sites.value().back());
}



class_edges_finite::env_pair<const Eigen::Tensor<class_edges_finite::Scalar, 3>>
    class_edges_finite::get_multisite_ene_blk(std::optional<std::list<size_t>> sites) const {
    const auto &envs = get_multisite_ene(std::move(sites));
    return {envs.L.block, envs.R.block};
}

class_edges_finite::env_pair<const Eigen::Tensor<class_edges_finite::Scalar, 4>>
    class_edges_finite::get_multisite_var_blk(std::optional<std::list<size_t>> sites) const {
    const auto &envs = get_multisite_var(std::move(sites));
    return {envs.L.block, envs.R.block};
}

class_edges_finite::env_pair<Eigen::Tensor<class_edges_finite::Scalar, 3>>
    class_edges_finite::get_multisite_ene_blk(std::optional<std::list<size_t>> sites) {
    auto envs = get_multisite_ene(std::move(sites));
    return {envs.L.block, envs.R.block};
}

class_edges_finite::env_pair<Eigen::Tensor<class_edges_finite::Scalar, 4>>
    class_edges_finite::get_multisite_var_blk(std::optional<std::list<size_t>> sites) {
    auto envs = get_multisite_var(std::move(sites));
    return {envs.L.block, envs.R.block};
}