//
// Created by david on 2020-05-12.
//

#include "class_edges_finite.h"
#include "class_env_ene.h"
#include "class_env_var.h"
#include <math/nmspc_math.h>
#include <tools/common/log.h>

class_edges_finite::class_edges_finite() = default; // Can't initialize lists since we don't know the model size yet

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_edges_finite::~class_edges_finite()                                   = default;            // default dtor
class_edges_finite::class_edges_finite(class_edges_finite &&other) noexcept = default;            // default move ctor
class_edges_finite &class_edges_finite::operator=(class_edges_finite &&other) noexcept = default; // default move assign

class_edges_finite::class_edges_finite(const class_edges_finite &other) : active_sites(other.active_sites) {
    eneL.clear();
    eneR.clear();
    varL.clear();
    varR.clear();
    for(const auto &other_eneL : other.eneL) eneL.emplace_back(std::make_unique<class_env_ene>(*other_eneL));
    for(const auto &other_eneR : other.eneR) eneR.emplace_back(std::make_unique<class_env_ene>(*other_eneR));
    for(const auto &other_varL : other.varL) varL.emplace_back(std::make_unique<class_env_var>(*other_varL));
    for(const auto &other_varR : other.varR) varR.emplace_back(std::make_unique<class_env_var>(*other_varR));
}

class_edges_finite &class_edges_finite::operator=(const class_edges_finite &other) {
    // check for self-assignment
    if(this != &other) {
        active_sites = other.active_sites;
        eneL.clear();
        eneR.clear();
        varL.clear();
        varR.clear();
        for(const auto &other_eneL : other.eneL) eneL.emplace_back(std::make_unique<class_env_ene>(*other_eneL));
        for(const auto &other_eneR : other.eneR) eneR.emplace_back(std::make_unique<class_env_ene>(*other_eneR));
        for(const auto &other_varL : other.varL) varL.emplace_back(std::make_unique<class_env_var>(*other_varL));
        for(const auto &other_varR : other.varR) varR.emplace_back(std::make_unique<class_env_var>(*other_varR));
    }
    return *this;
}

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

void class_edges_finite::initialize(size_t model_size) {
    for(size_t pos = 0; pos < model_size; pos++) {
        eneL.emplace_back(std::make_unique<class_env_ene>("L", pos));
        varL.emplace_back(std::make_unique<class_env_var>("L", pos));
        eneR.emplace_front(std::make_unique<class_env_ene>("R", pos));
        varR.emplace_front(std::make_unique<class_env_var>("R", pos));
    }
}

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

void class_edges_finite::eject_inactive_edges(std::optional<std::list<size_t>> sites) {
    if(not sites) sites = active_sites;
    if(not sites or sites->empty()) {
        // If there are no active sites we may just as well
        // eject everything and let the next rebuild take
        // care of it.
        eject_all_edges();
        return;
    }
    for(auto &env : eneL)
        if(env->get_position() > sites->front()) env->clear();
    for(auto &env : varL)
        if(env->get_position() > sites->front()) env->clear();
    for(auto &env : eneR)
        if(env->get_position() < sites->back()) env->clear();
    for(auto &env : varR)
        if(env->get_position() < sites->back()) env->clear();
}

void class_edges_finite::eject_all_edges() {
    for(auto &env : eneL) env->clear();
    for(auto &env : varL) env->clear();
    for(auto &env : eneR) env->clear();
    for(auto &env : varR) env->clear();
}

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

class_edges_finite::env_pair<Eigen::Tensor<class_edges_finite::Scalar, 3>> class_edges_finite::get_multisite_ene_blk(std::optional<std::list<size_t>> sites) {
    auto envs = get_multisite_ene(std::move(sites));
    return {envs.L.block, envs.R.block};
}

class_edges_finite::env_pair<Eigen::Tensor<class_edges_finite::Scalar, 4>> class_edges_finite::get_multisite_var_blk(std::optional<std::list<size_t>> sites) {
    auto envs = get_multisite_var(std::move(sites));
    return {envs.L.block, envs.R.block};
}