//
// Created by david on 2020-05-15.
//

#include "class_edges_infinite.h"
#include "class_env_ene.h"
#include "class_env_var.h"
#include <math/num.h>
#include <tools/common/log.h>

class_edges_infinite::class_edges_infinite()
    : eneL(std::make_unique<class_env_ene>()), eneR(std::make_unique<class_env_ene>()), varL(std::make_unique<class_env_var>()),
      varR(std::make_unique<class_env_var>()) {}

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_edges_infinite::~class_edges_infinite()                            = default;            // default dtor
class_edges_infinite::class_edges_infinite(class_edges_infinite &&other) = default;            // default move ctor
class_edges_infinite &class_edges_infinite::operator=(class_edges_infinite &&other) = default; // default move assign
class_edges_infinite::class_edges_infinite(const class_edges_infinite &other)
    : eneL(std::make_unique<class_env_ene>(*other.eneL)), eneR(std::make_unique<class_env_ene>(*other.eneR)),
      varL(std::make_unique<class_env_var>(*other.varL)), varR(std::make_unique<class_env_var>(*other.varR)) {}

class_edges_infinite &class_edges_infinite::operator=(const class_edges_infinite &other) {
    // check for self-assignment
    if(this != &other) {
        eneL = std::make_unique<class_env_ene>(*other.eneL);
        eneR = std::make_unique<class_env_ene>(*other.eneR);
        varL = std::make_unique<class_env_var>(*other.varL);
        varR = std::make_unique<class_env_var>(*other.varR);
    }
    return *this;
}

void class_edges_infinite::initialize() {
    eneL = std::make_unique<class_env_ene>("L", 0);
    eneR = std::make_unique<class_env_ene>("L", 0);
    varL = std::make_unique<class_env_var>("R", 1);
    varR = std::make_unique<class_env_var>("R", 1);
}

void class_edges_infinite::eject_edges() {
    eneL->clear();
    eneR->clear();
    varL->clear();
    varR->clear();
}

template<typename T, typename = std::void_t<>>
struct has_validity : public std::false_type {};
template<typename T>
struct has_validity<T, std::void_t<decltype(std::declval<T>().assertValidity())>> : public std::true_type {};
template<typename T>
inline constexpr bool has_validity_v = has_validity<T>::value;

template<typename env_type>
void class_edges_infinite::env_pair<env_type>::assert_validity() const {
    if constexpr(has_validity_v<env_type>) {
        L.assert_validity();
        R.assert_validity();
    }
}

size_t class_edges_infinite::get_length() const {
    if(not num::all_equal(eneL->get_sites(), eneR->get_sites(), varL->get_sites(), varR->get_sites()))
        throw std::runtime_error(fmt::format("Site mismatch in edges: eneL {} | eneR {} | varL {} | varR {}", eneL->get_sites(), eneR->get_sites(),
                                             varL->get_sites(), varR->get_sites()));
    return eneL->get_sites() + eneR->get_sites() + 2;
}

size_t class_edges_infinite::get_position() const {
    if(not num::all_equal(eneL->get_position(), varL->get_position()))
        throw std::runtime_error(fmt::format("Position mismatch in edges: eneL {} | varL {}", eneL->get_position(), varL->get_position()));
    if(not num::all_equal(eneR->get_position(), varR->get_position()))
        throw std::runtime_error(fmt::format("Position mismatch in edges: eneR {} | varR {}", eneR->get_position(), varR->get_position()));
    return eneL->get_position();
}

/* clang-format off */
bool class_edges_infinite::is_real() const {
    return eneL->is_real() and
           eneR->is_real() and
           varL->is_real() and
           varR->is_real();
}

bool class_edges_infinite::has_nan() const {
    return eneL->has_nan() or
           eneR->has_nan() or
           varL->has_nan() or
           varR->has_nan();
}

void class_edges_infinite::assert_validity() const {
    eneL->assert_validity();
    eneR->assert_validity();
    varL->assert_validity();
    varR->assert_validity();
}
/* clang-format on */

class_edges_infinite::env_pair<const class_env_ene> class_edges_infinite::get_ene() const { return {*eneL, *eneR}; }
class_edges_infinite::env_pair<const class_env_var> class_edges_infinite::get_var() const { return {*varL, *varR}; }
class_edges_infinite::env_pair<class_env_ene>       class_edges_infinite::get_ene() { return {*eneL, *eneR}; }
class_edges_infinite::env_pair<class_env_var>       class_edges_infinite::get_var() { return {*varL, *varR}; }

class_edges_infinite::env_pair<const Eigen::Tensor<class_edges_infinite::Scalar, 3>> class_edges_infinite::get_ene_blk() const {
    return {eneL->get_block(), eneR->get_block()};
}
class_edges_infinite::env_pair<const Eigen::Tensor<class_edges_infinite::Scalar, 3>> class_edges_infinite::get_var_blk() const {
    return {varL->get_block(), varR->get_block()};
}
class_edges_infinite::env_pair<Eigen::Tensor<class_edges_infinite::Scalar, 3>> class_edges_infinite::get_ene_blk() {
    return {eneL->get_block(), eneR->get_block()};
}
class_edges_infinite::env_pair<Eigen::Tensor<class_edges_infinite::Scalar, 3>> class_edges_infinite::get_var_blk() {
    return {varL->get_block(), varR->get_block()};
}