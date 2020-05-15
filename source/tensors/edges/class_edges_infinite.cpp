//
// Created by david on 2020-05-15.
//

#include "class_edges_infinite.h"
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
void class_edges_infinite::env_pair<env_type>::assert_validity() const {
    if constexpr(has_validity_v<env_type>) {
        L.assert_validity();
        R.assert_validity();
    }
}

size_t class_edges_infinite::get_length() const {
    if(not math::all_equal(eneL->sites, eneR->sites, varL->sites, varR->sites))
        throw std::runtime_error(
            fmt::format("Site mismatch in edges: eneL {} | eneR {} | varL {} | varR {}", eneL->sites, eneR->sites, varL->sites, varR->sites));
    return eneL->sites + eneR->sites + 2;
}

size_t class_edges_infinite::get_position() const {
    if(not math::all_equal(eneL->get_position(), varL->get_position()))
        throw std::runtime_error(fmt::format("Position mismatch in edges: eneL {} | varL {}", eneL->get_position(), varL->get_position()));
    if(not math::all_equal(eneR->get_position(), varR->get_position()))
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
class_edges_infinite::env_pair<class_env_ene> class_edges_infinite::get_ene() { return {*eneL, *eneR}; }
class_edges_infinite::env_pair<class_env_var> class_edges_infinite::get_var() { return {*varL, *varR}; }

class_edges_infinite::env_pair<const Eigen::Tensor<class_edges_infinite::Scalar,3>> class_edges_infinite::get_ene_blk() const { return {eneL->block, eneR->block}; }
class_edges_infinite::env_pair<const Eigen::Tensor<class_edges_infinite::Scalar,4>> class_edges_infinite::get_var_blk() const { return {varL->block, varR->block}; }
class_edges_infinite::env_pair<Eigen::Tensor<class_edges_infinite::Scalar,3>> class_edges_infinite::get_ene_blk() { return {eneL->block, eneR->block}; }
class_edges_infinite::env_pair<Eigen::Tensor<class_edges_infinite::Scalar,4>> class_edges_infinite::get_var_blk() { return {varL->block, varR->block}; }