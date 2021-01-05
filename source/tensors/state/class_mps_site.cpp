//
// Created by david on 2019-07-06.
//

#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>
// textra must appear first
#include "class_mps_site.h"
#include <math/hash.h>
#include <tools/common/fmt.h>
#include <utility>

using Scalar = class_mps_site::Scalar;

class_mps_site::class_mps_site() = default;
class_mps_site::class_mps_site(const Eigen::Tensor<Scalar, 3> &M_, const Eigen::Tensor<Scalar, 1> &L_, size_t pos, double error, std::string label_)
    : M(M_), L(L_), position(pos), truncation_error(error), label(std::move(label_)) {}

class_mps_site::class_mps_site(const Eigen::Tensor<Scalar, 3> &M_, std::optional<Eigen::Tensor<Scalar, 1>> L_, size_t pos, double error, std::string label_)
    : M(M_), L(std::move(L_)), position(pos), truncation_error(error), label(std::move(label_)) {}

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_mps_site::~class_mps_site()                      = default;            // default dtor
class_mps_site::class_mps_site(class_mps_site &&other) = default;            // default move ctor
class_mps_site &class_mps_site::operator=(class_mps_site &&other) = default; // default move assign
class_mps_site::class_mps_site(const class_mps_site &other)       = default;
class_mps_site &class_mps_site::operator=(const class_mps_site &other) = default;

bool class_mps_site::isCenter() const { return LC.has_value(); }

bool class_mps_site::has_L() const { return L.has_value(); }
bool class_mps_site::has_M() const { return M.has_value(); }
bool class_mps_site::has_LC() const { return LC.has_value(); }

Eigen::DSizes<long, 3> class_mps_site::dimensions() const { return Eigen::DSizes<long, 3>{spin_dim(), get_chiL(), get_chiR()}; }

bool class_mps_site::is_real() const { return Textra::isReal(get_M_bare(), "M_bare") and Textra::isReal(get_L(), "L"); }
bool class_mps_site::has_nan() const { return Textra::hasNaN(get_M_bare(), "M_bare") and Textra::hasNaN(get_L(), "L"); }
void class_mps_site::assert_validity() const {
    if(has_nan()) throw std::runtime_error(fmt::format("class_mps_site::assert_validity(): MPS (M or L) at position {} has NaN's", get_position()));
}

void class_mps_site::assert_dimensions() const {
    if(get_label() == "B" and get_chiR() != get_L().dimension(0))
        throw std::runtime_error(fmt::format("class_mps_site: Assert failed: Dimensions for B and L are incompatible at site {}: B {} | L {}", get_position(),
                                             get_M_bare().dimensions(), get_L().dimensions()));
    if(get_label() == "A" and get_chiL() != get_L().dimension(0))
        throw std::runtime_error(fmt::format("class_mps_site: Assert failed: Dimensions for A and L are incompatible at site {}: A {} | L {}", get_position(),
                                             get_M_bare().dimensions(), get_L().dimensions()));
    if(get_label() == "AC") {
        if(get_chiL() != get_L().dimension(0))
            throw std::runtime_error(fmt::format("class_mps_site: Assert failed: Dimensions for AC and L are incompatible at site {}: AC {} | L {}",
                                                 get_position(), get_M_bare().dimensions(), get_L().dimensions()));
        if(get_chiR() != get_LC().dimension(0))
            throw std::runtime_error(fmt::format("class_mps_site: Assert failed: Dimensions for AC and L are incompatible at site {}: AC {} | LC{}",
                                                 get_position(), get_M_bare().dimensions(), get_L().dimensions()));
    }
}

void class_mps_site::assert_identity() const {
    if(get_label() == "B") {
        Eigen::Tensor<Scalar, 2> id = get_M_bare().contract(get_M_bare().conjugate(), Textra::idx({0, 2}, {0, 2}));
        if(not Textra::TensorMatrixMap(id).isIdentity(1e-4)) {
            throw std::runtime_error(fmt::format("class_mps_site: B^dagger B is not identity at pos {}", get_position()));
        }
    } else {
        Eigen::Tensor<Scalar, 2> id = get_M_bare().contract(get_M_bare().conjugate(), Textra::idx({0, 1}, {0, 1}));
        if(not Textra::TensorMatrixMap(id).isIdentity(1e-4)) {
            throw std::runtime_error(fmt::format("class_mps_site: A^dagger A is not identity at pos {}", get_position()));
        }
    }
}

const Eigen::Tensor<Scalar, 3> &class_mps_site::get_M_bare() const {
    if(M) {
        if(M.value().size() == 0) throw std::runtime_error(fmt::format("class_mps_site::get_M_bare(): M has size 0 at position {}", get_position()));
        return M.value();
    } else
        throw std::runtime_error(fmt::format("class_mps_site::get_M_bare(): M has not been set at position {}", get_position()));
}
const Eigen::Tensor<Scalar, 3> &class_mps_site::get_M() const {
    if(isCenter()) {
        if(LC.value().dimension(0) != get_M_bare().dimension(2))
            throw std::runtime_error(fmt::format("class_mps_site::get_M(): M and LC dim mismatch: {} != {} at position {}", get_M_bare().dimension(2),
                                                 LC.value().dimension(0), get_position()));
        if(MC) {
            if(MC.value().size() == 0) throw std::runtime_error(fmt::format("class_mps_site::get_M(): MC has size 0 at position {}", get_position()));
            return MC.value();
        } else {
            MC                                   = Eigen::Tensor<Scalar, 3>(spin_dim(), get_chiL(), get_chiR());
            MC->device(Textra::omp::getDevice()) = get_M_bare().contract(Textra::asDiagonal(get_LC()), Textra::idx({2}, {0}));
            if(MC->size() == 0) throw std::runtime_error(fmt::format("class_mps_site::get_M(): built MC with size 0 at position {}", get_position()));
            return MC.value();
        }
    } else
        return get_M_bare();
}

const Eigen::Tensor<Scalar, 1> &class_mps_site::get_L() const {
    if(L) {
        if(L.value().size() == 0) throw std::runtime_error(fmt::format("class_mps_site::get_L(): L has size 0 at position {}", get_position()));
        return L.value();
    } else
        throw std::runtime_error(fmt::format("class_mps_site::get_L(): L has not been set at position {}", get_position()));
}
const Eigen::Tensor<Scalar, 1> &class_mps_site::get_LC() const {
    if(isCenter()) {
        if(LC.value().dimension(0) != get_M_bare().dimension(2))
            throw std::runtime_error(fmt::format("class_mps_site::get_LC(): M dimensions {} are incompatible with LC dimensions {} at position {}",
                                                 get_M_bare().dimensions(), LC.value().dimensions(), get_position()));
        if(LC.value().size() == 0) throw std::runtime_error(fmt::format("class_mps_site::get_LC(): LC has size 0 at position {}", get_position()));
        return LC.value();
    } else
        throw std::runtime_error(fmt::format("class_mps_site::get_LC(): Site at position {} is not a center", get_position()));
}

Eigen::Tensor<Scalar, 3> &class_mps_site::get_M_bare() { return const_cast<Eigen::Tensor<Scalar, 3> &>(std::as_const(*this).get_M_bare()); }
Eigen::Tensor<Scalar, 3> &class_mps_site::get_M() { return const_cast<Eigen::Tensor<Scalar, 3> &>(std::as_const(*this).get_M()); }
Eigen::Tensor<Scalar, 1> &class_mps_site::get_L() { return const_cast<Eigen::Tensor<Scalar, 1> &>(std::as_const(*this).get_L()); }
Eigen::Tensor<Scalar, 1> &class_mps_site::get_LC() { return const_cast<Eigen::Tensor<Scalar, 1> &>(std::as_const(*this).get_LC()); }
double                    class_mps_site::get_truncation_error() const { return truncation_error; }
double                    class_mps_site::get_truncation_error_LC() const { return truncation_error_LC; }
std::string               class_mps_site::get_label() const {
    if(label.empty()) throw std::runtime_error(fmt::format("No label found at position {}", get_position()));
    return label;
}
std::tuple<long, long, long> class_mps_site::get_dims() const { return {spin_dim(), get_chiL(), get_chiR()}; }
long                         class_mps_site::spin_dim() const { return get_M_bare().dimension(0); }
long                         class_mps_site::get_chiL() const { return get_M_bare().dimension(1); }
long                         class_mps_site::get_chiR() const { return get_M_bare().dimension(2); }

void class_mps_site::set_position(const long position_) {
    position = position_;
    MC.reset();
}

template<typename T>
T class_mps_site::get_position() const {
    if(position) {
        return static_cast<T>(position.value());
    } else {
        throw std::runtime_error("class_mps_site::get_position(): Position hasn't been set on mps site.");
    }
}
template size_t class_mps_site::get_position<size_t>() const;
template long   class_mps_site::get_position<long>() const;

void class_mps_site::set_mps(const Eigen::Tensor<Scalar, 3> &M_, const Eigen::Tensor<Scalar, 1> &L_, double error, const std::string &label_) {
    // M has to be a "bare" matrix, i.e. not an MC which would include LC.
    set_M(M_);
    set_L(L_);
    set_truncation_error(error);
    set_label(label_);
}

void class_mps_site::set_M(const Eigen::Tensor<Scalar, 3> &M_) {
    // M has to be a "bare" matrix, i.e. not an MC which would include LC.
    if(position) {
        M = M_;
        MC.reset();
        unique_id = std::nullopt;
    } else
        throw std::runtime_error("class_mps_site::set_M(const Eigen::Tensor<Scalar, 3> &): Can't set M: Position hasn't been set yet");
}
void class_mps_site::set_L(const Eigen::Tensor<Scalar, 1> &L_, double error) {
    if(position) {
        L                = L_;
        truncation_error = error;
        unique_id        = std::nullopt;
    } else
        throw std::runtime_error("Can't set L: Position hasn't been set yet");
}
void class_mps_site::set_L(const std::pair<Eigen::Tensor<Scalar, 1>, double> &L_and_error) { set_L(L_and_error.first, L_and_error.second); }

void class_mps_site::set_LC(const Eigen::Tensor<Scalar, 1> &LC_, double error) {
    if(position) {
        LC = LC_;
        MC.reset();
        truncation_error_LC = error;
        unique_id           = std::nullopt;
        //        tools::log->trace("Setting LC on site {} | size {}", get_position(), LC->dimensions());
        set_label("AC");
    } else
        throw std::runtime_error("Can't set LC: Position hasn't been set yet");
}

void class_mps_site::set_LC(const std::pair<Eigen::Tensor<Scalar, 1>, double> &LC_and_error) { set_LC(LC_and_error.first, LC_and_error.second); }

void class_mps_site::set_truncation_error(double error) { truncation_error = error; }
void class_mps_site::set_truncation_error_LC(double error) { truncation_error_LC = error; }
void class_mps_site::set_label(const std::string &label_) { label = label_; }
void class_mps_site::unset_LC() {
    LC.reset();
    MC.reset();
    unique_id = std::nullopt;
    //    tools::log->trace("Unset LC on site {}",get_position());
    if(label == "AC") label = "A";
}
void class_mps_site::merge_mps(const class_mps_site &other) {
    // This operation is done when merging mps after an svd split, for instance
    if(get_position() != other.get_position())
        throw std::runtime_error(fmt::format("class_mps_site::merge_mps(const class_mps_site &): position mismatch: this pos {} != other pos {}",
                                             get_position(), other.get_position()));

    if(not other.has_M()) throw std::runtime_error("class_mps_site::merge_mps(const class_mps_site &): Got mps site with undefined M");

    // We have to copy the bare "M", i.e. not MC, which would include LC.
    set_M(other.get_M_bare());

    // We should only copy the L matrix if it exists
    // If it does not exist it may just not have been set.
    // For instance, when splitting a multisite tensor into its contituent mps sites,
    // the edge mps sites have empty L matrices on them because these are not updated.
    // This is intended, but should be checked for.
    // Note that if there is no L on the other mps, there has to exist a previous one
    // on this mps -- we can't leave here without having defined L.
    if(other.has_L()) set_L(other.get_L(), other.get_truncation_error());
    else if(not has_L())
        throw std::runtime_error(fmt::format(
            "class_mps_site::merge_mps(const class_mps_site &): Got mps site with undefined L, and no L is defined on this mps either, at position {}",
            get_position()));

    set_label(other.get_label()); // Copy the A/B label

    // If the other is a center it should be fine to copy its contents
    // without checking its size. However for debugging purposes, since
    // that would be unexpected, just give a warning.
    if(other.isCenter()) {
        set_LC(other.get_LC(), other.get_truncation_error_LC());
        if(get_M_bare().dimension(2) != get_LC().dimension(0))
            throw std::runtime_error(
                fmt::format("Got center of wrong dimension. M dims {} are not compatible with LC dims {}", get_M_bare().dimensions(), get_LC().dimensions()));
    }
}

void class_mps_site::apply_mpo(const Eigen::Tensor<Scalar, 4> &mpo) {
    tools::log->trace("Applying mpo (dims {}) at position {} | isCenter: {}", mpo.dimensions(), get_position(), isCenter());
    long mpoDimL = mpo.dimension(0);
    long mpoDimR = mpo.dimension(1);
    if(mpoDimL != mpoDimR) throw std::logic_error("Can't apply mpo's with different L/R dims: not implemented yet");

    if(isCenter()) {
        Eigen::Tensor<Scalar, 1> LC_temp = get_LC().broadcast(Textra::array1{mpoDimR});
        set_LC(LC_temp);
    }
    Eigen::Tensor<Scalar, 3> M_bare_temp(Textra::array3{spin_dim(), get_chiL() * mpoDimL, get_chiR() * mpoDimR});
    M_bare_temp.device(Textra::omp::getDevice()) = get_M_bare()
                                                       .contract(mpo, Textra::idx({0}, {2}))
                                                       .shuffle(Textra::array5{4, 0, 2, 1, 3})
                                                       .reshape(Textra::array3{spin_dim(), get_chiL() * mpoDimL, get_chiR() * mpoDimR});
    Eigen::Tensor<Scalar, 1> L_temp = get_L().broadcast(Textra::array1{mpoDimL});

    set_L(L_temp);
    set_M(M_bare_temp);
}

void class_mps_site::apply_mpo(const Eigen::Tensor<Scalar, 2> &mpo) {
    Eigen::Tensor<Scalar, 3> M_bare_temp(mpo.dimension(1), get_chiL(), get_chiR());
    M_bare_temp.device(Textra::omp::getDevice()) = mpo.contract(get_M_bare(), Textra::idx({0}, {0}));
    set_M(M_bare_temp);
}

void class_mps_site::stash_U(const Eigen::Tensor<Scalar, 3> &U) const { U_stash = U; }
void class_mps_site::stash_S(const Eigen::Tensor<Scalar, 1> &S, double error) const {
    S_stash                  = S;
    truncation_error_S_stash = error;
}
void class_mps_site::stash_S(const std::pair<Eigen::Tensor<Scalar, 1>, double> &S_and_error) const {
    std::tie(S_stash, truncation_error_S_stash) = S_and_error;
}

void class_mps_site::stash_V(const Eigen::Tensor<Scalar, 3> &V) const { V_stash = V; }

bool class_mps_site::has_stash_U() const { return U_stash.has_value(); }
bool class_mps_site::has_stash_S() const { return S_stash.has_value() and truncation_error_S_stash.has_value(); }
bool class_mps_site::has_stash_V() const { return V_stash.has_value(); }

Eigen::Tensor<Scalar, 3> class_mps_site::unstash_U() const {
    if(not U_stash) throw std::runtime_error(fmt::format("No U was found stashed at position {}", get_position()));
    Eigen::Tensor<Scalar, 3> U_tmp = U_stash.value();
    U_stash                        = std::nullopt;
    return U_tmp;
}
Eigen::Tensor<Scalar, 3> class_mps_site::unstash_V() const {
    if(not V_stash) throw std::runtime_error(fmt::format("No V was found stashed at position {}", get_position()));
    Eigen::Tensor<Scalar, 3> V_tmp = V_stash.value();
    V_stash                        = std::nullopt;
    return V_tmp;
}

std::pair<Eigen::Tensor<Scalar, 1>, double> class_mps_site::unstash_S() const {
    if(not S_stash) throw std::runtime_error(fmt::format("No LC was found stashed at position {}", get_position()));
    if(not truncation_error_S_stash) throw std::runtime_error(fmt::format("No LC truncation error was found stashed at position {}", get_position()));
    auto tmp                 = std::make_pair(S_stash.value(), truncation_error_S_stash.value());
    S_stash                  = std::nullopt;
    truncation_error_S_stash = std::nullopt;
    return tmp;
}

void class_mps_site::unstash() const {
    U_stash                  = std::nullopt;
    S_stash                  = std::nullopt;
    V_stash                  = std::nullopt;
    truncation_error_S_stash = std::nullopt;
}

void class_mps_site::merge_stash(const class_mps_site &other) {
    if(get_position() == other.get_position() + 1) {
        /* Left-to-right move.
         * In this case there should be a "V" in "other" that should be absorbed into this site from the left.
         *
         *   1 --[New "M"]-- 2       1 --[V]-- 2   1 --[this "M"]-- 2
         *          |            =        |                |
         *          0                 0(dim=1)             0
         *
         *  In addition, if this is an A site we expect an "S" to come along, which didn't become an LC
         *  for the site on the left. Presumably the true LC is on some site further to the right.
         *  Here we simply set it as the new L of this site.
         */

        auto V = other.unstash_V();
        if(V.dimension(0) != 1)
            throw std::logic_error(fmt::format("Failed to merge stash from left site {} into {}: "
                                               "V has invalid dimensions {}. Dim at idx 0 should be == 1",
                                               other.get_position(), get_position(), V.dimensions()));
        if(V.dimension(2) != get_chiL())
            throw std::logic_error(fmt::format("Failed to merge stash from left site {} into {}: "
                                               "V dimensions {} | M dimensions {} |  Expected V(2) == M(1)",
                                               other.get_position(), get_position(), dimensions()));
        Textra::array3           dims = {spin_dim(), V.dimension(1), get_chiR()};
        Eigen::Tensor<Scalar, 3> VM(dims);
        VM.device(Textra::omp::getDevice()) = V.contract(get_M_bare(), Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(dims);
        set_M(VM);
        if(other.has_stash_S()) {
            if(label == "B") tools::log->warn("Setting L from from the left at pos {} with label {}", get_position(), get_label());
            set_L(other.unstash_S());
        }
    } else if(get_position() == other.get_position() - 1) {
        /* Right-to-left move.
         * In this case there should be a U and an S in "other".
         * U should be absorbed into this site from the right.
         * S becomes the new LC.
         *
         *    1 --[New "M"]-- 2       1 --[this "M"]-- 2   1 --[U]-- 2
         *           |            =          |                  |
         *           0                       0               0(dim=1)
         *
         *   If this site is already a center, or has label "A" or "AC" then we absorb S as a new LC
         *
         *   0--[New "LC"]--1     =  0--[S]--1
         *
         *   Otherwise the S becomes the new "L" for this B site.
         *
         */
        auto U = other.unstash_U();
        if(U.dimension(0) != 1)
            throw std::logic_error(fmt::format("Failed to merge stash from right site {} into {}: "
                                               "U has invalid dimensions {}. Dim at idx 0 should be == 1",
                                               other.get_position(), get_position(), U.dimensions()));

        if(U.dimension(1) != get_chiR())
            throw std::logic_error(fmt::format("Failed to merge stash from right site {} into {}: "
                                               "M dimensions {} | U dimensions {} | Expected M(2) == U(1)",
                                               other.get_position(), get_position(), dimensions()));
        Textra::array3           dims = {spin_dim(), get_chiL(), U.dimension(2)};
        Eigen::Tensor<Scalar, 3> MU(dims);
        MU.device(Textra::omp::getDevice()) = get_M_bare().contract(U, Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(dims);

        set_M(MU);
        if(isCenter() or label.find('A') != std::string::npos) set_LC(other.unstash_S());
        else
            set_L(other.unstash_S());
    } else
        throw std::logic_error(fmt::format("Failed to merge stash: Wrong positions: this {} | other {}", get_position(), other.get_position()));
}

std::size_t class_mps_site::get_unique_id() const {
    if(unique_id) return unique_id.value();
    unique_id = hash::hash_buffer(get_M().data(), static_cast<size_t>(get_M().size()));
    unique_id = hash::hash_buffer(get_L().data(), static_cast<size_t>(get_L().size()), unique_id.value());
    //    if(isCenter()) unique_id = hash::hash_buffer(get_LC().data(), static_cast<size_t>(get_LC().size()), unique_id.value());
    return unique_id.value();
}
