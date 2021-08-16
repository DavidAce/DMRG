#include <math/tenx.h>
// tenx must appear first
#include "MpsSite.h"
#include <config/debug.h>
#include <math/hash.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <utility>

using Scalar = MpsSite::Scalar;

MpsSite::MpsSite() = default;
MpsSite::MpsSite(const Eigen::Tensor<Scalar, 3> &M_, const Eigen::Tensor<Scalar, 1> &L_, size_t pos, double error, std::string_view label_)
//    : M(M_), L(L_), position(pos), truncation_error(error), label(std::move(label_))
{
    set_position(pos);
    set_label(label_);
    set_M(M_);
    set_L(L_);
    set_truncation_error(error);
}

MpsSite::MpsSite(const Eigen::Tensor<Scalar, 3> &M_, std::optional<Eigen::Tensor<Scalar, 1>> L_, size_t pos, double error, std::string_view label_)
//    : M(M_), L(std::move(L_)), position(pos), truncation_error(error), label(std::move(label_))
{
    set_position(pos);
    set_label(label_);
    set_M(M_);
    if(L_) set_L(L_.value());
    set_truncation_error(error);
}

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
MpsSite::~MpsSite()               = default;            // default dtor
MpsSite::MpsSite(MpsSite &&other) = default;            // default move ctor
MpsSite &MpsSite::operator=(MpsSite &&other) = default; // default move assign
MpsSite::MpsSite(const MpsSite &other)       = default;
MpsSite &MpsSite::operator=(const MpsSite &other) = default;

bool MpsSite::isCenter() const { return LC.has_value(); }

bool MpsSite::has_L() const { return L.has_value(); }
bool MpsSite::has_M() const { return M.has_value(); }
bool MpsSite::has_LC() const { return LC.has_value(); }

Eigen::DSizes<long, 3> MpsSite::dimensions() const { return Eigen::DSizes<long, 3>{spin_dim(), get_chiL(), get_chiR()}; }

bool MpsSite::is_real() const { return tenx::isReal(get_M_bare(), "M_bare") and tenx::isReal(get_L(), "L"); }
bool MpsSite::has_nan() const { return tenx::hasNaN(get_M_bare(), "M_bare") and tenx::hasNaN(get_L(), "L"); }
void MpsSite::assert_validity() const {
    if(has_nan()) throw std::runtime_error(fmt::format("MpsSite::assert_validity(): MPS (M or L) at position {} has NaN's", get_position()));
}

void MpsSite::assert_dimensions() const {
    if(get_label() == "B" and get_chiR() != get_L().dimension(0))
        throw std::runtime_error(fmt::format("MpsSite: Assert failed: Dimensions for B and L are incompatible at site {}: B {} | L {}", get_position(),
                                             get_M_bare().dimensions(), get_L().dimensions()));
    if(get_label() == "A" and get_chiL() != get_L().dimension(0))
        throw std::runtime_error(fmt::format("MpsSite: Assert failed: Dimensions for A and L are incompatible at site {}: A {} | L {}", get_position(),
                                             get_M_bare().dimensions(), get_L().dimensions()));
    if(get_label() == "AC") {
        if(get_chiL() != get_L().dimension(0))
            throw std::runtime_error(fmt::format("MpsSite: Assert failed: Dimensions for AC and L are incompatible at site {}: AC {} | L {}", get_position(),
                                                 get_M_bare().dimensions(), get_L().dimensions()));
        if(get_chiR() != get_LC().dimension(0))
            throw std::runtime_error(fmt::format("MpsSite: Assert failed: Dimensions for AC and L are incompatible at site {}: AC {} | LC{}", get_position(),
                                                 get_M_bare().dimensions(), get_L().dimensions()));
    }
}

void MpsSite::assert_identity() const {
    auto t_dbg = tid::tic_token("assert_identity");
    if(get_label() == "B") {
        Eigen::Tensor<Scalar, 2> id = get_M_bare().contract(get_M_bare().conjugate(), tenx::idx({0, 2}, {0, 2}));
        if(not tenx::MatrixMap(id).isIdentity(1e-10)) {
            throw std::runtime_error(fmt::format("MpsSite: {0}^dagger {0} is not identity at pos {1}: \n{2}", get_label(), get_position(), get_M_bare()));
        }
    } else {
        Eigen::Tensor<Scalar, 2> id = get_M_bare().contract(get_M_bare().conjugate(), tenx::idx({0, 1}, {0, 1}));
        if(not tenx::MatrixMap(id).isIdentity(1e-10)) {
            throw std::runtime_error(fmt::format("MpsSite: {0}^dagger {0} is not identity at pos {1}: \n{2}", get_label(), get_position(), get_M_bare()));
        }
    }
    if(isCenter() or get_label() == "AC") {
        Eigen::Tensor<Scalar, 0> MM   = get_M().contract(get_M().conjugate(), tenx::idx({0, 1, 2}, {0, 1, 2}));
        auto                     norm = std::real(MM(0));
        if(std::abs(norm - 1) > 1e-10)
            throw std::runtime_error(fmt::format("MpsSite: {0}^dagger {0} is not unity at pos {1}: {2:.16f}", get_label(), get_position(), norm));
    }
}

const Eigen::Tensor<Scalar, 3> &MpsSite::get_M_bare() const {
    if(M) {
        if(M.value().size() == 0) throw std::runtime_error(fmt::format("MpsSite::get_M_bare(): M has size 0 at position {}", get_position()));
        return M.value();
    } else
        throw std::runtime_error(fmt::format("MpsSite::get_M_bare(): M has not been set at position {}", get_position()));
}
const Eigen::Tensor<Scalar, 3> &MpsSite::get_M() const {
    if(isCenter()) {
        if(LC.value().dimension(0) != get_M_bare().dimension(2))
            throw std::runtime_error(fmt::format("MpsSite::get_M(): M and LC dim mismatch: {} != {} at position {}", get_M_bare().dimension(2),
                                                 LC.value().dimension(0), get_position()));
        if(MC) {
            if(MC.value().size() == 0) throw std::runtime_error(fmt::format("MpsSite::get_M(): MC has size 0 at position {}", get_position()));
            return MC.value();
        } else {
            MC                                 = Eigen::Tensor<Scalar, 3>(spin_dim(), get_chiL(), get_chiR());
            MC->device(tenx::omp::getDevice()) = get_M_bare().contract(tenx::asDiagonal(get_LC()), tenx::idx({2}, {0}));
            if(MC->size() == 0) throw std::runtime_error(fmt::format("MpsSite::get_M(): built MC with size 0 at position {}", get_position()));
            return MC.value();
        }
    } else
        return get_M_bare();
}

const Eigen::Tensor<Scalar, 1> &MpsSite::get_L() const {
    if(L) {
        if constexpr(settings::debug) {
            if(L->size() == 0) throw std::runtime_error(fmt::format("MpsSite::get_L(): L has size 0 at position {} | label {}", get_position(), get_label()));
            if(get_label() != "B" and L->size() != get_chiL())
                throw std::runtime_error(fmt::format("MpsSite::get_L(): L.size() {} != chiL {} at position {} | M dims {} | label {}", L->size(), get_chiL(),
                                                     get_position(), dimensions(), get_label()));
            if(get_label() == "B" and L->size() != get_chiR())
                throw std::runtime_error(fmt::format("MpsSite::get_L(): L.size() {} != chiR {} at position {} | M dims {} | label {}", L->size(), get_chiR(),
                                                     get_position(), dimensions(), get_label()));
        }
        return L.value();
    } else
        throw std::runtime_error(fmt::format("MpsSite::get_L(): L has not been set at position {}", get_position()));
}
const Eigen::Tensor<Scalar, 1> &MpsSite::get_LC() const {
    if(isCenter()) {
        if(LC->dimension(0) != get_M_bare().dimension(2))
            throw std::runtime_error(fmt::format("MpsSite::get_LC(): M dimensions {} are incompatible with LC dimensions {} at position {}",
                                                 get_M_bare().dimensions(), LC.value().dimensions(), get_position()));
        if(LC.value().size() == 0) throw std::runtime_error(fmt::format("MpsSite::get_LC(): LC has size 0 at position {}", get_position()));
        return LC.value();
    } else
        throw std::runtime_error(fmt::format("MpsSite::get_LC(): Site at position {} is not a center", get_position()));
}

Eigen::Tensor<Scalar, 3> &MpsSite::get_M_bare() { return const_cast<Eigen::Tensor<Scalar, 3> &>(std::as_const(*this).get_M_bare()); }
Eigen::Tensor<Scalar, 3> &MpsSite::get_M() { return const_cast<Eigen::Tensor<Scalar, 3> &>(std::as_const(*this).get_M()); }
Eigen::Tensor<Scalar, 1> &MpsSite::get_L() { return const_cast<Eigen::Tensor<Scalar, 1> &>(std::as_const(*this).get_L()); }
Eigen::Tensor<Scalar, 1> &MpsSite::get_LC() { return const_cast<Eigen::Tensor<Scalar, 1> &>(std::as_const(*this).get_LC()); }
double                    MpsSite::get_truncation_error() const { return truncation_error; }
double                    MpsSite::get_truncation_error_LC() const { return truncation_error_LC; }
std::string_view          MpsSite::get_label() const {
    if(label.empty()) throw std::runtime_error(fmt::format("No label found at position {}", get_position()));
    return label;
}
std::tuple<long, long, long> MpsSite::get_dims() const { return {spin_dim(), get_chiL(), get_chiR()}; }
long                         MpsSite::spin_dim() const { return get_M_bare().dimension(0); }
long                         MpsSite::get_chiL() const { return get_M_bare().dimension(1); }
long                         MpsSite::get_chiR() const { return get_M_bare().dimension(2); }

void MpsSite::set_position(const size_t position_) {
    position = position_;
    MC.reset();
}

template<typename T>
T MpsSite::get_position() const {
    if(position) {
        return static_cast<T>(position.value());
    } else {
        throw std::runtime_error("MpsSite::get_position(): Position hasn't been set on mps site.");
    }
}
template size_t MpsSite::get_position<size_t>() const;
template long   MpsSite::get_position<long>() const;

void MpsSite::set_mps(const Eigen::Tensor<Scalar, 3> &M_, const Eigen::Tensor<Scalar, 1> &L_, double error, std::string_view label_) {
    // M has to be a "bare" matrix, i.e. not an MC which would include LC.
    set_M(M_);
    set_L(L_);
    set_truncation_error(error);
    set_label(label_);
}

void MpsSite::set_M(const Eigen::Tensor<Scalar, 3> &M_) {
    // M has to be a "bare" matrix, i.e. not an MC which would include LC.
    if(position) {
        M = M_;
        MC.reset();
        unique_id = std::nullopt;
    } else
        throw std::runtime_error("MpsSite::set_M(const Eigen::Tensor<Scalar, 3> &): Can't set M: Position hasn't been set yet");
}
void MpsSite::set_L(const Eigen::Tensor<Scalar, 1> &L_, double error) {
    if constexpr(settings::debug) {
        auto norm = tenx::VectorMap(L_).norm();
        if(std::abs(norm - 1) > 1e-8) tools::log->warn("MpsSite::set_L(): Norm of L is too far from unity: {:.16f}", norm);
        //            throw std::runtime_error(fmt::format("MpsSite::set_L(): Can't set L: Norm of L is too far from unity: {:.16f}", norm));
        //        if(position and position.value() == 0 and get_label() != "B"){
        //            if(L_.size() != 1) throw std::logic_error("Left edge L should have size 1");
        //            if(std::abs(L_.coeff(0)) != 1.0) throw std::logic_error("Left edge L should be equal to 1");
        //        }
    }

    if(position) {
        L                = L_;
        truncation_error = error;
        unique_id        = std::nullopt;
    } else
        throw std::runtime_error("Can't set L: Position hasn't been set yet");
}
void MpsSite::set_L(const std::pair<Eigen::Tensor<Scalar, 1>, double> &L_and_error) { set_L(L_and_error.first, L_and_error.second); }

void MpsSite::set_LC(const Eigen::Tensor<Scalar, 1> &LC_, double error) {
    if constexpr(settings::debug) {
        auto norm = tenx::VectorMap(LC_).norm();
        if(std::abs(norm - 1) > 1e-8) tools::log->warn("MpsSite::set_LC(): Norm of LC is too far from unity: {:.16f}", norm);
        //            throw std::runtime_error(fmt::format("MpsSite::set_LC(): Can't set L: Norm of LC is too far from unity: {:.16f}", norm));
        if(position and position.value() == 0 and get_label() == "B") {
            throw std::logic_error("Only an A-site can become an AC site (not really true though");
        }
    }
    if(position) {
        LC = LC_;
        MC.reset();
        truncation_error_LC = error;
        unique_id           = std::nullopt;
        set_label("AC");
    } else
        throw std::runtime_error("Can't set LC: Position hasn't been set yet");
}

void MpsSite::set_LC(const std::pair<Eigen::Tensor<Scalar, 1>, double> &LC_and_error) { set_LC(LC_and_error.first, LC_and_error.second); }

void MpsSite::set_truncation_error(double error) { truncation_error = error; }
void MpsSite::set_truncation_error_LC(double error) { truncation_error_LC = error; }
void MpsSite::set_label(std::string_view label_) {
    // If we are flipping the kind of site from B to non B (or vice-versa), then the L matrix
    // should be removed since it changes side.
    if(not label.empty() and not label_.empty() and label.front() != label_.front()) unset_L();
    label = label_;
}
void MpsSite::unset_LC() {
    LC        = std::nullopt;
    MC        = std::nullopt;
    unique_id = std::nullopt;
    if(label == "AC") label = "A";
}

void MpsSite::unset_L() {
    L         = std::nullopt;
    unique_id = std::nullopt;
}

void MpsSite::fuse_mps(const MpsSite &other) {
    auto t_fuse = tid::tic_scope("fuse");
    // This operation is done when merging mps after an svd split, for instance
    if(get_position() != other.get_position())
        throw std::runtime_error(
            fmt::format("MpsSite::fuse_mps(const MpsSite &): position mismatch: this pos {} != other pos {}", get_position(), other.get_position()));

    if(not other.has_M()) throw std::runtime_error("MpsSite::fuse_mps(const MpsSite &): Got mps site with undefined M");
    if constexpr(settings::debug) {
        if(other.has_L())
            tools::log->trace("Merging {}[{}] | M dims {} | L dim {}", other.get_label(), other.get_position(), other.dimensions(), other.get_L().dimensions());
        else
            tools::log->trace("{}[{}] | M dims {}", other.get_label(), other.get_position(), other.dimensions());
    }

    // We have to copy the bare "M", i.e. not MC, which would include LC.
    set_M(other.get_M_bare());

    // Copy the A/AC/B label
    set_label(other.get_label());

    // If we are flipping the kind of site from B to non B (or vice versa), then the L matrix
    // should be removed since it changes side.
    if(get_label().front() != other.get_label().front()) unset_L();

    // We should only copy the L matrix if it exists
    // If it does not exist it may just not have been set.
    // For instance, when splitting a multisite tensor into its contituent mps sites,
    // the edge mps sites have empty L matrices on them because these are not updated.
    // This is intended, but should be checked for.
    // Note that if there is no L on the other mps, there has to exist a previous one
    // on this mps -- we can't leave here without having defined L.

    if(other.has_L())
        set_L(other.get_L(), other.get_truncation_error());
    else if(not has_L()) {
        // Now we have a situation where no L was given, and no L is found on this site.
        // For edge-sites this is simple: then the mps has edge-dimension == 1, so the only
        // way to have a normalized L is to set it to 1. Otherwise, this is a failure.
        auto one = Eigen::Tensor<Scalar, 1>(1).setConstant(1.0);
        if((label != "B" and get_chiL() == 1) or (label == "B" and get_chiR() == 1))
            set_L(one);
        else
            throw std::runtime_error(fmt::format("MpsSite::fuse_mps(const MpsSite &): position {} | label {} | "
                                                 "Got other mps without an L, and there is no preexisting L",
                                                 get_position(), get_label()));
    }

    // Copy the center LC if it's in other
    if(other.isCenter()) {
        set_LC(other.get_LC(), other.get_truncation_error_LC());
        if(get_M_bare().dimension(2) != get_LC().dimension(0))
            throw std::runtime_error(
                fmt::format("Got center of wrong dimension. M dims {} are not compatible with LC dims {}", get_M_bare().dimensions(), get_LC().dimensions()));
    }
}

void MpsSite::apply_mpo(const Eigen::Tensor<Scalar, 4> &mpo) {
    auto t_mpo = tid::tic_token("apply_mpo");
    tools::log->trace("Applying mpo (dims {}) at position {} | isCenter: {}", mpo.dimensions(), get_position(), isCenter());
    long mpoDimL = mpo.dimension(0);
    long mpoDimR = mpo.dimension(1);
    if(mpoDimL != mpoDimR) throw std::logic_error("Can't apply mpo's with different L/R dims: not implemented yet");

    if(isCenter()) {
        Eigen::Tensor<Scalar, 1> LC_temp = get_LC().broadcast(tenx::array1{mpoDimR});
        tenx::normalize(LC_temp);
        set_LC(LC_temp);
    }
    Eigen::Tensor<Scalar, 3> M_bare_temp(tenx::array3{spin_dim(), get_chiL() * mpoDimL, get_chiR() * mpoDimR});
    M_bare_temp.device(tenx::omp::getDevice()) = get_M_bare()
                                                     .contract(mpo, tenx::idx({0}, {2}))
                                                     .shuffle(tenx::array5{4, 0, 2, 1, 3})
                                                     .reshape(tenx::array3{spin_dim(), get_chiL() * mpoDimL, get_chiR() * mpoDimR});
    Eigen::Tensor<Scalar, 1> L_temp = get_L().broadcast(tenx::array1{mpoDimL});
    tenx::normalize(L_temp);
    set_L(L_temp);
    set_M(M_bare_temp);
}

void MpsSite::apply_mpo(const Eigen::Tensor<Scalar, 2> &mpo) {
    auto                     t_mpo = tid::tic_token("apply_mpo");
    Eigen::Tensor<Scalar, 3> M_bare_temp(mpo.dimension(1), get_chiL(), get_chiR());
    M_bare_temp.device(tenx::omp::getDevice()) = mpo.contract(get_M_bare(), tenx::idx({0}, {0}));
    set_M(M_bare_temp);
}

void MpsSite::stash_U(const Eigen::Tensor<Scalar, 3> &U, size_t dst) const { U_stash = stash<Eigen::Tensor<Scalar, 3>>{U, 0, dst}; }
void MpsSite::stash_S(const Eigen::Tensor<Scalar, 1> &S, double error, size_t dst) const { S_stash = stash<Eigen::Tensor<Scalar, 1>>{S, error, dst}; }
void MpsSite::stash_S(const std::pair<Eigen::Tensor<Scalar, 1>, double> &S_and_error, size_t dst) const {
    S_stash = stash<Eigen::Tensor<Scalar, 1>>{S_and_error.first, S_and_error.second, dst};
}
void MpsSite::stash_C(const Eigen::Tensor<Scalar, 1> &C, double error, size_t dst) const { C_stash = stash<Eigen::Tensor<Scalar, 1>>{C, error, dst}; }
void MpsSite::stash_C(const std::pair<Eigen::Tensor<Scalar, 1>, double> &C_and_error, size_t dst) const {
    C_stash = stash<Eigen::Tensor<Scalar, 1>>{C_and_error.first, C_and_error.second, dst};
}
void MpsSite::stash_V(const Eigen::Tensor<Scalar, 3> &V, size_t dst) const { V_stash = stash<Eigen::Tensor<Scalar, 3>>{V, 0, dst}; }

std::optional<stash<Eigen::Tensor<Scalar, 3>>> &MpsSite::get_U_stash() const { return U_stash; }
std::optional<stash<Eigen::Tensor<Scalar, 1>>> &MpsSite::get_S_stash() const { return S_stash; }
std::optional<stash<Eigen::Tensor<Scalar, 1>>> &MpsSite::get_C_stash() const { return C_stash; }
std::optional<stash<Eigen::Tensor<Scalar, 3>>> &MpsSite::get_V_stash() const { return V_stash; }

void MpsSite::drop_stash() const {
    if constexpr(settings::debug) {
        if(U_stash) tools::log->trace("Dropping U_stash");
        if(S_stash) tools::log->trace("Dropping S_stash");
        if(C_stash) tools::log->trace("Dropping C_stash");
        if(V_stash) tools::log->trace("Dropping V_stash");
    }
    U_stash = std::nullopt;
    S_stash = std::nullopt;
    C_stash = std::nullopt;
    V_stash = std::nullopt;
}

void MpsSite::take_stash(const MpsSite &other) {
    auto t_stash = tid::tic_token("take_stash");
    if(other.V_stash and other.V_stash->pos_dst == get_position()) {
        /* Left-to-right move.
         * In this case there is a "V" in "other" that should be absorbed into this site from the left.
         *
         *   1 --[New "M"]-- 2       1 --[V]-- 2   1 --[this "M"]-- 2
         *          |            =        |                |
         *          0                 0(dim=1)             0
         *
         *  In addition, if this is an A site we expect an "S" to come along, which didn't become an LC
         *  for the site on the left. Presumably the true LC is on some site further to the right.
         *  Here we simply set it as the new L of this site.
         */
        if constexpr(settings::debug) tools::log->trace("Merging V stash from site {} into {}", other.get_position(), get_position());

        if(get_position() != other.get_position() + 1)
            throw std::logic_error(fmt::format("Found V stash at the wrong position: This {} | Other {} | Stash destination {}", get_position(),
                                               other.get_position(), other.V_stash->pos_dst));

        auto &V = other.V_stash->data;
        if(V.dimension(0) != 1)
            throw std::logic_error(fmt::format("Failed to take stash from left site {} into {}: "
                                               "V has invalid dimensions {}. Dim at idx 0 should be == 1",
                                               other.get_position(), get_position(), V.dimensions()));
        if(V.dimension(2) != get_chiL())
            throw std::logic_error(fmt::format("Failed to take stash from left site {} into {}: "
                                               "V dimensions {} | M dimensions {} |  Expected V(2) == M(1)",
                                               other.get_position(), get_position(), V.dimensions(), dimensions()));
        tenx::array3             dims = {spin_dim(), V.dimension(1), get_chiR()};
        Eigen::Tensor<Scalar, 3> VM(dims);
        VM.device(tenx::omp::getDevice()) = V.contract(get_M_bare(), tenx::idx({2}, {1})).shuffle(tenx::array4{0, 2, 1, 3}).reshape(dims);
        set_M(VM);
        other.V_stash = std::nullopt; // Stash has been consumed
    }
    if(other.U_stash and other.U_stash->pos_dst == get_position()) {
        /* Right-to-left move.
         * In this case there should be a U and an S in "other".
         * U should be absorbed into this site from the right.
         * S becomes the new LC.
         *
         *    1 --[New "M"]-- 2       1 --[this "M"]-- 2   1 --[U]-- 2
         *           |            =          |                  |
         *           0                       0               0(dim=1)
         *
         */

        if constexpr(settings::debug) tools::log->trace("Merging U stash from site {} into {}", other.get_position(), get_position());

        if(get_position() != other.get_position() - 1)
            throw std::logic_error(fmt::format("Found V stash at the wrong position: This {} | Other {} | Stash destination {}", get_position(),
                                               other.get_position(), other.U_stash->pos_dst));

        auto &U = other.U_stash->data;
        if(U.dimension(0) != 1)
            throw std::logic_error(fmt::format("Failed to take stash from right site {} into {}: "
                                               "U has invalid dimensions {}. Dim at idx 0 should be == 1",
                                               other.get_position(), get_position(), U.dimensions()));

        if(U.dimension(1) != get_chiR())
            throw std::logic_error(fmt::format("Failed to take stash from right site {} into {}: "
                                               "M dimensions {} | U dimensions {} | Expected M(2) == U(1)",
                                               other.get_position(), get_position(), dimensions(), U.dimensions()));
        tenx::array3             dims = {spin_dim(), get_chiL(), U.dimension(2)};
        Eigen::Tensor<Scalar, 3> MU(dims);
        MU.device(tenx::omp::getDevice()) = get_M_bare().contract(U, tenx::idx({2}, {1})).shuffle(tenx::array4{0, 2, 1, 3}).reshape(dims);
        set_M(MU);
        other.U_stash = std::nullopt;
    }
    if(other.S_stash and other.S_stash->pos_dst == get_position()) {
        /* We get this stash when
         *      - This is the right-most last A-site and AC is further to the right: Then we get both S and V to absorb on this site.
         *      - This is being transformed from AC to a B-site. Then the old LC matrix is inherited as an L matrix.
         */
        if constexpr(settings::debug) tools::log->trace("Merging S stash from site {} into {}", other.get_position(), get_position());
        set_L(other.S_stash->data, other.S_stash->error);
        other.S_stash = std::nullopt;
    }
    if(other.C_stash and other.C_stash->pos_dst == get_position()) {
        /* We get this stash when
         *      - This is the right-most last A-site and supposed to become an AC: Then we get S, V and LC to absorb on this site.
         *      - This is being transformed from a B to an AC-site. Then the LC was just created in an SVD.
         */
        if constexpr(settings::debug) tools::log->trace("Merging C stash from site {} into {}", other.get_position(), get_position());
        if(label == "B") tools::log->warn("Setting LC at pos {} with label {}", get_position(), get_label());
        set_LC(other.C_stash->data, other.C_stash->error);
        other.C_stash = std::nullopt;
    }
}

std::size_t MpsSite::get_unique_id() const {
    if(unique_id) return unique_id.value();
    unique_id = hash::hash_buffer(get_M().data(), static_cast<size_t>(get_M().size()));
    unique_id = hash::hash_buffer(get_L().data(), static_cast<size_t>(get_L().size()), unique_id.value());
    //    if(isCenter()) unique_id = hash::hash_buffer(get_LC().data(), static_cast<size_t>(get_LC().size()), unique_id.value());
    return unique_id.value();
}