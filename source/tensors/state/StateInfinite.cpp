#include "tensors/state/StateInfinite.h"
#include "config/settings.h"
#include "math/tenx.h"
#include "tensors/site/mps/MpsSite.h"
#include "tools/common/contraction.h"
#include "tools/common/log.h"
#include "tools/common/split.h"
#include "tools/common/views.h"
#include "tools/infinite/measure.h"
#include "tools/infinite/mps.h"
using Scalar = StateInfinite::Scalar;

StateInfinite::StateInfinite() : MPS_A(std::make_unique<MpsSite>()), MPS_B(std::make_unique<MpsSite>()) { tools::log->trace("Constructing state"); }

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
StateInfinite::~StateInfinite()                     = default;            // default dtor
StateInfinite::StateInfinite(StateInfinite &&other) = default;            // default move ctor
StateInfinite &StateInfinite::operator=(StateInfinite &&other) = default; // default move assign

/* clang-format off */
StateInfinite::StateInfinite(const StateInfinite &other):
    MPS_A(std::make_unique<MpsSite>(*other.MPS_A)),
    MPS_B(std::make_unique<MpsSite>(*other.MPS_B)),
    swapped(other.swapped),
    cache(other.cache),
    name(other.name),
    algo(other.algo),
    measurements(other.measurements),
    lowest_recorded_variance(other.lowest_recorded_variance){

}
/* clang-format on */

StateInfinite &StateInfinite::operator=(const StateInfinite &other) {
    // check for self-assignment
    if(this != &other) {
        MPS_A                    = std::make_unique<MpsSite>(*other.MPS_A);
        MPS_B                    = std::make_unique<MpsSite>(*other.MPS_B);
        swapped                  = other.swapped;
        cache                    = other.cache;
        name                     = other.name;
        algo                     = other.algo;
        measurements             = other.measurements;
        lowest_recorded_variance = other.lowest_recorded_variance;
    }
    return *this;
}

void StateInfinite::initialize(ModelType model_type) {
    tools::log->trace("Initializing state");
    long spin_dim = 2;
    switch(model_type) {
        case ModelType::ising_tf_rf: spin_dim = settings::model::ising_tf_rf::spin_dim; break;
        case ModelType::ising_sdual: spin_dim = settings::model::ising_sdual::spin_dim; break;
        default: spin_dim = 2;
    }
    Eigen::Tensor<Scalar, 3> M(spin_dim, 1, 1);
    Eigen::Tensor<Scalar, 1> L(1);
    // Default is a product state, spins pointing up in z.
    M.setZero();
    M(0, 0, 0) = 1;
    L.setConstant(1.0);
    MPS_A = std::make_unique<MpsSite>(M, L, 0, 0, "A");
    MPS_B = std::make_unique<MpsSite>(M, L, 1, 0, "B");
    MPS_A->set_LC(L);
}

void        StateInfinite::set_name(std::string_view statename) { name = statename; }
std::string StateInfinite::get_name() const { return name; }

void          StateInfinite::set_algorithm(const AlgorithmType &algo_type) { algo = algo_type; }
AlgorithmType StateInfinite::get_algorithm() const { return algo; }

std::pair<size_t, size_t> StateInfinite::get_positions() { return std::make_pair(MPS_A->get_position(), MPS_B->get_position()); }
size_t                    StateInfinite::get_positionA() { return MPS_A->get_position(); }
size_t                    StateInfinite::get_positionB() { return MPS_B->get_position(); }

long StateInfinite::chiC() const { return MPS_A->get_LC().dimension(0); }
long StateInfinite::chiA() const { return MPS_A->get_L().dimension(0); }
long StateInfinite::chiB() const { return MPS_B->get_L().dimension(0); }

long StateInfinite::get_spin_dimA() const { return MPS_A->spin_dim(); }
long StateInfinite::get_spin_dimB() const { return MPS_B->spin_dim(); }

double StateInfinite::get_truncation_error() const {
    // Should get the the current limit on allowed bond dimension for the duration of the simulation
    return get_mps_siteA().get_truncation_error();
}

Eigen::DSizes<long, 3> StateInfinite::dimensions() const { return Eigen::DSizes<long, 3>{get_spin_dimA() * get_spin_dimB(), chiA(), chiB()}; }

const MpsSite &StateInfinite::get_mps_siteA() const { return *MPS_A; }
const MpsSite &StateInfinite::get_mps_siteB() const { return *MPS_B; }
MpsSite       &StateInfinite::get_mps_siteA() { return *MPS_A; }
MpsSite       &StateInfinite::get_mps_siteB() { return *MPS_B; }
const MpsSite &StateInfinite::get_mps_site(size_t pos) const {
    if(pos == 0) return *MPS_A;
    if(pos == 1) return *MPS_B;
    throw except::runtime_error("Got wrong site position {}. Expected 0 or 1", pos);
}
MpsSite &StateInfinite::get_mps_site(size_t pos) {
    if(pos == 0) return *MPS_A;
    if(pos == 1) return *MPS_B;
    throw except::runtime_error("Got wrong site position {}. Expected 0 or 1", pos);
}
const MpsSite &StateInfinite::get_mps_site(std::string_view pos) const {
    if(pos == "A") return *MPS_A;
    if(pos == "B") return *MPS_B;
    throw except::runtime_error("Got wrong site position {}. Expected 0 or 1", pos);
}
MpsSite &StateInfinite::get_mps_site(std::string_view pos) {
    if(pos == "A") return *MPS_A;
    if(pos == "B") return *MPS_B;
    throw except::runtime_error("Got wrong site position {}. Expected 0 or 1", pos);
}

/* clang-format off */
const Eigen::Tensor<Scalar, 3> &StateInfinite::A_bare() const { return MPS_A->get_M_bare(); }
const Eigen::Tensor<Scalar, 3> &StateInfinite::A() const { return MPS_A->get_M(); }
const Eigen::Tensor<Scalar, 3> &StateInfinite::B() const { return MPS_B->get_M(); }
const Eigen::Tensor<Scalar, 2> &StateInfinite::LC_diag() const { if(cache.LC_diag) return cache.LC_diag.value(); else cache.LC_diag = tenx::asDiagonal(MPS_A->get_LC()); return cache.LC_diag.value();}
const Eigen::Tensor<Scalar, 2> &StateInfinite::LA_diag() const { if(cache.LA_diag) return cache.LA_diag.value(); else cache.LA_diag = tenx::asDiagonal(MPS_A->get_L()); return cache.LA_diag.value();}
const Eigen::Tensor<Scalar, 2> &StateInfinite::LB_diag() const { if(cache.LB_diag) return cache.LB_diag.value(); else cache.LB_diag = tenx::asDiagonal(MPS_B->get_L()); return cache.LB_diag.value();}
const Eigen::Tensor<Scalar, 2> &StateInfinite::LC_diag_inv() const { if(cache.LC_diag_inv) return cache.LC_diag_inv.value(); else cache.LC_diag_inv = tenx::asDiagonalInversed(MPS_A->get_LC()); return cache.LC_diag_inv.value();}
const Eigen::Tensor<Scalar, 2> &StateInfinite::LA_diag_inv() const { if(cache.LA_diag_inv) return cache.LA_diag_inv.value(); else cache.LA_diag_inv = tenx::asDiagonalInversed(MPS_A->get_L()); return cache.LA_diag_inv.value();}
const Eigen::Tensor<Scalar, 2> &StateInfinite::LB_diag_inv() const { if(cache.LB_diag_inv) return cache.LB_diag_inv.value(); else cache.LB_diag_inv = tenx::asDiagonalInversed(MPS_B->get_L()); return cache.LB_diag_inv.value();}
const Eigen::Tensor<Scalar, 1> &StateInfinite::LC() const { return MPS_A->get_LC(); }
const Eigen::Tensor<Scalar, 1> &StateInfinite::LA() const { return MPS_A->get_L(); }
const Eigen::Tensor<Scalar, 1> &StateInfinite::LB() const { return MPS_B->get_L(); }

const Eigen::Tensor<Scalar, 3> &StateInfinite::GA() const {
    if(cache.GA) return cache.GA.value();
    else {
        Eigen::Tensor<Scalar,1> L_inv = MPS_A->get_L().inverse();
        cache.GA = tools::common::contraction::contract_bnd_mps_temp(L_inv, MPS_A->get_M_bare());
    }
    return cache.GA.value();
}
const Eigen::Tensor<Scalar, 3> &StateInfinite::GB() const {
    if(cache.GB) return cache.GB.value();
    else{
        Eigen::Tensor<Scalar,1> L_inv = MPS_B->get_L().inverse();
        cache.GB = tools::common::contraction::contract_mps_bnd_temp(MPS_B->get_M_bare(), L_inv);
    }
    return cache.GB.value();
}

/* clang-format on */

const Eigen::Tensor<Scalar, 3> &StateInfinite::get_2site_mps(Scalar norm) const {
    /*!
 * Returns a two-site tensor
     @verbatim
        1--[ LA ]--[ GA ]--[ LC ]-- [ GB ] -- [ LB ]--3
                     |                 |
                     0                 2
       which in our notation is simply A * B  (because A = LA * GA * LC and B = GB * LB),
       becomes

        1--[ 2site_tensor ]--2
                  |
                  0

     @endverbatim
 */

    if(cache.twosite_mps) return cache.twosite_mps.value();
    cache.twosite_mps = tools::common::contraction::contract_mps_mps_temp(A(), B()) / norm;
    return cache.twosite_mps.value();
}

void StateInfinite::assert_validity() const {
    MPS_A->assert_validity();
    MPS_B->assert_validity();
}
bool StateInfinite::is_real() const { return MPS_A->is_real() and MPS_B->is_real(); }
bool StateInfinite::has_nan() const { return MPS_A->has_nan() or MPS_B->has_nan(); }
//
// template<typename T>
// Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> StateInfinite::get_H_local_matrix() const {
//    Eigen::Tensor<T, 5> tempL;
//    Eigen::Tensor<T, 5> tempR;
//    if constexpr(std::is_same<T, double>::value) {
//        if(not Lblock->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from Lblock when building H_local");
//        }
//        if(not Rblock->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from Rblock when building H_local");
//        }
//        if(not HA->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from MPO A when building H_local");
//        }
//        if(not HB->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from MPO B when building H_local");
//        }
//        tempL = Lblock->block.contract(HA->MPO(), tenx::idx({2}, {0})).real().shuffle(tenx::array5{4, 1, 3, 0, 2}).real();
//        tempR = Rblock->block.contract(HB->MPO(), tenx::idx({2}, {1})).real().shuffle(tenx::array5{4, 1, 3, 0, 2}).real();
//    } else {
//        tempL = Lblock->block.contract(HA->MPO(), tenx::idx({2}, {0})).shuffle(tenx::array5{4, 1, 3, 0, 2});
//        tempR = Rblock->block.contract(HB->MPO(), tenx::idx({2}, {1})).shuffle(tenx::array5{4, 1, 3, 0, 2});
//    }
//    long                shape   = mps_sites->chiA() * mps_sites->spindim() * mps_sites->chiB() * mps_sites->spin_dim_A();
//    Eigen::Tensor<T, 8> H_local = tempL.contract(tempR, tenx::idx({4}, {4})).shuffle(tenx::array8{0, 1, 4, 5, 2, 3, 6, 7});
//    return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(H_local.data(), shape, shape);
//}

// template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>               StateInfinite::get_H_local_matrix<double>() const;
// template Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> StateInfinite::get_H_local_matrix<std::complex<double>>() const;
//
// template<typename T>
// Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> StateInfinite::get_H_local_sq_matrix() const {
//    Eigen::Tensor<T, 6> tempL;
//    Eigen::Tensor<T, 6> tempR;
//    if constexpr(std::is_same<T, double>::value) {
//        if(not Lblock2->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from Lblock2 when building H_local_sq");
//        }
//        if(not Rblock2->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from Rblock2 when building H_local_sq");
//        }
//        if(not HA->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from MPO A when building H_local_sq");
//        }
//        if(not HB->isReal()) {
//            throw std::runtime_error("Discarding imaginary data from MPO B when building H_local_sq");
//        }
//        tempL = Lblock2->block.contract(HA->MPO(), tenx::idx({2}, {0}))
//                    .contract(HA->MPO(), tenx::idx({2, 5}, {0, 2}))
//                    .real()
//                    .shuffle(tenx::array6{5, 1, 3, 0, 2, 4});
//        tempR = Rblock2->block.contract(HB->MPO(), tenx::idx({2}, {1}))
//                    .contract(HB->MPO(), tenx::idx({2, 5}, {1, 2}))
//                    .real()
//                    .shuffle(tenx::array6{5, 1, 3, 0, 2, 4});
//    } else {
//        tempL = Lblock2->block.contract(HA->MPO(), tenx::idx({2}, {0}))
//                    .contract(HA->MPO(), tenx::idx({2, 5}, {0, 2}))
//                    .shuffle(tenx::array6{5, 1, 3, 0, 2, 4});
//        tempR = Rblock2->block.contract(HB->MPO(), tenx::idx({2}, {1}))
//                    .contract(HB->MPO(), tenx::idx({2, 5}, {1, 2}))
//                    .shuffle(tenx::array6{5, 1, 3, 0, 2, 4});
//    }
//
//    long                shape   = mps_sites->chiA() * mps_sites->spindim() * mps_sites->chiB() * mps_sites->spin_dim_A();
//    Eigen::Tensor<T, 8> H_local = tempL.contract(tempR, tenx::idx({4, 5}, {4, 5})).shuffle(tenx::array8{0, 1, 4, 5, 2, 3, 6, 7});
//    return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(H_local.data(), shape, shape);
//}

// template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>               StateInfinite::get_H_local_sq_matrix<double>() const;
// template Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> StateInfinite::get_H_local_sq_matrix<std::complex<double>>() const;

// void StateInfinite::enlarge_environment(int direction) {
//    assert_positions();
//    auto position_A_new = mps_sites->MPS_A->get_position() + 1;
//    *Lblock             = Lblock->enlarge(*mps_sites->MPS_A, *HA);
//    *Rblock             = Rblock->enlarge(*mps_sites->MPS_B, *HB);
//    *Lblock2            = Lblock2->enlarge(*mps_sites->MPS_A, *HA);
//    *Rblock2            = Rblock2->enlarge(*mps_sites->MPS_B, *HB);
//    set_positions(position_A_new);
//    if(direction != 0) throw std::runtime_error("Ooops, direction != 0");
//    if (direction == 1){
//        *Lblock  = Lblock->enlarge(*mps_sites->MPS_A,  *HA);
//        *Lblock2 = Lblock2->enlarge(*mps_sites->MPS_A, *HA);
//        Lblock->set_position (HB->get_position());
//        Lblock2->set_position(HB->get_position());
//    }else if (direction == -1){
//        *Rblock  = Rblock->enlarge(*mps_sites->MPS_B,  *HB);
//        *Rblock2 = Rblock2->enlarge(*mps_sites->MPS_B, *HB);
//        Rblock->set_position (HA->get_position());
//        Rblock2->set_position(HA->get_position());
//    }else if(direction == 0){
//
//
//    }
//    clear_measurements();
//    assert_positions();
//}

void StateInfinite::set_positions(size_t position) {
    MPS_A->set_position(position);
    MPS_B->set_position(position + 1);
}

void StateInfinite::set_mps(const Eigen::Tensor<Scalar, 3> &twosite_tensor, std::optional<svd::config> svd_cfg) {
    tools::infinite::mps::merge_twosite_tensor(*this, twosite_tensor, svd_cfg);
}

void StateInfinite::set_mps(const std::vector<MpsSite> &mps_list) {
    if(mps_list.size() != 2) throw except::runtime_error("Expected 2 sites, got: {}", mps_list.size());
    const auto &mpsA = *std::next(mps_list.begin(), 0);
    const auto &mpsB = *std::next(mps_list.begin(), 1);
    set_mps(mpsA, mpsB);
    clear_cache();
}

void StateInfinite::set_mps(const MpsSite &mpsA, const MpsSite &mpsB) {
    if(not mpsA.isCenter()) throw std::runtime_error("Given mps for site A is not a center");
    MPS_A = std::make_unique<MpsSite>(mpsA);
    MPS_B = std::make_unique<MpsSite>(mpsB);
    clear_cache();
}

void StateInfinite::set_mps(const Eigen::Tensor<Scalar, 3> &MA, const Eigen::Tensor<Scalar, 1> &LC, const Eigen::Tensor<Scalar, 3> &MB) {
    MPS_A->set_M(MA);
    MPS_A->set_LC(LC);
    MPS_B->set_M(MB);
    clear_cache();
}
void StateInfinite::set_mps(const Eigen::Tensor<Scalar, 1> &LA, const Eigen::Tensor<Scalar, 3> &MA, const Eigen::Tensor<Scalar, 1> &LC,
                            const Eigen::Tensor<Scalar, 3> &MB, const Eigen::Tensor<Scalar, 1> &LB) {
    MPS_A->set_mps(MA, LA, 0, "A");
    MPS_A->set_LC(LC);
    MPS_B->set_mps(MB, LB, 0, "B");
    clear_cache();
}

bool StateInfinite::is_bond_limited(long bond_lim, double truncation_threshold) const {
    bool has_truncated   = get_truncation_error() > truncation_threshold;
    bool has_reached_max = chiC() >= bond_lim;
    return has_truncated or has_reached_max;
}

void StateInfinite::clear_cache() const { cache = Cache(); }

void StateInfinite::clear_measurements() const {
    measurements                              = MeasurementsStateInfinite();
    tools::common::views::components_computed = false;
}

void StateInfinite::do_all_measurements() const { tools::infinite::measure::do_all_measurements(*this); }

void StateInfinite::swap_AB() {
    tools::log->trace("Swapping AB");
    swapped = !swapped;
    // Store the positions
    auto position_left  = MPS_A->get_position();
    auto position_right = MPS_B->get_position();

    // Swap Gamma
    Eigen::Tensor<Scalar, 1> LC = MPS_A->get_LC();
    MPS_A->unset_LC();
    MPS_B->unset_LC();
    MPS_A.swap(MPS_B);
    MPS_A->set_LC(MPS_A->get_L());
    MPS_A->set_L(LC);
    MPS_B->set_L(LC);
    MPS_A->set_position(position_left);
    MPS_B->set_position(position_right);
    clear_cache();
    clear_measurements();
}
