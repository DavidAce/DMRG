//
// Created by david on 2017-11-13.
//

#include <config/nmspc_settings.h>
#include <iomanip>
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <tensors/state/class_mps_2site.h>
#include <tensors/state/class_mps_site.h>
#include <tools/common/log.h>
using namespace std;
using namespace Textra;

using Scalar = class_mps_2site::Scalar;

class_mps_2site::class_mps_2site(const class_mps_2site &other) : swapped(other.swapped), MPS_A(copy_unique(other.MPS_A)), MPS_B(copy_unique(other.MPS_B)) {}

void class_mps_2site::assert_validity() const {
    MPS_A->assert_validity();
    MPS_B->assert_validity();
}
bool class_mps_2site::is_real() const { return MPS_A->is_real() and MPS_B->is_real(); }
bool class_mps_2site::has_nan() const { return MPS_A->has_nan() or MPS_B->has_nan(); }
long class_mps_2site::chiA() const { return MPS_A->get_L().dimension(0); }
long class_mps_2site::chiB() const { return MPS_B->get_L().dimension(0); }
long class_mps_2site::chiC() const { return MPS_A->get_LC().dimension(0); }
long class_mps_2site::spin_dim_A() const { return MPS_A->spin_dim(); }
long class_mps_2site::spin_dim_B() const { return MPS_B->spin_dim(); }

const Eigen::Tensor<Scalar, 3> &class_mps_2site::A_bare() const { return MPS_A->get_M_bare(); }
const Eigen::Tensor<Scalar, 3> &class_mps_2site::A() const { return MPS_A->get_M(); }
const Eigen::Tensor<Scalar, 3> &class_mps_2site::B() const { return MPS_B->get_M(); }
Eigen::Tensor<Scalar, 2>        class_mps_2site::LC() const { return Textra::asDiagonal(MPS_A->get_LC()); }

Eigen::Tensor<Scalar, 3> class_mps_2site::GA() const {
    return MPS_A->get_M_bare().contract(Textra::asDiagonalInversed(MPS_A->get_L()), Textra::idx({1}, {1})).shuffle(Textra::array3{0, 2, 1});
}
Eigen::Tensor<Scalar, 3> class_mps_2site::GB() const { return MPS_B->get_M_bare().contract(Textra::asDiagonalInversed(MPS_B->get_L()), Textra::idx({2}, {0})); }
Eigen::Tensor<Scalar, 2> class_mps_2site::LA() const { return Textra::asDiagonal(MPS_A->get_L()); }
Eigen::Tensor<Scalar, 2> class_mps_2site::LB() const { return Textra::asDiagonal(MPS_B->get_L()); }

void class_mps_2site::set_mps(const Eigen::Tensor<Scalar, 3> &MA, const Eigen::Tensor<Scalar, 1> &LC, const Eigen::Tensor<Scalar, 3> &MB) {
    MPS_A->set_M(MA);
    MPS_A->set_LC(LC);
    MPS_B->set_M(MB);
}

void class_mps_2site::set_mps(const Eigen::Tensor<Scalar, 1> &LA, const Eigen::Tensor<Scalar, 3> &MA, const Eigen::Tensor<Scalar, 1> &LC,
                              const Eigen::Tensor<Scalar, 3> &MB, const Eigen::Tensor<Scalar, 1> &LB) {
    MPS_A->set_mps(MA, LA);
    MPS_A->set_LC(LC);
    MPS_B->set_mps(MB, LB);
}

Eigen::DSizes<long, 3> class_mps_2site::dimensions() const { return Eigen::DSizes<long, 3>{spin_dim_A()*spin_dim_B(), chiA(), chiB()}; }

void class_mps_2site::swap_AB() {
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
}

Eigen::Tensor<class_mps_2site::Scalar, 3> class_mps_2site::get_mps(Scalar norm) const
/*!
 * Returns a two-site MPS
     @verbatim
        1--[ LB ]--[ GA ]--[ LA ]-- [ GB ] -- [ LB ]--3
                     |                 |
                     0                 2
       Becomes

        1--[ mps ]--2
              |
              0

     @endverbatim
 */
{
    long dim0 = MPS_A->spin_dim() * MPS_B->spin_dim();
    long dim1 = MPS_A->get_chiL();
    long dim2 = MPS_B->get_chiR();
    return A().contract(B(), Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(Textra::array3{dim0, dim1, dim2}) / norm;
}
