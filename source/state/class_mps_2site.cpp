//
// Created by david on 2017-11-13.
//


#include <state/class_mps_2site.h>
#include <state/class_mps_site.h>
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <simulation/nmspc_settings.h>
#include <iomanip>
#include <tools/common/log.h>
using namespace std;
using namespace Textra;

using Scalar = class_mps_2site::Scalar;

// Function definitions

//class_mps_2site::class_mps_2site(){
//    tools::log->trace("Constructing 2site MPS");

//size_t spin_dim = 2;
//if(model_type_str == "ising_tf_rf_nn")
//spin_dim = settings::model::ising_tf_rf_nn::spin_dim;
//else if(model_type_str == "ising_selfdual_tf_rf_nn")
//spin_dim = settings::model::ising_selfdual_tf_rf_nn::spin_dim;
//else
//throw std::runtime_error(fmt::format("Unexpected model string: {}", model_type_str));
////
//    Eigen::Tensor<Scalar,3> A(spin_dimension,1,1);
//    Eigen::Tensor<Scalar,3> B(spin_dimension,1,1);
//    Eigen::Tensor<Scalar,1> LA(1);
//    Eigen::Tensor<Scalar,1> LC(1);
//    Eigen::Tensor<Scalar,1> LB(1);
//
//    // Default is a product state, spins pointing up in z.
//    A.setZero();
//    B.setZero();
//    LA.setConstant(1.0);
//    LB.setConstant(1.0);
//    LC.setConstant(1.0);
//
//    A(0,0,0) = 1;
//    B(0,0,0) = 1;
//    MPS_A = std::make_unique<class_mps_site>(class_mps_site(A,LA,0));
//    MPS_B = std::make_unique<class_mps_site>(class_mps_site(B,LB,1));
//}

class_mps_2site::class_mps_2site(const class_mps_2site &other)
    : swapped(other.swapped), MPS_A(copy_unique(other.MPS_A)), MPS_B(copy_unique(other.MPS_B))
{
}


//void class_mps_2site::initialize(int spin_dim){
//    spin_dimension = spin_dim;
//    Eigen::Tensor<Scalar,3> UA(spin_dimension,1,1);
//    Eigen::Tensor<Scalar,3> VB(spin_dimension,1,1);
//    Eigen::Tensor<Scalar,1> LA(1);
//    Eigen::Tensor<Scalar,1> LC(1);
//    Eigen::Tensor<Scalar,1> LB(1);
//
//    // Default is a product state, spins pointing up in z.
//    UA.setZero();
//    VB.setZero();
//    LA.setConstant(1.0);
//    LB.setConstant(1.0);
//    LC.setConstant(1.0);
//
//    UA(0,0,0) = 1;
//    VB(0,0,0) = 1;
//
//    set_mps(LA,UA,LC,VB,LB);
//    swapped = false;
//}


bool class_mps_2site::isReal()const {return MPS_A->isReal() and MPS_B->isReal();}
bool class_mps_2site::hasNaN()const {return MPS_A->hasNaN() or MPS_B->hasNaN();}
long class_mps_2site::chiA () const {return MPS_A->get_L().dimension(0);}
long class_mps_2site::chiB () const {return MPS_B->get_L().dimension(0);}
long class_mps_2site::chiC () const {return MPS_A->get_LC().dimension(0);}
long class_mps_2site::spindim() const {
    auto spindimA = MPS_A->get_M().dimension(0);
    auto spindimB = MPS_B->get_M().dimension(0);
    if(spindimA  != spindimB)
        throw std::runtime_error(fmt::format("Mismatch in spin dimension on A and B sites: {} != {}", spindimA,spindimB ));
    return spindimA;
}
const Eigen::Tensor<Scalar,3> & class_mps_2site::A_bare() const {return MPS_A->get_M_bare();}
const Eigen::Tensor<Scalar,3> & class_mps_2site::A() const {return MPS_A->get_M();}
const Eigen::Tensor<Scalar,3> & class_mps_2site::B() const {return MPS_B->get_M();}
Eigen::Tensor<Scalar,2> class_mps_2site::LC() const {return Textra::asDiagonal(MPS_A->get_LC());}

Eigen::Tensor<Scalar,3> class_mps_2site::GA() const {
    return MPS_A->get_M_bare()
    .contract(Textra::asDiagonalInversed(MPS_A->get_L()), Textra::idx({1},{1}))
    .shuffle(Textra::array3{0,2,1});
}
Eigen::Tensor<Scalar,3> class_mps_2site::GB() const {
    return MPS_B->get_M_bare()
    .contract(Textra::asDiagonalInversed(MPS_B->get_L()), Textra::idx({2},{0}));
}
Eigen::Tensor<Scalar,2> class_mps_2site::LA() const {return Textra::asDiagonal(MPS_A->get_L());}
Eigen::Tensor<Scalar,2> class_mps_2site::LB() const {return Textra::asDiagonal(MPS_B->get_L());}

void class_mps_2site::set_mps(const Eigen::Tensor<Scalar,1> &LA,
                              const Eigen::Tensor<Scalar,3> &A,
                              const Eigen::Tensor<Scalar,1> &LC,
                              const Eigen::Tensor<Scalar,3> &B,
                              const Eigen::Tensor<Scalar,1> &LB)
{
    MPS_A->set_mps(A, LA);
    MPS_A->set_LC(LC);
    MPS_B->set_mps(B, LB);
}

Eigen::DSizes<long,4> class_mps_2site::dimensions() const {return Eigen::DSizes<long,4>{spindim(),chiA(), spindim(),chiB()};}





void class_mps_2site::swap_AB() {
    tools::log->trace("Swapping AB");
    swapped = !swapped;
    // Store the positions
    auto position_left  = MPS_A->get_position();
    auto position_right = MPS_B->get_position();

    //Swap Gamma
    Eigen::Tensor<Scalar,1> LC = MPS_A->get_LC();
    MPS_A->unset_LC();
    MPS_B->unset_LC();
    MPS_A.swap(MPS_B);
    MPS_A->set_LC(MPS_A->get_L());
    MPS_A->set_L(LC);
    MPS_B->set_L(LC);
    MPS_A->set_position(position_left);
    MPS_B->set_position(position_right);

}


Eigen::Tensor<class_mps_2site::Scalar,4> class_mps_2site::get_theta(Scalar norm) const
/*!
 * Returns a two-site MPS
     @verbatim
        1--[ LB ]--[ GA ]--[ LA ]-- [ GB ] -- [ LB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return
            A().contract(B(), Textra::idx({2},{1})) / norm;
//    return
//            A().contract(LC(), idx({2},{0}))
//               .contract(B(), idx({2},{1})) / norm;
}
