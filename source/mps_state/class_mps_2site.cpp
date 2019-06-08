//
// Created by david on 2017-11-13.
//


#include <mps_routines/class_mps_2site.h>
#include <general/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <iomanip>

using namespace std;
using namespace Textra;

using Scalar = class_mps_2site::Scalar;

// Function definitions

class_mps_2site::class_mps_2site(){
    MPS_A = std::make_unique<class_vidal_mps>();
    MPS_B = std::make_unique<class_vidal_mps>();
}

class_mps_2site::class_mps_2site(const class_mps_2site &other)
    : spin_dimension(other.spin_dimension), swapped(other.swapped), MPS_A(copy_unique(other.MPS_A)), MPS_B(copy_unique(other.MPS_B)),  LC(other.LC)
{
}


void class_mps_2site::initialize(int spin_dim){
    spin_dimension = spin_dim;
    Eigen::Tensor<Scalar,3> GA(spin_dimension,1,1);
    Eigen::Tensor<Scalar,3> GB(spin_dimension,1,1);
    Eigen::Tensor<Scalar,1> LA(1);
    Eigen::Tensor<Scalar,1> LC(1);
    Eigen::Tensor<Scalar,1> LB(1);

    // Default is a product state, spins pointing up in z.
    GA.setZero();
    GB.setZero();
    LA.setConstant(1.0);
    LB.setConstant(1.0);
    LC.setConstant(1.0);

    GA(0,0,0) = 1;
    GB(0,0,0) = 1;

    set_mps(LA,GA,LC,GB,LB);
    swapped = false;
//    theta = get_theta();
}


void class_mps_2site::swap_AB() {
    swapped = !swapped;

    //Swap Gamma
    tmp3 = MPS_A->get_G();
    MPS_A->set_G(MPS_B->get_G());
    MPS_B->set_G(tmp3);


    tmp1 = LC;
    LC = MPS_B->get_L();
    MPS_A->set_L(tmp1);
    MPS_B->set_L(tmp1);
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
            A().contract(C(), idx({2},{0}))
               .contract(B(), idx({2},{1})) / norm;
}
