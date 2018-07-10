//
// Created by david on 2017-11-13.
//


#include <mps_routines/class_mps_2site.h>
#include <general/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <general/class_arpack_eigsolver.h>
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
    theta = get_theta();
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



void class_mps_2site::compute_mps_components(){
//    int chiA2 = (int)(chiA()*chiA());
    int chiB2 = (int)(chiB()*chiB());
    int chiC2 = (int)(chiC()*chiC());


    Eigen::Tensor<Scalar,2> theta_evn_transfer_mat   = get_transfer_matrix_theta_evn().reshape(array2{chiB2,chiB2});
    Eigen::Tensor<Scalar,2> theta_odd_transfer_mat   = get_transfer_matrix_theta_odd().reshape(array2{chiC2,chiC2});
    class_arpack_eigsolver<Scalar, Form::GENERAL> solver;

    int ncvC = std::min(16, chiC2);
    int ncvB = std::min(16, chiB2);
    solver.eig(theta_evn_transfer_mat.data(), chiB2, 1, ncvB, Ritz::LM, Side::R, true, true);
    auto[eigvec_R_evn,eigval_R_evn]  = solver.get_eig_vecs_vals();
    solver.eig(theta_evn_transfer_mat.data(), chiB2, 1, ncvB, Ritz::LM, Side::L, true, true);
    auto[eigvec_L_evn,eigval_L_evn]  = solver.get_eig_vecs_vals();
    solver.eig(theta_odd_transfer_mat.data(), chiC2, 1, ncvC, Ritz::LM, Side::R, true, true);
    auto[eigvec_R_odd,eigval_R_odd]  = solver.get_eig_vecs_vals();
    solver.eig(theta_odd_transfer_mat.data(), chiC2, 1, ncvC, Ritz::LM, Side::L, true, true);
    auto[eigvec_L_odd,eigval_L_odd]  = solver.get_eig_vecs_vals();

    MatrixType<Scalar> eigvec_R_evn_map = Eigen::Map<const MatrixType<Scalar>>(eigvec_R_evn.data(), eigvec_R_evn.size(),1);
    MatrixType<Scalar> eigvec_L_evn_map = Eigen::Map<const MatrixType<Scalar>>(eigvec_L_evn.data(), eigvec_L_evn.size(),1);
    MatrixType<Scalar> eigvec_R_odd_map = Eigen::Map<const MatrixType<Scalar>>(eigvec_R_odd.data(), eigvec_R_odd.size(),1);
    MatrixType<Scalar> eigvec_L_odd_map = Eigen::Map<const MatrixType<Scalar>>(eigvec_L_odd.data(), eigvec_L_odd.size(),1);

    Scalar normalization_evn = sqrt((eigvec_L_evn_map.transpose() * eigvec_R_evn_map).sum());
    Scalar normalization_odd = sqrt((eigvec_L_odd_map.transpose() * eigvec_R_odd_map).sum());

    r_evn = Matrix_to_Tensor2(eigvec_R_evn_map).reshape(array2{chiB(),chiB()})/normalization_evn;
    l_evn = Matrix_to_Tensor2(eigvec_L_evn_map).reshape(array2{chiB(),chiB()})/normalization_evn;
    r_odd = Matrix_to_Tensor2(eigvec_R_odd_map).reshape(array2{chiC(),chiC()})/normalization_odd;
    l_odd = Matrix_to_Tensor2(eigvec_L_odd_map).reshape(array2{chiC(),chiC()})/normalization_odd;

    theta                = get_theta();
    theta_sw             = get_theta_swapped();
    theta_evn_normalized = get_theta_evn(sqrt(eigval_R_evn[0]));
    theta_odd_normalized = get_theta_odd(sqrt(eigval_R_odd[0]));

    LBGA                 = A();// / (T) sqrt(eigval_R_LBGA(0));
    LAGB                 = asDiagonal(LC).contract(MPS_B->get_G(), idx({1},{1})).shuffle(array3{1,0,2});// / (T) sqrt(eigval_R_LAGB(0));

    transfer_matrix_evn    = theta_evn_normalized.contract(theta_evn_normalized.conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3});
    transfer_matrix_odd    = theta_odd_normalized.contract(theta_odd_normalized.conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3});
    Eigen::Tensor<Scalar,4> transfer_matrix_LAGB_unnormalized = LAGB.contract(LAGB.conjugate(), idx({0},{0})).shuffle(array4{0,2,1,3});
    Eigen::Tensor<Scalar,4> transfer_matrix_LBGA_unnormalized = LBGA.contract(LBGA.conjugate(), idx({0},{0})).shuffle(array4{0,2,1,3});

    Eigen::Tensor<Scalar,0> l_evn_LBGA_r_odd = l_evn.contract(transfer_matrix_LBGA_unnormalized, idx({0,1},{0,1})).contract(r_odd, idx({0,1},{0,1}));
    Eigen::Tensor<Scalar,0> l_odd_LAGB_r_evn = l_odd.contract(transfer_matrix_LAGB_unnormalized, idx({0,1},{0,1})).contract(r_evn, idx({0,1},{0,1}));

    transfer_matrix_LAGB = transfer_matrix_LAGB_unnormalized /l_odd_LAGB_r_evn(0);
    transfer_matrix_LBGA = transfer_matrix_LBGA_unnormalized /l_evn_LBGA_r_odd(0);
    LAGB = LAGB / sqrt(l_odd_LAGB_r_evn(0));
    LBGA = LBGA / sqrt(l_evn_LBGA_r_odd(0));

//    std::cout << "Check:" << setprecision(10) <<  std::endl;
//    std::cout << " l_odd_LAGB_r_evn          = " << l_odd_LAGB_r_evn(0) << std::endl;
//    std::cout << " l_evn_LBGA_r_odd          = " << l_evn_LBGA_r_odd(0) << std::endl;
//    std::cout << " < l_evn | r_evn >         = " << l_evn.contract(r_evn, idx({0,1},{0,1})) << std::endl;
//    std::cout << " < l_odd | r_odd >         = " << l_odd.contract(r_odd, idx({0,1},{0,1})) << std::endl;
//    std::cout << " < l_evn | LBGA  | r_odd > = " << l_evn.contract(transfer_matrix_LBGA, idx({0,1},{0,1})).contract(r_odd, idx({0,1},{0,1})) << std::endl;
//    std::cout << " < l_odd | LAGB  | r_evn > = " << l_odd.contract(transfer_matrix_LAGB, idx({0,1},{0,1})).contract(r_evn, idx({0,1},{0,1})) << std::endl;
//    std::cout << " < theta     | theta >     = " << theta.contract(theta.conjugate(), idx({1,3,0,2},{1,3,0,2})) << std::endl;
//    std::cout << " < theta_evn_normalized | theta_evn_normalized > = " << theta_evn_normalized.contract(theta_evn_normalized.conjugate(), idx({0,2},{0,2})).contract(l_evn, idx({0,2},{0,1})).contract(r_evn,idx({0,1},{0,1})) << std::endl;
//    std::cout << " < theta_evn_normalized | theta_evn_normalized > = " << transfer_matrix_evn.contract(l_evn, idx({0,1},{0,1})).contract(r_evn,idx({0,1},{0,1})) << std::endl;
//    std::cout << " < theta_odd_normalized | theta_odd_normalized > = " << theta_odd_normalized.contract(theta_odd_normalized.conjugate(), idx({0,2},{0,2})).contract(l_odd, idx({0,2},{0,1})).contract(r_odd,idx({0,1},{0,1})) << std::endl;
//    std::cout << " < theta_odd_normalized | theta_odd_normalized > = " << transfer_matrix_odd.contract(l_odd, idx({0,1},{0,1})).contract(r_odd,idx({0,1},{0,1})) << std::endl;
}




Eigen::Tensor<Scalar,4> class_mps_2site::get_theta(Scalar norm) const
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
            A().contract(asDiagonal(LC), idx({2},{0}))
               .contract(B(), idx({2},{1})) / norm;
}



Eigen::Tensor<Scalar,4> class_mps_2site::get_theta_swapped(Scalar norm) const
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LA ]--[ GB ]--[ LB ]-- [ GA ] -- [ LA ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return  asDiagonal(LC) //whatever L_A was in the previous step
            .contract(B(),            idx({1},{1}))
            .contract(MPS_A->get_G(), idx({2},{1}))
            .contract(asDiagonal(LC), idx({3},{0}))
            .shuffle(array4{1,0,2,3})
            /norm;
}





Eigen::Tensor<Scalar,4> class_mps_2site::get_theta_evn(Scalar norm) const
/*!
 * Returns a right normalized two-site MPS
     @verbatim
        1--[ LB ]--[ GA ]-- [ LC ] -- [ GB ]--3
                     |                 |
                     0                 2
     @endverbatim
 */

{
    return  A()
            .contract(asDiagonal(LC),  idx({2},{0}))
            .contract(MPS_B->get_G(),  idx({2},{1}))
//            .shuffle(array4{1,0,2,3})
            /norm;
};

Eigen::Tensor<Scalar,4> class_mps_2site::get_theta_odd(Scalar norm) const
/*!
 * Returns a two-site MPS with A and B swapped
     @verbatim
        1--[ LA ]--[ GB ]-- [ LB ] -- [ GA ]--3
                     |                 |
                     0                 2
     @endverbatim
 */
{
    return  asDiagonal(LC)
            .contract(MPS_B->get_G(),         idx({1},{1}))
            .contract(A(),                    idx({2},{1}))
            .shuffle(array4{1,0,2,3})
            /norm;
}


Eigen::Tensor<Scalar,4> class_mps_2site::get_transfer_matrix_zero() const {
    Eigen::Tensor<Scalar,1> I = LC;
    I.setConstant(1.0);
    Eigen::array<Eigen::IndexPair<long>,0> pair = {};

    return asDiagonal(I).contract(asDiagonal(I), pair ).shuffle(array4{0,2,1,3});
};



Eigen::Tensor<Scalar,4> class_mps_2site::get_transfer_matrix_LBGA(Scalar norm) const {
    return A().contract( A().conjugate() , idx({0},{0}))
              .shuffle(array4{0,3,1,2})
            /norm;
}


Eigen::Tensor<Scalar,4> class_mps_2site::get_transfer_matrix_GALC(Scalar norm) const {
    return asDiagonal(LC)
            .contract(MPS_A->get_G(),               idx({2},{0}))
            .contract(MPS_A->get_G().conjugate(),   idx({0},{0}))
            .contract(asDiagonal(LC),               idx({3},{0}) )
            .shuffle(array4{0,2,1,3})
           /norm;
}

Eigen::Tensor<Scalar,4> class_mps_2site::get_transfer_matrix_GBLB(Scalar norm) const {
    return B().contract(B().conjugate() ,   idx({0},{0}))
              .shuffle(array4{0,2,1,3})
              /norm;
}


Eigen::Tensor<Scalar,4> class_mps_2site::get_transfer_matrix_LCGB(Scalar norm) const {
    return  asDiagonal(LC)
            .contract(MPS_B->get_G(),               idx({1},{1}))
            .contract(MPS_B->get_G().conjugate(),   idx({1},{0}))
            .contract(asDiagonal(LC),               idx({2},{1}) )
            .shuffle(array4{0,3,1,2})
            /norm;
}


Eigen::Tensor<Scalar,4> class_mps_2site::get_transfer_matrix_theta_evn(Scalar norm) const {
    return get_theta_evn().contract(get_theta_evn().conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3}) / norm;
}

Eigen::Tensor<Scalar,4> class_mps_2site::get_transfer_matrix_theta_odd(Scalar norm) const {
    return get_theta_odd().contract(get_theta_odd().conjugate(), idx({0,2},{0,2})).shuffle(array4{0,2,1,3}) / norm;
}


Eigen::Tensor<Scalar,4> class_mps_2site::get_transfer_matrix_AB(int p) const {
    Eigen::Tensor<Scalar,4> temp = get_transfer_matrix_zero();
    Eigen::Tensor<Scalar,4> temp2;
    for (int i = 0; i < p-2; i++){
        if(Math::mod(i,2) == 0){
            temp2 = temp.contract(get_transfer_matrix_LBGA(), idx({2,3},{0,1}));

        }else{
            temp2 = temp.contract(get_transfer_matrix_LCGB(), idx({2,3},{0,1}));
        }
        temp = temp2;


    }
    return temp;
}
