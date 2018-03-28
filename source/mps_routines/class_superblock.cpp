//
// Created by david on 7/22/17.
//

#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_mps.h>
#include <mps_routines/class_mpo.h>
#include <iomanip>
#include <SymEigsSolver.h>
#include <MatOp/SparseSymMatProd.h>
#include <general/class_svd_wrapper.h>
#include <mps_routines/class_custom_contraction.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include "arscomp.h"


using namespace std;
using namespace Textra;
using Scalar = class_superblock::Scalar;

class_superblock::class_superblock():
        MPS(std::make_shared<class_mps>()),
        H(std::make_shared<class_mpo>()),
        Lblock(std::make_shared<class_environment>("L")),
        Rblock(std::make_shared<class_environment>("R")),
        Lblock2(std::make_shared<class_environment_var>("L")),
        Rblock2(std::make_shared<class_environment_var>("R")),
        SVD(std::make_shared<class_SVD<Scalar>>())
{
    t_eig.set_properties(true, 5," Time: ");
    MPS->initialize(H->local_dimension);
    M_minus_E  = H->compute_M();
//    MM_minus_E = H->compute_MM();
    chain_length = H->mps_sites;
    set_current_dimensions();
}


//============================================================================//
// Find smallest eigenvalue using Arpack.
//============================================================================//

Textra::Tensor<Scalar,4> class_superblock::optimize_MPS(Textra::Tensor<Scalar, 4> &theta, int eigSteps, double eigThreshold){
    class_custom_contraction<Scalar>  esp(Lblock->block, Rblock->block, H->M, shape4);
    int ncv = std::min(settings::precision::eig_max_ncv,(int)theta.size());
    int nev = 1;
    int dim = esp.cols();
    t_eig.tic();
    //The operator A in "Ax = Ex" should be real and symmetric, so the eigenvector should be real as well, but might come out with a complex phase factor that can be removed.
    ARCompStdEig<double, class_custom_contraction<Scalar>> eig (dim, nev, &esp, &class_custom_contraction<Scalar>::MultMv, "SR", ncv,eigThreshold,eigSteps, theta.data());
    eig.ChangeTol(1e-12);
    eig.FindEigenvectors();
    int rows = eig.GetN();
    int cols = std::min(1, eig.GetNev());
    Textra::Tensor<Scalar, 2> eigvecs = TensorMap<Tensor<Scalar,2>> (eig.RawEigenvectors(), rows,cols);
    Textra::Tensor<Scalar, 1> eigvals = TensorMap<Tensor<Scalar,1>> (eig.RawEigenvalues(), cols);
    Scalar inv_phase = -1.0i * std::arg(eigvecs(0,0));
    eigvecs = (eigvecs *  std::exp(inv_phase)); //Rotate the phase of the eigenvector back so that it becomes real.
    eigvecs =  eigvecs.real().cast<Scalar>();    //Remove the, now negligible, imaginary part
    t_eig.toc();

    using namespace chrono;
    energy = std::real(eigvals(0));

    std::cout << "Time: " << duration_cast<duration<double>>(t_eig.measured_time).count() << " ";
    t_eig.print_delta();
    std::cout << " iter: " << eig.GetIter() << " energy: " << energy << "\n";
    return eigvecs.reshape(shape4);
}





//============================================================================//
// Do unitary evolution on an MPS
//============================================================================//
Textra::Tensor<Scalar, 4> class_superblock::evolve_MPS(const Textra::Tensor<Scalar, 4> &U)
/*!
@verbatim
  1--[ Θ ]--3
     |   |
     0   2
                   1--[ Θ ]--3
     0   1   --->     |   |
     |   |            0   2
     [ U ]
     |   |
     2   3
@endverbatim
*/

{
    return U.contract(MPS->get_theta(), idx({0,1},{0,2}))
            .shuffle(array4{0,2,1,3});
}

Textra::Tensor<Scalar, 4> class_superblock::evolve_MPS(const Textra::Tensor<Scalar, 4> &theta, const Textra::Tensor<Scalar, 4> &U)
/*!
@verbatim
  1--[ Θ ]--3
     |   |
     0   2
                   1--[ Θ ]--3
     0   1   --->     |   |
     |   |            0   2
     [ U ]
     |   |
     2   3
@endverbatim
*/
{
    return U.contract(theta, idx({0,1},{0,2}))
            .shuffle(array4{0,2,1,3});
}

//============================================================================//
// Do SVD decomposition, truncation and normalization of the MPS->
//============================================================================//
Textra::Tensor<Scalar,4> class_superblock::truncate_MPS(const Textra::Tensor<Scalar, 4> &theta,long chi_max_, double SVDThreshold){
    SVD->setThreshold(SVDThreshold);
    auto[U, S, V] = SVD->schmidt(theta, d, chiL, chi_max_, chiR);
    MPS->truncation_error = SVD->get_truncation_error();
    MPS->LA  = S;
    MPS->GA  = asDiagonalInversed(MPS->L_tail).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
    MPS->GB  = V.contract(asDiagonalInversed(MPS->LB), idx({2},{0}));
    return MPS->get_theta();
}

void class_superblock::truncate_MPS(const Textra::Tensor<Scalar, 4> &theta, const std::shared_ptr<class_mps> &MPS_given,long chi_max_, double SVDThreshold){
    SVD->setThreshold(SVDThreshold);
    auto[U, S, V] = SVD->schmidt(theta, d, chiL, chi_max_, chiR);
    MPS_given->truncation_error = SVD->get_truncation_error();
    MPS_given->GA  = asDiagonalInversed(MPS_given->L_tail).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
    MPS_given->GB  = V.contract(asDiagonalInversed(MPS_given->LB), idx({2},{0}));
}



void class_superblock::enlarge_environment(int direction){
    Tensor<Scalar, 0>  E_two_site =
            Lblock->block
            .contract(MPS->theta,               idx({0},{1}))
            .contract(H->M,                     idx({1,2},{0,2}))
            .contract(H->M,                     idx({3,1},{0,2}))
            .contract(MPS->theta.conjugate(),   idx({0,2,4},{1,0,2}))
            .contract(Rblock->block,            idx({0,2,1},{0,1,2}));

//    Tensor<Scalar, 0>  Etest2 =
//            Lblock->block
//                    .contract(MPS->theta,               idx({0},{1}))
//                    .contract(M_minus_E,                     idx({1,2},{0,2}))
//                    .contract(M_minus_E,                     idx({3,1},{0,2}))
//                    .contract(MPS->theta.conjugate(),   idx({0,2,4},{1,0,2}))
//                    .contract(Rblock->block,            idx({0,2,1},{0,1,2}));
//    Tensor<Scalar, 0>  Etest3 =
//            Lblock->block
//                    .contract(asDiagonal(MPS->LA), idx({0},{0}))
//                    .contract(asDiagonal(MPS->LA), idx({0},{0}))
//                    .contract(Rblock->block,       idx({1,2,0},{0,1,2}));

    double E_one_site = E_two_site(0).real()/2.0;



    M_minus_E    = H->compute_M   (E_one_site);
    if (direction == 1){
        Lblock->enlarge(MPS, M_minus_E);
        Lblock2->enlarge(MPS, M_minus_E);
    }else if (direction == -1){
        Rblock->enlarge(MPS, M_minus_E);
        Rblock2->enlarge(MPS, M_minus_E);
    }else if(direction == 0){
        Lblock->enlarge(MPS,  M_minus_E);
        Rblock->enlarge(MPS,  M_minus_E);
        Lblock2->enlarge(MPS, M_minus_E);
        Rblock2->enlarge(MPS, M_minus_E);
        chain_length += 2;
//
//        Tensor<Scalar, 0>  Etest =
//                Lblock->block
//                        .contract(asDiagonal(MPS->LA), idx({0},{0}))
//                        .contract(asDiagonal(MPS->LA), idx({0},{0}))
//                        .contract(Rblock->block,       idx({1,2,0},{0,1,2}));
//
//        std::cout << Etest << std::endl;
//        std::cout << Etest2 << std::endl;
//        std::cout << Etest3 << std::endl;
    }
}


void class_superblock::set_current_dimensions(){
    d       = H->local_dimension;
    chi     = MPS->GA.dimension(2);
    chiL    = MPS->GA.dimension(1);
    chiR    = MPS->GB.dimension(2);
    shape1  = {d * chiL * d * chiR};
    shape2  = {d * chiL , d * chiR};
    shape4  = {d , chiL , d , chiR};
}


void  class_superblock::swap_AB(){
    MPS->swap_AB();
}

