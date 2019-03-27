//
// Created by david on 7/22/17.
//

//#include <mps_routines/class_optimize_mps.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_mps_2site.h>
#include <mps_routines/nmspc_mps_tools.h>
#include <iomanip>
#include <general/class_svd_wrapper.h>
#include <general/class_eigsolver.h>
#include <general/arpack_extra/matrix_product_hamiltonian.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <model/class_hamiltonian_factory.h>
#define profile_optimization 0

using namespace std;
using namespace Textra;
using Scalar = class_superblock::Scalar;

class_superblock::class_superblock(SimulationType sim_type_):
        sim_type(sim_type_),
        MPS(std::make_shared<class_mps_2site>()),
        HA (class_hamiltonian_factory::create_mpo(settings::model::model_type)),
        HB (class_hamiltonian_factory::create_mpo(settings::model::model_type)),
        Lblock(std::make_shared<class_environment>("L")),
        Rblock(std::make_shared<class_environment>("R")),
        Lblock2(std::make_shared<class_environment_var>("L")),
        Rblock2(std::make_shared<class_environment_var>("R")),
        SVD(std::make_shared<class_SVD<Scalar>>())
{
    HA->set_position(0);
    HB->set_position(1);
    t_eig.set_properties(profile_optimization, 10,"Time optimizing ");
    spin_dimension = HA->get_spin_dimension();
    MPS->initialize(spin_dimension);
    set_profiling_labels();
}




// We need to make a destructor manually for the enclosing class "class_superblock"
// that encloses "class_hamiltonian_base". Otherwise unique_ptr will forcibly inline its
// own default deleter.
// This allows us to forward declare the abstract base class "class_hamiltonian_base"
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
//class_superblock::~class_superblock()=default;




void class_superblock::clear(){
    *this = class_superblock(sim_type);
}


size_t class_superblock::get_length() const {
    size_t length = Lblock->size + Rblock->size + 2;
    assert(length == environment_size + 2 and "ERROR: System length mismatch!");
    return length;

}

size_t class_superblock::get_position() const { return Lblock->size;}

size_t class_superblock::get_chi() const {return MPS->chiC();}

Eigen::DSizes<long,4> class_superblock::dimensions()const{return MPS->dimensions();}

Eigen::Tensor<Scalar, 4> class_superblock::get_theta() const { return MPS->get_theta();}



//============================================================================//
// Find smallest eigenvalue using Arpack.
//============================================================================//



Eigen::Tensor<Scalar,4> class_superblock::optimize_MPS(Eigen::Tensor<Scalar, 4> &theta, eigutils::eigSetting::Ritz ritz){
    std::array<long,4> shape_theta4 = theta.dimensions();
    std::array<long,4> shape_mpo4   = HA->MPO.dimensions();

    t_eig.tic();
    int nev = 1;
    using namespace settings::precision;
    using namespace eigutils::eigSetting;
    DenseHamiltonianProduct<Scalar>  matrix (Lblock->block.data(), Rblock->block.data(), HA->MPO.data(), HB->MPO.data(), shape_theta4, shape_mpo4);
    class_eigsolver solver;
    solver.eigs_dense(matrix,nev,eigMaxNcv,NAN,Form::SYMMETRIC,ritz,Side::R,true,true);

    auto eigvals           = Eigen::TensorMap<const Eigen::Tensor<double,1>>  (solver.solution.get_eigvals<Form::SYMMETRIC>().data() ,solver.solution.meta.cols);
    auto eigvecs           = Eigen::TensorMap<const Eigen::Tensor<Scalar,1>>  (solver.solution.get_eigvecs<Type::CPLX, Form::SYMMETRIC>().data(),solver.solution.meta.rows);

    t_eig.toc();
    t_eig.print_delta();

    E_optimal = std::real(eigvals(0));
    return eigvecs.reshape(theta.dimensions());
}




//============================================================================//
// Do unitary evolution on an MPS
//============================================================================//
Eigen::Tensor<Scalar, 4> class_superblock::evolve_MPS(const Eigen::Tensor<Scalar, 4> &U)
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

Eigen::Tensor<Scalar, 4> class_superblock::evolve_MPS(const Eigen::Tensor<Scalar, 4> &theta, const Eigen::Tensor<Scalar, 4> &U)
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
Eigen::Tensor<Scalar,4> class_superblock::truncate_MPS(const Eigen::Tensor<Scalar, 4> &theta,long chi_, double SVDThreshold){
    SVD->setThreshold(SVDThreshold);
    auto[U, S, V] = SVD->schmidt(theta,chi_);
    MPS->truncation_error         = SVD->get_truncation_error();
    measurements.truncation_error = SVD->get_truncation_error();
    MPS->LC  = S;
    Eigen::Tensor<Scalar,3> L_U = asDiagonalInversed(MPS->MPS_A->get_L()).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
    Eigen::Tensor<Scalar,3> V_L = V.contract(asDiagonalInversed(MPS->MPS_B->get_L()), idx({2},{0}));
    MPS->MPS_A->set_G(L_U);
    MPS->MPS_B->set_G(V_L);
    return get_theta();
}

void class_superblock::truncate_MPS(const Eigen::Tensor<Scalar, 4> &theta, const std::shared_ptr<class_mps_2site> &MPS_out,long chi_, double SVDThreshold){
    SVD->setThreshold(SVDThreshold);
    auto[U, S, V] = SVD->schmidt(theta, chi_);
    MPS_out->truncation_error = SVD->get_truncation_error();
    MPS_out->LC  = S;
    Eigen::Tensor<Scalar,3> L_U = asDiagonalInversed(MPS_out->MPS_A->get_L()).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
    Eigen::Tensor<Scalar,3> V_L = V.contract(asDiagonalInversed(MPS_out->MPS_B->get_L()), idx({2},{0}));
    MPS_out->MPS_A->set_G(L_U);
    MPS_out->MPS_B->set_G(V_L);

}





Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic> class_superblock::get_H_local_matrix () const{
    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
    Eigen::Tensor<Scalar,5>tempL   = Lblock->block.contract(HA->MPO,Textra::idx({2},{0})).shuffle(Textra::array5{4,1,3,0,2});
    Eigen::Tensor<Scalar,5>tempR   = Rblock->block.contract(HB->MPO,Textra::idx({2},{1})).shuffle(Textra::array5{4,1,3,0,2});
    Eigen::Tensor<Scalar,8>H_local = tempL.contract(tempR, Textra::idx({4},{4})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
//    tempL.resize(0,0,0,0,0);
//    tempR.resize(0,0,0,0,0);
//    tempH = tempH.shuffle(Textra::array8{0,1,4,5,2,3,6,7}).eval();
    return Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(H_local.data(),shape,shape);
}

template<typename Scalar, auto rank>
bool has_imaginary_part(const Eigen::Tensor<Scalar,rank> &tensor, double threshold = 1e-14) {
    Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>> vector (tensor.data(),tensor.size());
    if constexpr (std::is_same<Scalar, std::complex<double>>::value){
        return vector.imag().cwiseAbs().sum() > threshold;
    }else{
        return false;
    }
}
Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> class_superblock::get_H_local_matrix_real ()const{
    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
    if(has_imaginary_part(Lblock->block)){throw std::runtime_error("Discarding imaginary data from Lblock when building H_local");}
    if(has_imaginary_part(Rblock->block)){throw std::runtime_error("Discarding imaginary data from Rblock when building H_local");}
    if(has_imaginary_part(HA->MPO))      {throw std::runtime_error("Discarding imaginary data from MPO A when building H_local");}
    if(has_imaginary_part(HB->MPO))      {throw std::runtime_error("Discarding imaginary data from MPO B when building H_local");}

    Eigen::Tensor<double,5>tempL   = Lblock->block.contract(HA->MPO,Textra::idx({2},{0})).real().shuffle(Textra::array5{4,1,3,0,2});
    Eigen::Tensor<double,5>tempR   = Rblock->block.contract(HB->MPO,Textra::idx({2},{1})).real().shuffle(Textra::array5{4,1,3,0,2});
    Eigen::Tensor<double,8>H_local = tempL.contract(tempR, Textra::idx({4},{4})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
    return Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>>(H_local.data(),shape,shape);
}


Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic> class_superblock::get_H_local_sq_matrix ()const{
    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
    Eigen::Tensor<Scalar,6>tempL = Lblock2->block.contract(HA->MPO,Textra::idx({2},{0}))
                                                 .contract(HA->MPO,Textra::idx({2,5},{0,2}))
                                                 .shuffle(Textra::array6{5,1,3,0,2,4});
    Eigen::Tensor<Scalar,6>tempR = Rblock2->block.contract(HB->MPO,Textra::idx({2},{1}))
                                                 .contract(HB->MPO,Textra::idx({2,5},{1,2}))
                                                 .shuffle(Textra::array6{5,1,3,0,2,4});

    Eigen::Tensor<Scalar,8> H_local = tempL.contract(tempR, Textra::idx({4,5},{4,5})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
    return Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(H_local.data(),shape,shape);
}

Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> class_superblock::get_H_local_sq_matrix_real ()const{
    if(has_imaginary_part(Lblock2->block)){throw std::runtime_error("Discarding imaginary data from Lblock when building H_local");}
    if(has_imaginary_part(Rblock2->block)){throw std::runtime_error("Discarding imaginary data from Rblock when building H_local");}
    if(has_imaginary_part(HA->MPO))      {throw std::runtime_error("Discarding imaginary data from MPO A when building H_local");}
    if(has_imaginary_part(HB->MPO))      {throw std::runtime_error("Discarding imaginary data from MPO B when building H_local");}
    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
    Eigen::Tensor<double,6>tempL = Lblock2->block.contract(HA->MPO,Textra::idx({2},{0}))
                                                 .contract(HA->MPO,Textra::idx({2,5},{0,2}))
                                                 .real().shuffle(Textra::array6{5,1,3,0,2,4});
    Eigen::Tensor<double,6>tempR = Rblock2->block.contract(HB->MPO,Textra::idx({2},{1}))
                                                 .contract(HB->MPO,Textra::idx({2,5},{1,2}))
                                                 .real().shuffle(Textra::array6{5,1,3,0,2,4});

    Eigen::Tensor<double,8>H_local = tempL.contract(tempR, Textra::idx({4,5},{4,5})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
    return Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>>(H_local.data(),shape,shape);
}


Textra::SparseMatrixType<Scalar> class_superblock::get_H_local_sparse_matrix (double prune)const{
    return Textra::Tensor2_to_SparseMatrix(get_H_local_rank2(),prune);
}
Textra::SparseMatrixType<Scalar> class_superblock::get_H_local_sq_sparse_matrix (double prune)const{
    return Textra::Tensor2_to_SparseMatrix(get_H_local_sq_rank2(),prune);
}


Eigen::Tensor<Scalar,2> class_superblock::get_H_local_rank2 ()const{
    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
    Eigen::Tensor<Scalar,6>tempL = Lblock2->block.contract(HA->MPO,Textra::idx({2},{0}))
            .contract(HA->MPO,Textra::idx({2,5},{0,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});
    Eigen::Tensor<Scalar,6>tempR = Rblock2->block.contract(HB->MPO,Textra::idx({2},{1}))
            .contract(HB->MPO,Textra::idx({2,5},{1,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});

    return tempL.contract(tempR, Textra::idx({4,5},{4,5})).shuffle(Textra::array8{0,1,4,5,2,3,6,7}).reshape(Textra::array2{shape, shape});
}

Eigen::Tensor<Scalar,8> class_superblock::get_H_local_rank8 ()const{
    Eigen::Tensor<Scalar,6>tempL = Lblock2->block.contract(HA->MPO,Textra::idx({2},{0}))
            .contract(HA->MPO,Textra::idx({2,5},{0,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});
    Eigen::Tensor<Scalar,6>tempR = Rblock2->block.contract(HB->MPO,Textra::idx({2},{1}))
            .contract(HB->MPO,Textra::idx({2,5},{1,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});

    return tempL.contract(tempR, Textra::idx({4,5},{4,5})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
}

Eigen::Tensor<Scalar,2> class_superblock::get_H_local_sq_rank2 ()const{

    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
    Eigen::Tensor<Scalar,6>tempL = Lblock2->block.contract(HA->MPO,Textra::idx({2},{0}))
            .contract(HA->MPO,Textra::idx({2,5},{0,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});
    Eigen::Tensor<Scalar,6>tempR = Rblock2->block.contract(HB->MPO,Textra::idx({2},{1}))
            .contract(HB->MPO,Textra::idx({2,5},{1,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});

    return tempL.contract(tempR, Textra::idx({4,5},{4,5})).shuffle(Textra::array8{0,1,4,5,2,3,6,7}).reshape(Textra::array2{shape, shape});
}

Eigen::Tensor<Scalar,8> class_superblock::get_H_local_sq_rank8 ()const{
    Eigen::Tensor<Scalar,6>tempL = Lblock2->block.contract(HA->MPO,Textra::idx({2},{0}))
            .contract(HA->MPO,Textra::idx({2,5},{0,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});
    Eigen::Tensor<Scalar,6>tempR = Rblock2->block.contract(HB->MPO,Textra::idx({2},{1}))
            .contract(HB->MPO,Textra::idx({2,5},{1,2}))
            .shuffle(Textra::array6{5,1,3,0,2,4});

    return tempL.contract(tempR, Textra::idx({4,5},{4,5})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
}





void class_superblock::enlarge_environment(int direction){
    if (direction == 1){
        assert(Lblock->get_position()  == HA->get_position());
        assert(Lblock2->get_position() == HA->get_position());
        Lblock->enlarge(*MPS,  HA->MPO_reduced_view());
        Lblock2->enlarge(*MPS, HA->MPO_reduced_view());
        Lblock->set_position (HB->get_position());
        Lblock2->set_position(HB->get_position());
    }else if (direction == -1){
        assert(Rblock->get_position()  == HB->get_position());
        assert(Rblock2->get_position() == HB->get_position());
        Rblock->enlarge(*MPS,  HB->MPO_reduced_view());
        Rblock2->enlarge(*MPS, HB->MPO_reduced_view());
        Rblock->set_position (HA->get_position());
        Rblock2->set_position(HA->get_position());
    }else if(direction == 0){
        Lblock->enlarge(*MPS,  HA->MPO_reduced_view());
        Rblock->enlarge(*MPS,  HB->MPO_reduced_view());
        Lblock2->enlarge(*MPS, HA->MPO_reduced_view());
        Rblock2->enlarge(*MPS, HB->MPO_reduced_view());

        Lblock->set_position (HB->get_position());
        Rblock->set_position (HB->get_position()+1);
        Lblock2->set_position(HB->get_position());
        Rblock2->set_position(HB->get_position()+1);
        environment_size = Lblock->size + Rblock->size;
    }
}


void class_superblock::set_positions(int position){
    MPS->MPS_A->set_position(position);
    MPS->MPS_B->set_position(position+1);
    Lblock->set_position(position);
    Rblock->set_position(position+1);
    Lblock2->set_position(position);
    Rblock2->set_position(position+1);
    HA->set_position(position);
    HB->set_position(position+1);
}


void class_superblock::set_measured_false(){
    has_been_measured = false;
    has_been_written  = false;
    MPS_Tools::Common::Views::components_computed = false;
}

void class_superblock::do_all_measurements() {
    using namespace MPS_Tools::Common;
    if (has_been_measured){return;}
    measurements.length                         = Measure::length(*this);
    measurements.bond_dimension                 = Measure::bond_dimension(*this);
    measurements.norm                           = Measure::norm(*this);
//    measurements.truncation_error               = Measure::truncation_error(*this);
    measurements.energy_mpo                     = Measure::energy_mpo(*this);  //This number is needed for variance calculation!
    measurements.energy_per_site_mpo            = Measure::energy_per_site_mpo(*this);
    measurements.energy_per_site_ham            = Measure::energy_per_site_ham(*this);
    measurements.energy_per_site_mom            = Measure::energy_per_site_mom(*this);
    measurements.energy_variance_mpo            = Measure::energy_variance_mpo(*this, measurements.energy_mpo);
    measurements.energy_variance_per_site_mpo   = Measure::energy_variance_per_site_mpo(*this, measurements.energy_mpo);
    measurements.energy_variance_per_site_ham   = Measure::energy_variance_per_site_ham(*this);
    measurements.energy_variance_per_site_mom   = Measure::energy_variance_per_site_mom(*this);
    measurements.current_entanglement_entropy   = Measure::current_entanglement_entropy(*this);
    has_been_measured = true;
}



void  class_superblock::swap_AB(){
    MPS->swap_AB();
}


// Profiling

void class_superblock::set_profiling_labels() {
    using namespace settings::profiling;
    t_ene_mpo.set_properties(on, precision,"↳ Energy (MPO)           ");
    t_ene_ham.set_properties(on, precision,"↳ Energy (HAM)           ");
    t_ene_mom.set_properties(on, precision,"↳ Energy (MOM)           ");
    t_var_mpo.set_properties(on, precision,"↳ Variance (MPO)         ");
    t_var_ham.set_properties(on, precision,"↳ Variance (HAM)         ");
    t_var_mom.set_properties(on, precision,"↳ Variance (MOM)         ");
    t_entropy.set_properties(on, precision,"↳ Ent. Entropy           ");
    t_temp1.set_properties(on, precision,  "↳ Temp1                  ");
    t_temp2.set_properties(on, precision,  "↳ Temp2                  ");
    t_temp3.set_properties(on, precision,  "↳ Temp3                  ");
    t_temp4.set_properties(on, precision,  "↳ Temp4                  ");

}

void class_superblock::print_profiling(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        std::cout << "\nComputing observables breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_ene_mpo.print_time_w_percent(t_parent);
        t_ene_ham.print_time_w_percent(t_parent);
        t_ene_mom.print_time_w_percent(t_parent);
        t_var_mpo.print_time_w_percent(t_parent);
        t_var_ham.print_time_w_percent(t_parent);
        t_var_mom.print_time_w_percent(t_parent);
        t_entropy.print_time_w_percent(t_parent);
        t_temp1.print_time_w_percent(t_parent);
        t_temp2.print_time_w_percent(t_parent);
        t_temp3.print_time_w_percent(t_parent);
        t_temp4.print_time_w_percent(t_parent);
    }
}
