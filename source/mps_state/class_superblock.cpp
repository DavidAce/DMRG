//
// Created by david on 7/22/17.
//

//#include <mps_state/class_optimize_mps.h>
#include <mps_state/class_superblock.h>
#include <mps_state/class_environment.h>
#include <mps_state/class_mps_2site.h>
#include <mps_tools/nmspc_mps_tools.h>
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

class_superblock::class_superblock(SimulationType sim_type_,std::string sim_name_):
        sim_type(sim_type_),
        sim_name(sim_name_),
        MPS(std::make_shared<class_mps_2site>()),
        HA (class_hamiltonian_factory::create_mpo(0,settings::model::model_type)),
        HB (class_hamiltonian_factory::create_mpo(1,settings::model::model_type)),
        Lblock(std::make_shared<class_environment>("L")),
        Rblock(std::make_shared<class_environment>("R")),
        Lblock2(std::make_shared<class_environment_var>("L")),
        Rblock2(std::make_shared<class_environment_var>("R"))
{
    log = Logger::setLogger(sim_name,settings::console::verbosity);
    t_eig.set_properties(profile_optimization, 10,"Time optimizing ");
    spin_dimension = HA->get_spin_dimension();
    MPS->initialize(spin_dimension);
}




// We need to make a destructor manually for the enclosing class "class_superblock"
// that encloses "class_hamiltonian_base". Otherwise unique_ptr will forcibly inline its
// own default deleter.
// This allows us to forward declare the abstract base class "class_hamiltonian_base"
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
//class_superblock::~class_superblock()=default;




void class_superblock::clear(){
    log->trace("Copying state to superblock");
    *this = class_superblock(sim_type,sim_name);
}


size_t class_superblock::get_length() const {
    size_t length = Lblock->sites + Rblock->sites + 2;
    assert(length == environment_size + 2 and "ERROR: System length mismatch!");
    return length;

}

size_t class_superblock::get_position() const { return Lblock->sites;}

size_t class_superblock::get_chi() const {return MPS->chiC();}

Eigen::DSizes<long,4> class_superblock::dimensions()const{return MPS->dimensions();}

Eigen::Tensor<Scalar, 4> class_superblock::get_theta() const { return MPS->get_theta();}



//============================================================================//
// Find smallest eigenvalue using Arpack.
//============================================================================//



Eigen::Tensor<Scalar,4> class_superblock::optimize_MPS(Eigen::Tensor<Scalar, 4> &theta, eigutils::eigSetting::Ritz ritz){
    std::array<long,4> shape_theta4 = theta.dimensions();
    std::array<long,4> shape_mpo4   = HA->MPO().dimensions();

    t_eig.tic();
    int nev = 1;
    using namespace settings::precision;
    using namespace eigutils::eigSetting;
    DenseHamiltonianProduct<Scalar>  matrix (Lblock->block.data(), Rblock->block.data(), HA->MPO().data(), HB->MPO().data(), shape_theta4, shape_mpo4);
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
    class_SVD SVD;
    SVD.setThreshold(SVDThreshold);
    auto[U, S, V] = SVD.schmidt(theta,chi_);
    MPS->truncation_error         = SVD.get_truncation_error();
    MPS->LC  = S;
    Eigen::Tensor<Scalar,3> L_U = asDiagonalInversed(MPS->MPS_A->get_L()).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
    Eigen::Tensor<Scalar,3> V_L = V.contract(asDiagonalInversed(MPS->MPS_B->get_L()), idx({2},{0}));
    MPS->MPS_A->set_G(L_U);
    MPS->MPS_B->set_G(V_L);
    return get_theta();
}

void class_superblock::truncate_MPS(const Eigen::Tensor<Scalar, 4> &theta, const std::shared_ptr<class_mps_2site> &MPS_out,long chi_, double SVDThreshold){
    class_SVD SVD;
    SVD.setThreshold(SVDThreshold);
    auto[U, S, V] = SVD.schmidt(theta, chi_);
    MPS_out->truncation_error = SVD.get_truncation_error();
    MPS_out->LC  = S;
    Eigen::Tensor<Scalar,3> L_U = asDiagonalInversed(MPS_out->MPS_A->get_L()).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
    Eigen::Tensor<Scalar,3> V_L = V.contract(asDiagonalInversed(MPS_out->MPS_B->get_L()), idx({2},{0}));
    MPS_out->MPS_A->set_G(L_U);
    MPS_out->MPS_B->set_G(V_L);
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

bool class_superblock::isReal() const {
    bool MPS_isReal         = MPS->isReal();
    bool HA_isReal          = HA->isReal();
    bool HB_isReal          = HB->isReal();
    bool Lblock_isReal      = Lblock->isReal();
    bool Rblock_isReal      = Rblock->isReal();
    bool Lblock2_isReal     = Lblock2->isReal();
    bool Rblock2_isReal     = Rblock2->isReal();

    return
        MPS_isReal
    and HA_isReal
    and HB_isReal
    and Lblock_isReal
    and Rblock_isReal
    and Lblock2_isReal
    and Rblock2_isReal;
}

template<typename T>
Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic> class_superblock::get_H_local_matrix ()const{
    Eigen::Tensor<T,5>tempL;
    Eigen::Tensor<T,5>tempR;
    if constexpr (std::is_same<T,double>::value){
        if(not Lblock->isReal()){throw std::runtime_error("Discarding imaginary data from Lblock when building H_local");}
        if(not Rblock->isReal()){throw std::runtime_error("Discarding imaginary data from Rblock when building H_local");}
        if(not HA->isReal())    {throw std::runtime_error("Discarding imaginary data from MPO A when building H_local");}
        if(not HB->isReal())    {throw std::runtime_error("Discarding imaginary data from MPO B when building H_local");}
        tempL   = Lblock->block.contract(HA->MPO(),Textra::idx({2},{0})).real().shuffle(Textra::array5{4,1,3,0,2}).real();
        tempR   = Rblock->block.contract(HB->MPO(),Textra::idx({2},{1})).real().shuffle(Textra::array5{4,1,3,0,2}).real();
    }else{
        tempL   = Lblock->block.contract(HA->MPO(),Textra::idx({2},{0})).shuffle(Textra::array5{4,1,3,0,2});
        tempR   = Rblock->block.contract(HB->MPO(),Textra::idx({2},{1})).shuffle(Textra::array5{4,1,3,0,2});
    }
    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
    Eigen::Tensor<T,8>H_local = tempL.contract(tempR, Textra::idx({4},{4})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
    return Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>>(H_local.data(),shape,shape);
}

template Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic>                class_superblock::get_H_local_matrix<double>() const;
template Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic>  class_superblock::get_H_local_matrix<std::complex<double>>() const;


template<typename T>
Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic> class_superblock::get_H_local_sq_matrix ()const{
    Eigen::Tensor<T,6>tempL;
    Eigen::Tensor<T,6>tempR;
    if constexpr(std::is_same<T,double>::value){
        if(not Lblock2->isReal()){throw std::runtime_error("Discarding imaginary data from Lblock2 when building H_local_sq");}
        if(not Rblock2->isReal()){throw std::runtime_error("Discarding imaginary data from Rblock2 when building H_local_sq");}
        if(not HA->isReal())     {throw std::runtime_error("Discarding imaginary data from MPO A when building H_local_sq");}
        if(not HB->isReal())     {throw std::runtime_error("Discarding imaginary data from MPO B when building H_local_sq");}
        tempL = Lblock2->block
                .contract(HA->MPO(),Textra::idx({2},{0}))
                .contract(HA->MPO(),Textra::idx({2,5},{0,2}))
                .real()
                .shuffle(Textra::array6{5,1,3,0,2,4});
        tempR = Rblock2->block
                .contract(HB->MPO(),Textra::idx({2},{1}))
                .contract(HB->MPO(),Textra::idx({2,5},{1,2}))
                .real()
                .shuffle(Textra::array6{5,1,3,0,2,4});
    }else{
        tempL = Lblock2->block
                .contract(HA->MPO(),Textra::idx({2},{0}))
                .contract(HA->MPO(),Textra::idx({2,5},{0,2}))
                .shuffle(Textra::array6{5,1,3,0,2,4});
        tempR = Rblock2->block
                .contract(HB->MPO(),Textra::idx({2},{1}))
                .contract(HB->MPO(),Textra::idx({2,5},{1,2}))
                .shuffle(Textra::array6{5,1,3,0,2,4});
    }

    long shape = MPS->chiA() * spin_dimension * MPS->chiB() * spin_dimension;
    Eigen::Tensor<T,8> H_local = tempL.contract(tempR, Textra::idx({4,5},{4,5})).shuffle(Textra::array8{0,1,4,5,2,3,6,7});
    return Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>>(H_local.data(),shape,shape);
}


template Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic>                class_superblock::get_H_local_sq_matrix<double>() const ;
template Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic>  class_superblock::get_H_local_sq_matrix<std::complex<double>>() const;



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
        environment_size = Lblock->sites + Rblock->sites;
    }
    unset_measurements();
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


void class_superblock::unset_measurements() const {
    measurements = Measurements();
    mpstools::common::views::components_computed = false;
}

void class_superblock::do_all_measurements() const {
    using namespace mpstools::common;
    measurements.length                         = measure::length(*this);
    measurements.bond_dimension                 = measure::bond_dimension(*this);
    measurements.norm                           = measure::norm(*this);
    measurements.truncation_error               = measure::truncation_error(*this);
    measurements.energy_mpo                     = measure::energy_mpo(*this);  //This number is needed for variance calculation!
    measurements.energy_per_site_mpo            = measure::energy_per_site_mpo(*this);
    measurements.energy_per_site_ham            = measure::energy_per_site_ham(*this);
    measurements.energy_per_site_mom            = measure::energy_per_site_mom(*this);
    measurements.energy_variance_mpo            = measure::energy_variance_mpo(*this);
    measurements.energy_variance_per_site_mpo   = measure::energy_variance_per_site_mpo(*this);
    measurements.energy_variance_per_site_ham   = measure::energy_variance_per_site_ham(*this);
    measurements.energy_variance_per_site_mom   = measure::energy_variance_per_site_mom(*this);
    measurements.current_entanglement_entropy   = measure::current_entanglement_entropy(*this);
}



void  class_superblock::swap_AB(){
    MPS->swap_AB();
}
