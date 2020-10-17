//
// Created by david on 2020-05-14.
//

#include "class_model_infinite.h"
#include "class_mpo_factory.h"
#include <general/nmspc_tensor_extra.h>
#include <math/num.h>

class_model_infinite::class_model_infinite() = default;

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_model_infinite::~class_model_infinite() = default;                                                    // default dtor
class_model_infinite::class_model_infinite(class_model_infinite &&other) = default;               // default move ctor
class_model_infinite &class_model_infinite::operator=(class_model_infinite &&other) = default;     // default move assign

class_model_infinite::class_model_infinite(const class_model_infinite &other):
    cache(other.cache),
    HA(other.HA->clone()),
    HB(other.HB->clone()),
    model_type(other.model_type)
{}

class_model_infinite &class_model_infinite::operator=(const class_model_infinite &other) {
    // check for self-assignment
    if(this != &other) {
        cache = other.cache;
        HA = other.HA->clone();
        HB = other.HB->clone();
        model_type = other.model_type;
    }
    return *this;
}


void class_model_infinite::initialize(ModelType model_type_){
    tools::log->trace("Initializing model");
    //Generate MPO
    model_type = model_type_;
    HA = class_mpo_factory::create_mpo(0,model_type);
    HB = class_mpo_factory::create_mpo(1,model_type);
}

void class_model_infinite::randomize() {
    tools::log->trace("Randomizing model mpos");
    HA->randomize_hamiltonian();
    HB->randomize_hamiltonian();
    std::vector<class_mpo_site::TableMap> all_params;

    all_params.push_back(HA->get_parameters());
    all_params.push_back(HB->get_parameters());

    HA->set_averages(all_params);
    HB->set_averages(all_params);

}

bool class_model_infinite::is_real() const{
    return HA->is_real() and HB->is_real();
}
bool class_model_infinite::has_nan() const{
    return HA->has_nan() or HB->has_nan();
}
void class_model_infinite::assert_validity() const{
    HA->assert_validity();
    HB->assert_validity();
}


const class_mpo_site &class_model_infinite::get_mpo_siteA() const{return *HA;}
const class_mpo_site &class_model_infinite::get_mpo_siteB() const{return *HB;}
class_mpo_site &      class_model_infinite::get_mpo_siteA(){return *HA;}
class_mpo_site &      class_model_infinite::get_mpo_siteB(){return *HB;}

Eigen::DSizes<long, 4> class_model_infinite::dimensions() const {
    long dim0 = HA->MPO().dimension(0);
    long dim1 = HB->MPO().dimension(1);
    long dim2 = HA->MPO().dimension(2) * HB->MPO().dimension(2);
    long dim3 = HA->MPO().dimension(3) * HB->MPO().dimension(3);
    return Eigen::DSizes<long, 4> { dim0, dim1, dim2, dim3 };
}



bool   class_model_infinite::is_reduced() const{
    return HA->is_reduced() and HB->is_reduced();
}


double class_model_infinite::get_energy_per_site_reduced() const{
    if(not num::all_equal(HA->get_reduced_energy(),HB->get_reduced_energy()))
        throw std::runtime_error(fmt::format("Reduced energy mismatch: HA {:.16f} != HB {:.16f}",HA->get_reduced_energy(),HB->get_reduced_energy()));
    return HA->get_reduced_energy();
}


void   class_model_infinite::set_reduced_energy_per_site(double site_energy){
    HA->set_reduced_energy(site_energy);
    HB->set_reduced_energy(site_energy);
}



const Eigen::Tensor<class_model_infinite::Scalar, 4> &class_model_infinite::get_2site_mpo() const {
    if(cache.twosite_tensor) return cache.twosite_tensor.value();
    long dim0 = HA->MPO().dimension(0);
    long dim1 = HB->MPO().dimension(1);
    long dim2 = HA->MPO().dimension(2) * HB->MPO().dimension(2);
    long dim3 = HA->MPO().dimension(3) * HB->MPO().dimension(3);
    cache.twosite_tensor = HA->MPO().contract(HB->MPO(), Textra::idx({1}, {0})).shuffle(Textra::array6{0, 3, 1, 4, 2, 5}).reshape(Textra::array4{dim0, dim1, dim2, dim3});
    return cache.twosite_tensor.value();
}


void class_model_infinite::clear_cache() { cache = Cache(); }
