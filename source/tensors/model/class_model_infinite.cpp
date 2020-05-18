//
// Created by david on 2020-05-14.
//

#include <general/nmspc_tensor_extra.h>
#include "class_model_infinite.h"
#include "class_mpo_factory.h"

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
class_model_infinite::class_model_infinite(class_model_infinite &&other)  noexcept = default;               // default move ctor
class_model_infinite &class_model_infinite::operator=(class_model_infinite &&other) noexcept = default;     // default move assign

class_model_infinite::class_model_infinite(const class_model_infinite &other):
    cache(other.cache),
    HA(other.HA->clone()),
    HB(other.HB->clone()),
    model_type(other.model_type),
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


Eigen::DSizes<long, 4> class_model_infinite::dimensions() const {
    long dim0 = HA->MPO().dimension(0);
    long dim1 = HB->MPO().dimension(1);
    long dim2 = HA->MPO().dimension(2) * HB->MPO().dimension(2);
    long dim3 = HA->MPO().dimension(3) * HB->MPO().dimension(3);
    return Eigen::DSizes<long, 4> { dim0, dim1, dim2, dim3 };
}

const class_mpo_base &class_model_infinite::get_mpo_siteA() const{return *HA;}
const class_mpo_base &class_model_infinite::get_mpo_siteB() const{return *HB;}
class_mpo_base &      class_model_infinite::get_mpo_siteA(){return *HA;}
class_mpo_base &      class_model_infinite::get_mpo_siteB(){return *HB;}


const Eigen::Tensor<class_model_infinite::Scalar, 4> &class_model_infinite::get_2site_tensor() const {
    if(cache.twosite_tensor) return cache.twosite_tensor.value();
    long dim0 = HA->MPO().dimension(0);
    long dim1 = HB->MPO().dimension(1);
    long dim2 = HA->MPO().dimension(2) * HB->MPO().dimension(2);
    long dim3 = HA->MPO().dimension(3) * HB->MPO().dimension(3);
    cache.twosite_tensor = HA->MPO().contract(HB->MPO(), Textra::idx({1}, {0})).shuffle(Textra::array6{0, 3, 1, 4, 2, 5}).reshape(Textra::array4{dim0, dim1, dim2, dim3});
    return cache.twosite_tensor.value();
}


void class_model_infinite::clear_cache() { cache = Cache(); }
