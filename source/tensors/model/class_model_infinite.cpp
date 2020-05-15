//
// Created by david on 2020-05-14.
//

#include "class_model_infinite.h"
#include "class_mpo_factory.h"
#include <general/nmspc_tensor_extra.h>

// We need to make a destructor manually for the enclosing class "class_model_infinite"
// that encloses "class_model_base". Otherwise unique_ptr will forcibly inline its
// own default deleter.
// This allows us to forward declare the abstract base class "class_model_base"
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_model_infinite::~class_model_infinite() = default;

class_model_infinite::class_model_infinite(const class_model_infinite &other) { *this = other; }

class_model_infinite &class_model_infinite::operator=(const class_model_infinite &other) {
    // check for self-assignment
    if(&other == this) return *this;
    // The MPO's are special and the whole point of doing this manually
    HA          = other.HA->clone();
    HB          = other.HB->clone();
    this->cache = other.cache;
    return *this;
}

Eigen::DSizes<long, 4> class_model_infinite::dimensions() const {
    long dim0 = HA->MPO().dimension(0);
    long dim1 = HB->MPO().dimension(1);
    long dim2 = HA->MPO().dimension(2) * HB->MPO().dimension(2);
    long dim3 = HA->MPO().dimension(3) * HB->MPO().dimension(3);
    return Eigen::DSizes<long, 4> { dim0, dim1, dim2, dim3 };
}

const Eigen::Tensor<class_model_infinite::Scalar, 4> &class_model_infinite::get_mpo() const {
    if(cache.mpo) return cache.mpo.value();
    long dim0 = HA->MPO().dimension(0);
    long dim1 = HB->MPO().dimension(1);
    long dim2 = HA->MPO().dimension(2) * HB->MPO().dimension(2);
    long dim3 = HA->MPO().dimension(3) * HB->MPO().dimension(3);
    cache.mpo = HA->MPO().contract(HB->MPO(), Textra::idx({1}, {0})).shuffle(Textra::array6{0, 3, 1, 4, 2, 5}).reshape(Textra::array4{dim0, dim1, dim2, dim3});
    return cache.mpo.value();
}

void class_model_infinite::clear_cache() { cache = Cache(); }
