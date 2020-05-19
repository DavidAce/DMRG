//
// Created by david on 2018-07-04.
//

#include "class_mpo_factory.h"
#include "class_ising_sdual.h"
#include "class_ising_tf_rf.h"
#include "class_mpo_site.h"

std::unique_ptr<class_mpo_site> class_mpo_factory::create_mpo(size_t position, ModelType model_type) {
    switch(model_type){
        case ModelType::ising_tf_rf: return std::make_unique<class_ising_tf_rf>(model_type,position);
        case ModelType::ising_sdual: return std::make_unique<class_ising_sdual>(model_type,position);
        default: throw std::runtime_error(fmt::format("Wrong model type: [{}]", enum2str(model_type)));
    }
}

std::unique_ptr<class_mpo_site> class_mpo_factory::clone(std::unique_ptr<class_mpo_site> other) { return other->clone(); }
