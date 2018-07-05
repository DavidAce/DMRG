//
// Created by david on 2018-07-04.
//



#include "class_hamiltonian_factory.h"
#include "class_tf_ising.h"
#include "class_hamiltonian_base.h"


std::unique_ptr<class_hamiltonian_base> class_hamiltonian_factory::create_mpo(ModelType model_type)
{
    switch (model_type){
        case ModelType::tf_ising:
            using namespace settings::model::tf_ising;
            return std::make_unique<class_tf_ising>();
            break;
        case ModelType::tf_nn_ising:
            using namespace settings::model::tf_nn_ising;
            return std::make_unique<class_tf_ising>();
            break;
    }
}


std::unique_ptr<class_hamiltonian_base> class_hamiltonian_factory::create_mpo(std::string model_type_str){

    if (model_type_str == std::string("tf_ising")){
        return create_mpo(ModelType::tf_ising);
    }
    else
    if (model_type_str == std::string("tf_nn_ising")){
        return create_mpo(ModelType::tf_nn_ising);
    }
    else{
        std::cerr << "Wrong model: [ " << model_type_str << " ]" <<std::endl;
        exit(1);

    }
}

std::unique_ptr<class_hamiltonian_base> class_hamiltonian_factory::clone(const std::unique_ptr<class_hamiltonian_base> &other){
    return other->clone();
}