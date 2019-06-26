//
// Created by david on 2018-07-04.
//



#include <general/nmspc_math.h>
#include "class_hamiltonian_factory.h"
#include "class_tf_ising.h"
#include "class_hamiltonian_base.h"
#include "class_hamiltonian_h5tables.h"
#include "class_selfdual_tf_rf_ising.h"

std::unique_ptr<class_hamiltonian_base> class_hamiltonian_factory::create_mpo(size_t position, std::string model_type_str){

    if (model_type_str == std::string("tf_ising")){
        return std::make_unique<class_tf_ising>(position,model_type_str);
    }
    else
    if (model_type_str == std::string("tf_nn_ising")){
        return std::make_unique<class_tf_ising>(position,model_type_str);
    }
    else
    if (model_type_str == std::string("selfdual_tf_rf_ising")){
        return std::make_unique<class_selfdual_tf_rf_ising>(position,model_type_str);
    }
    else{
        throw std::runtime_error("Wrong model: [ "  + model_type_str + " ]");
    }
}



std::unique_ptr<class_hamiltonian_base> class_hamiltonian_factory::clone(std::unique_ptr<class_hamiltonian_base> other){
    return other->clone();
}

