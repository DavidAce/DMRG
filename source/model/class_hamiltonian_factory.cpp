//
// Created by david on 2018-07-04.
//



#include <general/nmspc_math.h>
#include "class_hamiltonian_factory.h"
#include "class_tf_ising.h"
#include "class_hamiltonian_base.h"
#include "class_hamiltonian_h5tables.h"
#include "class_selfdual_tf_rf_ising.h"


std::unique_ptr<class_hamiltonian_base> class_hamiltonian_factory::create_mpo(std::string model_type_str){

    if (model_type_str == std::string("tf_ising")){
        return std::make_unique<class_tf_ising>();
    }
    else
    if (model_type_str == std::string("tf_nn_ising")){
        return std::make_unique<class_tf_ising>();
    }
    else
    if (model_type_str == std::string("selfdual_tf_rf_ising")){
        return std::make_unique<class_selfdual_tf_rf_ising>();
    }
    else{
        std::cerr << "Wrong model: [ " << model_type_str << " ]" <<std::endl;
        exit(1);
    }
}


//std::unique_ptr<class_hamiltonian_h5table_base> class_hamiltonian_factory::create_table(std::string model_type_str){
//
//    if (model_type_str == std::string("tf_ising")){
//        return std::make_unique<class_tf_ising>();
//    }
//    else
//    if (model_type_str == std::string("tf_nn_ising")){
//        return std::make_unique<class_tf_ising>();
//    }
//    else
//    if (model_type_str == std::string("selfdual_tf_rf_ising")){
//        return std::make_unique<class_selfdual_tf_rf_ising>();
//    }
//    else{
//        std::cerr << "Wrong model: [ " << model_type_str << " ]" <<std::endl;
//        exit(1);
//    }
//}


std::unique_ptr<class_hamiltonian_base> class_hamiltonian_factory::clone(const std::unique_ptr<class_hamiltonian_base> &other){
    return other->clone();
}

