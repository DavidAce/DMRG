//
// Created by david on 2018-07-04.
//

#ifndef CLASS_HAMILTONIAN_FACTORY_H
#define CLASS_HAMILTONIAN_FACTORY_H

#include <iostream>
#include <memory>
#include "class_hamiltonian_base.h"
#include "class_hamiltonian_h5tables.h"

class class_hamiltonian_factory{
public:
    static std::unique_ptr<class_hamiltonian_base>         create_mpo(size_t position,std::string model_type_str);

    template <typename... T>
    static std::unique_ptr<class_hamiltonian_base>         create_mpo(size_t position,std::string model_type_str, T... args){
        auto mpo = create_mpo(position,model_type_str);
        mpo->set_hamiltonian(args...);
        return mpo;
    }
//    static std::unique_ptr<class_hamiltonian_h5table_base> create_table(std::string model_type_str);
    static std::unique_ptr<class_hamiltonian_base>         clone(std::unique_ptr<class_hamiltonian_base> other);
};



#endif //CLASS_HAMILTONIAN_FACTORY_H
