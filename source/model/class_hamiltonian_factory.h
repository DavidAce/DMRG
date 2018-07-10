//
// Created by david on 2018-07-04.
//

#ifndef CLASS_HAMILTONIAN_FACTORY_H
#define CLASS_HAMILTONIAN_FACTORY_H

#include <iostream>
#include <memory>
#include "class_hamiltonian_base.h"

class class_hamiltonian_factory{
public:
    static std::unique_ptr<class_hamiltonian_base> create_mpo(std::string model_type_str);
    static std::unique_ptr<class_hamiltonian_base> clone(const std::unique_ptr<class_hamiltonian_base> &other);
};



#endif //CLASS_HAMILTONIAN_FACTORY_H
