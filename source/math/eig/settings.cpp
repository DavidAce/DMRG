#include "settings.h"
#include <iostream>


void eig::settings::clear(){
    *this = settings();
}

std::string eig::settings::get_ritz_string() const{
    if(not ritz) throw std::runtime_error("Ritz has not been set");
    return std::string(eig::RitzToString(ritz.value()));}

void eig::settings::checkRitz()
// Checks that ritz is valid for the given problem.
// The valid ritzes are stated in the arpack++ manual page 78.
{
    std::string suggestion;
    if(not ritz) throw std::runtime_error("Ritz has not been set");
    if (type==Type::CPLX or form==Form::NSYM){
        if (ritz.value() == Ritz::LA or
            ritz.value() == Ritz::SA or
            ritz.value() == Ritz::BE)
        {
            std::cerr << "WARNING: Invalid ritz for cplx or nonsym problem: " << RitzToString(ritz.value()) << std::endl;
            if (ritz.value() == Ritz::LA) suggestion = "LR";
            if (ritz.value() == Ritz::SA) suggestion = "SR";
            if (ritz.value() == Ritz::BE) suggestion = "LM";
            std::cerr << "         Suggested ritz: " << suggestion << std::endl;
        }
    }else if (type==Type::REAL and form==Form::SYMM) {
        if (ritz.value() == Ritz::LR or
            ritz.value() == Ritz::SR or
            ritz.value() == Ritz::LI or
            ritz.value() == Ritz::SI)
        {
            std::cerr << "WARNING: Invalid ritz for real and sym problem: " << RitzToString(ritz.value())<< std::endl;
            if (ritz.value() == Ritz::LR) suggestion = "LA";
            if (ritz.value() == Ritz::SR) suggestion = "SA";
            if (ritz.value() == Ritz::LI) suggestion = "LM";
            if (ritz.value() == Ritz::SI) suggestion = "SM";
            std::cerr << "         Suggested ritz: " << suggestion << std::endl;
        }
    }
    if(not suggestion.empty())
        throw std::runtime_error("Invalid ritz");
}