#include "settings.h"
#include <string>

void eig::settings::clear() { *this = settings(); }

std::string eig::settings::get_ritz_string() const {
    if(not ritz) throw std::runtime_error("Ritz has not been set");
    return std::string(eig::RitzToString(ritz.value()));
}

void eig::settings::checkRitz()
// Checks that ritz is valid for the given problem.
// The valid ritzes are stated in the arpack++ manual page 78.
{
    std::string error_mesg;
    if(not ritz) throw std::runtime_error("Ritz has not been set");
    if(type == Type::CPLX or form == Form::NSYM) {
        if(ritz.value() == Ritz::LA or ritz.value() == Ritz::SA or ritz.value() == Ritz::BE) {
            error_mesg.append("Invalid ritz for cplx or nonsym problem: ").append(RitzToString(ritz.value()));
            if(ritz.value() == Ritz::LA) error_mesg.append(" | Suggested ritz: LR");
            if(ritz.value() == Ritz::SA) error_mesg.append(" | Suggested ritz: SR");
            if(ritz.value() == Ritz::BE) error_mesg.append(" | Suggested ritz: LM");
        }
    } else if(type == Type::REAL and form == Form::SYMM) {
        if(ritz.value() == Ritz::LR or ritz.value() == Ritz::SR or ritz.value() == Ritz::LI or ritz.value() == Ritz::SI) {
            error_mesg.append("Invalid ritz for real and sym problem: ").append(RitzToString(ritz.value()));
            if(ritz.value() == Ritz::LR) error_mesg.append(" | Suggested ritz: LA");
            if(ritz.value() == Ritz::SR) error_mesg.append(" | Suggested ritz: SA");
            if(ritz.value() == Ritz::LI) error_mesg.append(" | Suggested ritz: LM");
            if(ritz.value() == Ritz::SI) error_mesg.append(" | Suggested ritz: SM");
        }
    }
    if(not error_mesg.empty()) throw std::runtime_error(error_mesg);
}