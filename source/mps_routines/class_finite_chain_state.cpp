//
// Created by david on 2019-01-29.
//


#include "class_finite_chain_state.h"

void class_finite_chain_state::set_max_length(size_t max_length_) {
    max_sites = max_length_;
//    MPS_L.resize(max_sites);
//    MPS_R.resize(max_sites);
//    ENV_L.resize(max_sites);
//    ENV_R.resize(max_sites);
//    ENV2_L.resize(max_sites,);
//    ENV2_R.resize(max_sites);
//    MPO_L.resize(max_sites);
//    MPO_R.resize(max_sites);
    max_length_is_set = true;

}

size_t class_finite_chain_state::get_length()    const {return MPS_L.size() + MPS_R.size();}
size_t class_finite_chain_state::get_position()  const {return max_length_is_set ? MPS_L.size() - 1 : 0 ;}

size_t class_finite_chain_state::get_sweeps()    const {return num_sweeps;}
size_t class_finite_chain_state::reset_sweeps()  {num_sweeps = 0; return num_sweeps;}

int  class_finite_chain_state::get_direction() const {return direction;}
void class_finite_chain_state::flip_direction() {direction *= -1;}


bool class_finite_chain_state::position_is_the_middle() const {
    return max_length_is_set ? (unsigned) get_position() + 1 == (unsigned)(max_sites / 2.0) and direction == 1: true ;
}
bool class_finite_chain_state::position_is_the_middle_any_direction() const {
    return max_length_is_set ? (unsigned) get_position() + 1 == (unsigned)(max_sites / 2.0) : true ;
}

bool class_finite_chain_state::position_is_the_left_edge() const {
    return get_position() == 0;
}

bool class_finite_chain_state::position_is_the_right_edge() const {
    return get_position() == max_sites - 2;
}


void class_finite_chain_state::set_measured_true(){
   energy_has_been_measured   = true;
   variance_has_been_measured = true;
   entropy_has_been_measured  = true;
}
void class_finite_chain_state::set_measured_false(){
    energy_has_been_measured   = false;
    variance_has_been_measured = false;
    entropy_has_been_measured  = false;
}