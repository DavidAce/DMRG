//
// Created by david on 2017-11-13.
//

#include "class_finite_chain.h"
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>

using namespace std;
using namespace Textra;
using Scalar = class_finite_chain::Scalar;

class_finite_chain::class_finite_chain(
        int max_length_,
        std::shared_ptr<class_superblock> superblock_,
        std::shared_ptr<class_hdf5_file>  hdf5_,
        SimulationType sim_type_,
        std::string    sim_name_
        )
{
    set_max_length(max_length_);
    set_superblock(std::move(superblock_));
    set_hdf5_file(std::move(hdf5_));
    sim_type = sim_type_;
    sim_name = sim_name_;

}


void class_finite_chain::set_max_length(int max_length_) {
    max_length = max_length_;
    max_length_is_set = true;
}
void class_finite_chain::set_superblock(std::shared_ptr<class_superblock> superblock_) {
    superblock = std::move(superblock_);
    superblock_is_set = true;
}



void class_finite_chain::update_current_length() {
    current_length = ENV_L.back().size + ENV_R.front().size + 2ul;
    assert(current_length == superblock->environment_size + 2ul);
}


int class_finite_chain::insert(){

    if(!max_length_is_set){print_error_and_exit(10);}
    if(!superblock_is_set){print_error_and_exit(11);}
    if(!hdf5_file_is_set){print_error_and_exit(12);}

    assert(ENV_L.size() + ENV_L.size() <= max_length);
    assert(MPS_L.size() + MPS_R.size() <= max_length);
    assert(ENV_L.size()        == superblock->Lblock->size);
    assert(ENV_R.size()        == superblock->Rblock->size);
    assert(ENV2_L.size()       == superblock->Lblock2->size);
    assert(ENV2_R.size()       == superblock->Rblock2->size);


    MPS_L.emplace_back (*superblock->MPS->MPS_A);
    MPS_R.emplace_front(*superblock->MPS->MPS_B);
    MPS_C = superblock->MPS->LC;

    ENV_L.emplace_back(*superblock->Lblock);
    ENV_R.emplace_front(*superblock->Rblock);
    ENV2_L.emplace_back(*superblock->Lblock2);
    ENV2_R.emplace_front(*superblock->Rblock2);
    MPO_L.emplace_back      (superblock->HA->clone());
    MPO_R.emplace_front     (superblock->HB->clone());
    update_current_length();


    int pos = 0;
    for (auto &MPS: MPS_L){MPS.set_position(pos++);}
    for (auto &MPS: MPS_R){MPS.set_position(pos++);}
    pos = 0;
    for (auto &ENV: ENV_L){ENV.set_position(pos++);}
    for (auto &ENV: ENV_R){ENV.set_position(pos++);}
    pos = 0;
    for (auto &ENV2: ENV2_L){ENV2.set_position(pos++);}
    for (auto &ENV2: ENV2_R){ENV2.set_position(pos++);}
    pos = 0;
    for (auto &MPO : MPO_L){MPO->set_position(pos++);}
    for (auto &MPO : MPO_R){MPO->set_position(pos++);}

    superblock->MPS->MPS_A->set_position(MPS_L.back().get_position());
    superblock->MPS->MPS_B->set_position(MPS_R.front().get_position());
    superblock->Lblock->set_position(ENV_L.back().get_position());
    superblock->Rblock->set_position(ENV_R.front().get_position());
    superblock->Lblock2->set_position(ENV2_L.back().get_position());
    superblock->Rblock2->set_position(ENV2_R.front().get_position());
    superblock->HA->set_position(MPO_L.back()->get_position());
    superblock->HB->set_position(MPO_R.front()->get_position());

    assert(ENV_L.back().size + ENV_R.front().size == superblock->environment_size);
    assert(ENV_L.back().size   == superblock->Lblock->size);
    assert(ENV_R.front().size  == superblock->Rblock->size);
    assert(ENV2_L.back().size  == superblock->Lblock2->size);
    assert(ENV2_R.front().size == superblock->Rblock2->size);

    full_mps_has_been_written = false;
    full_mpo_has_been_written = false;

    return (int)MPS_L.size();
}



int class_finite_chain::load(){
    return (int)MPS_L.size();
}



void class_finite_chain::overwrite_local_MPS() {
//    std::cout << "Overwriting positions: "  << MPS_L.size() << " and " <<  MPS_R.size()  << std::endl;
    assert(!MPS_L.empty() and !MPS_R.empty());
    assert(MPS_L.size() + MPS_R.size() <= max_length);
    assert(MPS_L.back().get_position()   == superblock->MPS->MPS_A->get_position());
    assert(MPS_R.front().get_position()  == superblock->MPS->MPS_B->get_position());


    MPS_L.back()    = *superblock->MPS->MPS_A;
    MPS_C           = superblock->MPS->LC;
    MPS_R.front()   = *superblock->MPS->MPS_B;
    full_mps_has_been_written = false;
}


void class_finite_chain::overwrite_local_MPO(){
//    std::cout << "Current state -- Overwrite: " << std::endl;
//    std::cout << "HA: " << superblock->HA->get_position() << " MPO_L back : " << MPO_L.back()->get_position() << std::endl;
//    std::cout << "HB: " << superblock->HB->get_position() << " MPO_R front: " << MPO_R.front()->get_position() << std::endl;
    assert(superblock->HA->get_position() == MPO_L.back()->get_position());
    assert(superblock->HB->get_position() == MPO_R.front()->get_position());
    MPO_L.back()    = superblock->HA->clone();
    MPO_R.front()   = superblock->HB->clone();
    assert(MPO_L.size() + MPO_R.size() == max_length);
    full_mpo_has_been_written = false;
}

void class_finite_chain::overwrite_local_ENV(){
    assert(superblock->Lblock->get_position() == ENV_L.back().get_position());
    assert(superblock->Rblock->get_position() == ENV_R.front().get_position());

    assert(superblock->Lblock2->get_position() == ENV2_L.back().get_position());
    assert(superblock->Rblock2->get_position() == ENV2_R.front().get_position());

    ENV_L.back()    = *superblock->Lblock;
    ENV2_L.back()   = *superblock->Lblock2;

    ENV_R.front()   = *superblock->Rblock;
    ENV2_R.front()  = *superblock->Rblock2;
    assert(ENV_L.size() + ENV_R.size() == max_length);
}


int class_finite_chain::move(){
    //Take current MPS and generate an Lblock one larger and store it in list for later loading
//    std::cout << "Current state -- Direction: " << direction << std::endl;
//    std::cout << "HA: " << superblock->HA->get_position() << " MPO_L back : " << MPO_L.back()->get_position() << std::endl;
//    std::cout << "HB: " << superblock->HB->get_position() << " MPO_R front: " << MPO_R.front()->get_position() << std::endl;
//

    assert(!MPS_L.empty() and !MPS_R.empty());
    assert(MPS_L.size() + MPS_R.size() == max_length);
    assert(ENV_L.size() + ENV_R.size() == max_length);
    assert(ENV_L.back().size + ENV_R.front().size == max_length - 2);
    assert(ENV_L.back().size + ENV_R.front().size == superblock->environment_size);

    assert(MPS_L.back().get_position()   == superblock->MPS->MPS_A->get_position());
    assert(MPS_R.front().get_position()  == superblock->MPS->MPS_B->get_position());


    assert(MPO_L.back()->get_position()  == superblock->HA->get_position() );
    assert(MPO_R.front()->get_position() == superblock->HB->get_position() );

    if (direction == 1){
        //Note that Lblock must just have grown!!
//        assert(MPS_R.size() > 1);
        assert(ENV_L.back().size   + 1  == superblock->Lblock->size);
        assert(ENV2_L.back().size  + 1  == superblock->Lblock2->size);
        assert(ENV_R.front().size       == superblock->Rblock->size);
        assert(ENV2_R.front().size      == superblock->Rblock2->size);
        assert(ENV_L.back().get_position()  + 1 == superblock->Lblock->get_position());
        assert(ENV2_L.back().get_position() + 1 == superblock->Lblock2->get_position());
        assert(ENV_R.front().get_position()     == superblock->Rblock->get_position());
        assert(ENV2_R.front().get_position()    == superblock->Rblock2->get_position());


        MPS_L.emplace_back(*superblock->MPS->MPS_B);
        MPO_L.emplace_back (superblock->HB->clone());
        ENV_L.emplace_back (*superblock->Lblock);
        ENV2_L.emplace_back (*superblock->Lblock2);


        MPS_R.pop_front();
        MPO_R.pop_front();
        ENV_R.pop_front();
        ENV2_R.pop_front();


        superblock->MPS->MPS_A->set_L(superblock->MPS->LC);
        superblock->MPS->MPS_A->set_G(superblock->MPS->MPS_B->get_G());
        superblock->MPS->MPS_A->set_position(MPS_L.back().get_position());
        superblock->MPS->LC = superblock->MPS->MPS_B->get_L();
        superblock->MPS->MPS_B->set_G(MPS_R.front().get_G());
        superblock->MPS->MPS_B->set_L(MPS_R.front().get_L());
        superblock->MPS->MPS_B->set_position(MPS_R.front().get_position());
        MPS_C  = superblock->MPS->LC;

        superblock->HA = MPO_L.back()->clone();
        superblock->HB = MPO_R.front()->clone();

        *superblock->Rblock  = ENV_R.front();
        *superblock->Rblock2 = ENV2_R.front();

    }else{
        //Note that Rblock must just have grown!!
//        assert(MPS_L.size() > 1);
        assert(ENV_R.front().size  + 1  == superblock->Rblock->size);
        assert(ENV2_R.front().size + 1  == superblock->Rblock2->size);
        assert(ENV_L.back().size        == superblock->Lblock->size);
        assert(ENV2_L.back().size       == superblock->Lblock2->size);
        assert(ENV_L.back().get_position()      == superblock->Lblock->get_position());
        assert(ENV2_L.back().get_position()     == superblock->Lblock2->get_position());
        assert(ENV_R.front().get_position()  -1 == superblock->Rblock->get_position());
        assert(ENV2_R.front().get_position() -1 == superblock->Rblock2->get_position());
        MPS_R.emplace_front (*superblock->MPS->MPS_A);
        MPO_R.emplace_front (superblock->HA->clone());
        ENV_R.emplace_front (*superblock->Rblock);
        ENV2_R.emplace_front(*superblock->Rblock2);

        MPS_L.pop_back();
        MPO_L.pop_back();
        ENV_L.pop_back();
        ENV2_L.pop_back();

        superblock->MPS->MPS_B->set_L(superblock->MPS->LC);
        superblock->MPS->MPS_B->set_G(superblock->MPS->MPS_A->get_G());
        superblock->MPS->MPS_B->set_position(MPS_R.front().get_position());
        superblock->MPS->LC = superblock->MPS->MPS_A->get_L();
        superblock->MPS->MPS_A->set_G(MPS_L.back().get_G());
        superblock->MPS->MPS_A->set_L(MPS_L.back().get_L());
        superblock->MPS->MPS_A->set_position(MPS_L.back().get_position());

        MPS_C                       = superblock->MPS->LC;

        superblock->HA = MPO_L.back()->clone();
        superblock->HB = MPO_R.front()->clone();

        *superblock->Lblock  = ENV_L.back();
        *superblock->Lblock2 = ENV2_L.back();
    }

    assert(superblock->MPS->MPS_A->get_G().dimension(2) == superblock->MPS->MPS_B->get_G().dimension(1));
    assert(MPO_L.size() + MPO_R.size() == max_length);
    assert(superblock->HA->get_position() + 1 == MPO_L.size());
    assert(superblock->HA->get_position() + 1 == max_length - MPO_R.size());
    assert(superblock->HB->get_position() + 1 == MPO_L.size() + 1);
    assert(superblock->HB->get_position() + 1 == max_length - MPO_R.size() + 1);
    assert(MPS_L.back().get_position()   == superblock->MPS->MPS_A->get_position());
    assert(MPS_R.front().get_position()  == superblock->MPS->MPS_B->get_position());

    assert(ENV_L.back().get_position()   == superblock->Lblock->get_position());
    assert(ENV_R.front().get_position()  == superblock->Rblock->get_position());
    assert(ENV2_L.back().get_position()  == superblock->Lblock2->get_position());
    assert(ENV2_R.front().get_position() == superblock->Rblock2->get_position());

    assert(MPO_L.back()->get_position()  == superblock->HA->get_position() );
    assert(MPO_R.front()->get_position() == superblock->HB->get_position() );

    //    Check edge
    if (position_is_the_left_edge() or position_is_the_right_edge()) {
        direction *= -1;
    }
    if (position_is_the_left_edge()){
        sweeps++;
    }
    full_mps_has_been_written = false;
    return get_sweeps();
}


void class_finite_chain::print_storage(){
    int i = 0;
    std::cout << setprecision(10);
    std::cout << "Current environment size: " << superblock->environment_size
              << " | Storage length: " << MPS_L.size() + MPS_R.size()
              << " | Particles in environment: " << ENV_L.back().size + ENV_R.front().size
              << std::endl;
    if(!ENV_L.empty()){std::cout << "ENV_L[" <<setw(3) << ENV_L.size()-1 << "]: " << ENV_L.back().block.dimensions() << " Particles: " << ENV_L.back().size << "  <--- Also current environment L" << std::endl;}

    for(auto &it : MPS_L){
        std::cout << "L[" << setw(3) << i  <<  "]: " << it.get_L().dimensions()  << " "
                  << "G[" << setw(3) << i  <<  "]: " << it.get_G().dimensions()  << " ";
        if(&it == &MPS_L.back()){
            std::cout << " <--- Position A";
        }
        std::cout << std::endl;
        i++;
    }
    std::cout << "L[" << setw(3) << '*'  <<  "]: " << MPS_C.dimensions() << "                    <--- Center" << std::endl;
    for(auto &it : MPS_R){
        std::cout << "G[" << setw(3) << i  <<  "]: " << it.get_G().dimensions()  << " "
                  << "L[" << setw(3) << i  <<  "]: " << it.get_L().dimensions()  << " ";
        if(&it == &MPS_R.front()){
            std::cout << " <--- Position B" ;
        }
        std::cout << std::endl;
        i++;
    }
    if(!ENV_R.empty()){std::cout << "ENV_R[" << setw(3) << ENV_R.size()-1 << "]: " << ENV_R.front().block.dimensions() << " Particles: " << ENV_R.front().size << " <--- Also current environment R"  << std::endl;}

}


void class_finite_chain::print_storage_compact(){

    std::cout << setprecision(10);
    std::cout << "Current environment size: " << superblock->environment_size
              << " | Storage length: " << MPS_L.size() + MPS_R.size()
              << " | Particles in environment: " << ENV_L.back().size + ENV_R.front().size
              << std::endl;
    if(!ENV_L.empty()){std::cout << "ENV_L[" <<setw(3) << ENV_L.size()-1 << "]: " << ENV_L.back().block.dimensions() << " Particles: " << ENV_L.back().size << "  <--- Also current environment L" << std::endl;}
    if(!MPS_L.empty()){std::cout << "MPS_L[" <<setw(3) << MPS_L.size()-1 << "]: " << MPS_L.back().get_G().dimensions() <<  "   <--- Also current iteration A" << std::endl;}
    std::cout << "L[" << setw(3) << '*'  <<  "]: " << MPS_C.dimensions() << "                    <--- Center" << std::endl;
    if(!MPS_R.empty()){std::cout << "MPS_R[" <<setw(3) << MPS_R.size()-1 << "]: " << MPS_R.front().get_G().dimensions() << "   <--- Also current iteration B" << std::endl;}
    if(!ENV_R.empty()){std::cout << "ENV_R[" <<setw(3) << ENV_R.size()-1 << "]: " << ENV_R.front().block.dimensions() << " Particles: " << ENV_R.front().size << " <--- Also current environment R"  << std::endl;}
 }

void class_finite_chain::print_hamiltonians() {
    MPO_L.begin()->get()->print_parameter_names();
    for(auto &it : MPO_L){
        it->print_parameter_values();
    }
    for(auto &it : MPO_R){
        it->print_parameter_values();
    }
}

int class_finite_chain::reset_sweeps()  {sweeps = 0; return sweeps;}
int class_finite_chain::get_direction() const {return direction;}
int class_finite_chain::get_sweeps()    const {return sweeps;}
int class_finite_chain::get_length()    const {return (int)(MPS_L.size() + MPS_R.size());}
int class_finite_chain::get_position()  const {return max_length_is_set ? (int)(MPS_L.size() - 1) : 0 ;}
bool class_finite_chain::position_is_the_middle() {
    return max_length_is_set ? (unsigned) get_position() + 1 == (unsigned)(max_length / 2.0) and direction == 1: true ;
}
bool class_finite_chain::position_is_the_middle_any_direction(){
    return max_length_is_set ? (unsigned) get_position() + 1 == (unsigned)(max_length / 2.0) : true ;
}

bool class_finite_chain::position_is_the_left_edge(){
    return get_position() == 0;
}

bool class_finite_chain::position_is_the_right_edge(){
    return get_position() == (int)max_length - 2;
}


void class_finite_chain::print_error_and_exit(int error_type){
    cerr << "Pointers for sweep env_storage has not been set!" << '\n';
    cerr << "try calling class_fDMRG_storage::set_XXXXX(__) before using it." << '\n';
    exit(error_type);
}

