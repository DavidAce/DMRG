//
// Created by david on 2017-11-13.
//

#include "class_finite_chain_sweeper.h"
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include "../../cmake-modules/unused/class_mpo.h"

using namespace std;
using namespace Textra;
using Scalar = class_finite_chain_sweeper::Scalar;

class_finite_chain_sweeper::class_finite_chain_sweeper(
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

};




void class_finite_chain_sweeper::set_max_length(int max_length_) {
    max_length = max_length_;
    max_length_is_set = true;
}
void class_finite_chain_sweeper::set_superblock(std::shared_ptr<class_superblock> superblock_) {
    superblock = std::move(superblock_);
    superblock_is_set = true;
}

void class_finite_chain_sweeper::set_hdf5_file(std::shared_ptr<class_hdf5_file> hdf5_){
    hdf5 = std::move(hdf5_);
    hdf5_file_is_set = true;
}

void class_finite_chain_sweeper::update_current_length() {
    current_length = ENV_L.back().size + ENV_R.front().size + 2ul;
    assert(current_length == superblock->environment_size + 2ul);
}

//int class_finite_chain_sweeper::insert_edges(){
//    if(ENV_L.empty() and ENV_R.empty() and ENV2_L.empty() and ENV2_R.empty()) {
//        //Create one
//        ENV_L.emplace_back(*superblock->Lblock);
//        ENV_R.emplace_front(*superblock->Rblock);
//        ENV2_L.emplace_back(*superblock->Lblock2);
//        ENV2_R.emplace_front(*superblock->Rblock2);
//    }else{
//        //Replace existing
//        ENV_L.front()  = *superblock->Lblock;
//        ENV_R.back()   = *superblock->Rblock;
//        ENV2_L.front() = *superblock->Lblock2;
//        ENV2_R.back()  = *superblock->Rblock2;
//    }
//    return 0;
//};




int class_finite_chain_sweeper::insert(){

    if(!max_length_is_set){print_error_and_exit(10);}
    if(!superblock_is_set){print_error_and_exit(11);}
    if(!hdf5_file_is_set){print_error_and_exit(12);}

    assert(ENV_L.size() + ENV_L.size() <= max_length);
    assert(MPS_L.size() + MPS_R.size() <= max_length);
    assert(ENV_L.size()        == superblock->Lblock->size);
    assert(ENV_R.size()        == superblock->Rblock->size);
    assert(ENV2_L.size()       == superblock->Lblock2->size);
    assert(ENV2_R.size()       == superblock->Rblock2->size);


    MPS_L.emplace_back(*superblock->MPS->MPS_A);
    MPS_R.emplace_front(*superblock->MPS->MPS_B);
//
//    MPS_L.emplace_back (std::make_tuple(superblock->MPS->LB_left, superblock->MPS->GA));
//    MPS_R.emplace_front(std::make_tuple(superblock->MPS->GB     , superblock->MPS->LB));
    MPS_C = superblock->MPS->LC;

    ENV_L.emplace_back(*superblock->Lblock);
    ENV_R.emplace_front(*superblock->Rblock);
    ENV2_L.emplace_back(*superblock->Lblock2);
    ENV2_R.emplace_front(*superblock->Rblock2);
    MPO_L.emplace_back      (superblock->HA->clone());
    MPO_R.emplace_front     (superblock->HB->clone());

    update_current_length();

//    std::cout << "Inserted -- New state reflects current superblock: " << std::endl;
//    print_storage();
//    print_hamiltonians();
    assert(ENV_L.back().size + ENV_R.front().size == superblock->environment_size);
    assert(ENV_L.back().size   == superblock->Lblock->size);
    assert(ENV_R.front().size  == superblock->Rblock->size);
    assert(ENV2_L.back().size  == superblock->Lblock2->size);
    assert(ENV2_R.front().size == superblock->Rblock2->size);
    return (int)MPS_L.size();
}



int class_finite_chain_sweeper::load(){
    return (int)MPS_L.size();
}



void class_finite_chain_sweeper::overwrite_local_MPS() {
//    std::cout << "Overwriting positions: "  << MPS_L.size() << " and " <<  MPS_R.size()  << std::endl;
    assert(!MPS_L.empty() and !MPS_R.empty());
    assert(ENV_L.size() + ENV_R.size() <= max_length);
    assert(MPS_L.size() + MPS_R.size() <= max_length);
    assert(ENV_L.back().size   == superblock->Lblock->size);
    assert(ENV_R.front().size  == superblock->Rblock->size);
    assert(ENV2_L.back().size  == superblock->Lblock2->size);
    assert(ENV2_R.front().size == superblock->Rblock2->size);

    MPS_L.back()    = *superblock->MPS->MPS_A;
    MPS_C           = superblock->MPS->LC;
    MPS_R.front()   = *superblock->MPS->MPS_B;

//    MPS_L.back()    = std::make_tuple(superblock->MPS->LB_left, superblock->MPS->GA);
//    MPS_C           = superblock->MPS->LA;
//    MPS_R.front()   = std::make_tuple(superblock->MPS->GB     , superblock->MPS->LB);

}


void class_finite_chain_sweeper::overwrite_local_MPO(){
    MPO_L.back()    = superblock->HA->clone();
    MPO_R.front()   = superblock->HB->clone();
}

void class_finite_chain_sweeper::overwrite_local_ENV(){
    ENV_L.back()    = *superblock->Lblock;
    ENV2_L.back()   = *superblock->Lblock2;

    ENV_R.front()   = *superblock->Rblock;
    ENV2_R.front()  = *superblock->Rblock2;
}


int class_finite_chain_sweeper::move(){
    //Take current MPS and generate an Lblock one larger and store it in list for later loading

//    std::cout << "Starting move in direction: " << direction << " MPS_L size: " << MPS_L.size() << " current e " << superblock->HA->get_energy_reduced()<< std::endl;
//    std::cout << "Starting move in direction: " << direction << " MPS_R size: " << MPS_R.size() << " current e " << superblock->HB->get_energy_reduced()<< std::endl;
    assert(!MPS_L.empty() and !MPS_R.empty());
    assert(MPS_L.size() + MPS_R.size() == max_length);
    assert(ENV_L.size() + ENV_R.size() == max_length);
    assert(ENV_L.back().size + ENV_R.front().size == max_length - 2);
    assert(ENV_L.back().size + ENV_R.front().size == superblock->environment_size);
    if (direction == 1){
        //Note that Lblock must just have grown!!
//        assert(MPS_R.size() > 1);
        assert(ENV_L.back().size   + 1  == superblock->Lblock->size);
        assert(ENV2_L.back().size  + 1  == superblock->Lblock2->size);
        assert(ENV_R.front().size       == superblock->Rblock->size);
        assert(ENV2_R.front().size      == superblock->Rblock2->size);
        MPS_L.emplace_back(*superblock->MPS->MPS_A);
//        MPS_L.emplace_back(std::make_tuple(superblock->MPS->LB_left, superblock->MPS->GA));
        MPO_L.emplace_back (superblock->HA->clone());
        ENV_L.emplace_back (*superblock->Lblock);
        ENV2_L.emplace_back (*superblock->Lblock2);


        MPS_R.pop_front();
        MPO_R.pop_front();
        ENV_R.pop_front();
        ENV2_R.pop_front();


        superblock->MPS->MPS_A->set_L(superblock->MPS->LC);
        superblock->MPS->MPS_A->set_G(superblock->MPS->MPS_B->get_G());
        superblock->MPS->LC = superblock->MPS->MPS_B->get_L();
        superblock->MPS->MPS_B->set_G(MPS_R.front().get_G());
        superblock->MPS->MPS_B->set_L(MPS_R.front().get_L());
        MPS_C  = superblock->MPS->LC;

        *superblock->HA = *superblock->HB;
        *superblock->HB = *MPO_R.front();

        *superblock->Rblock  = ENV_R.front();
        *superblock->Rblock2 = ENV2_R.front();

    }else{
        //Note that Rblock must just have grown!!
//        assert(MPS_L.size() > 1);
        assert(ENV_R.front().size  + 1  == superblock->Rblock->size);
        assert(ENV2_R.front().size + 1  == superblock->Rblock2->size);
        assert(ENV_L.back().size        == superblock->Lblock->size);
        assert(ENV2_L.back().size       == superblock->Lblock2->size);
        MPS_R.emplace_front (*superblock->MPS->MPS_B);
        MPO_R.emplace_front (superblock->HB->clone());
        ENV_R.emplace_front (*superblock->Rblock);
        ENV2_R.emplace_front (*superblock->Rblock2);


        MPS_L.pop_back();
        MPO_L.pop_back();
        ENV_L.pop_back();
        ENV2_L.pop_back();


        superblock->MPS->MPS_B->set_L(superblock->MPS->LC);
        superblock->MPS->MPS_B->set_G(superblock->MPS->MPS_A->get_G());
        superblock->MPS->LC = superblock->MPS->MPS_A->get_L();
        superblock->MPS->MPS_A->set_G(MPS_L.back().get_G());
        superblock->MPS->MPS_A->set_L(MPS_L.back().get_L());
        MPS_C                       = superblock->MPS->LC;
        *superblock->HB = *superblock->HA;
        *superblock->HA = *MPO_L.back();

        *superblock->Lblock  = ENV_L.back();
        *superblock->Lblock2 = ENV2_L.back();
    }
    assert(superblock->MPS->MPS_A->get_G().dimension(2) == superblock->MPS->MPS_B->get_G().dimension(1));
//    assert(superblock->MPS->GA.dimension(2) == superblock->MPS->LA.dimension(0));
//    assert(superblock->MPS->GB.dimension(1) == superblock->MPS->LA.dimension(0));
//    assert(superblock->MPS->GB.dimension(2) == superblock->MPS->LB.dimension(0));

//    Check edge
    if (position_is_the_left_edge() or position_is_the_right_edge()) {
        direction *= -1;
    }
    if (position_is_the_left_edge()){
        sweeps++;
    }
    return get_sweeps();
}


void class_finite_chain_sweeper::write_chain_to_file() {
    unsigned long counter = 0;

    write_list_to_file(get_MPS_L(),sim_name + "/chain/MPS/", counter);
    hdf5->write_dataset(get_MPS_C(),  sim_name + "/chain/MPS/L_" + "_" + to_string(counter++));
    write_list_to_file(get_MPS_R(),sim_name + "/chain/MPS/", counter);

    counter = 0;
    write_list_to_file(get_MPO_L(),sim_name + "/chain/MPO/H", counter);
    write_list_to_file(get_MPO_R(),sim_name + "/chain/MPO/H", counter);

    counter = 0;
    write_list_to_file(get_ENV_L(),sim_name + "/chain/ENV/L", counter);
    counter = 0;
    write_list_to_file(get_ENV_R(),sim_name + "/chain/ENV/R", counter);

    counter = 0;
    write_list_to_file(get_ENV2_L(),sim_name + "/chain/ENV2/L", counter);
    counter = 0;
    write_list_to_file(get_ENV2_R(),sim_name + "/chain/ENV2/R", counter);
}


void class_finite_chain_sweeper::print_storage(){
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


void class_finite_chain_sweeper::print_storage_compact(){

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

void class_finite_chain_sweeper::print_hamiltonians() {
    MPO_L.begin()->get()->print_parameter_names();
    for(auto &it : MPO_L){
        it.get()->print_parameter_values();
    }
    for(auto &it : MPO_R){
        it.get()->print_parameter_values();
    }
}

int class_finite_chain_sweeper::reset_sweeps()  {sweeps = 0; return sweeps;}
int class_finite_chain_sweeper::get_direction() const {return direction;}
int class_finite_chain_sweeper::get_sweeps()    const {return sweeps;}
int class_finite_chain_sweeper::get_length()    const {return (int)(MPS_L.size() + MPS_R.size());}
int class_finite_chain_sweeper::get_position()  const {return max_length_is_set ? (int)(MPS_L.size() - 1) : 0 ;}
bool class_finite_chain_sweeper::position_is_the_middle() {
    return max_length_is_set ? get_position() + 1 == max_length / 2 and direction == 1: true ;
}
bool class_finite_chain_sweeper::position_is_the_left_edge(){
    return get_position() == 0;
}

bool class_finite_chain_sweeper::position_is_the_right_edge(){
    return get_position() == (int)max_length - 2;
}


void class_finite_chain_sweeper::print_error_and_exit(int error_type){
    cerr << "Pointers for sweep env_storage has not been set!" << '\n';
    cerr << "try calling class_fDMRG_storage::set_XXXXX(__) before using it." << '\n';
    exit(error_type);
}

