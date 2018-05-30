//
// Created by david on 2017-11-13.
//

#include "class_finite_chain_sweeper.h"
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_mps.h>
#include <mps_routines/class_mpo.h>

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


    MPS_L.emplace_back (std::make_tuple(superblock->MPS->LB_left, superblock->MPS->GA));
    MPS_R.emplace_front(std::make_tuple(superblock->MPS->GB     , superblock->MPS->LB));
    MPS_C = superblock->MPS->LA;

    ENV_L.emplace_back(*superblock->Lblock);
    ENV_R.emplace_front(*superblock->Rblock);
    ENV2_L.emplace_back(*superblock->Lblock2);
    ENV2_R.emplace_front(*superblock->Rblock2);
    MPO_L.emplace_back      (*superblock->HA);
    MPO_R.emplace_front     (*superblock->HB);

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

void class_finite_chain_sweeper::overwrite_MPS() {
//    std::cout << "Overwriting positions: "  << MPS_L.size() << " and " <<  MPS_R.size()  << std::endl;
    assert(!MPS_L.empty() and !MPS_R.empty());
    assert(ENV_L.size() + ENV_R.size() <= max_length);
    assert(MPS_L.size() + MPS_R.size() <= max_length);
    assert(ENV_L.back().size   == superblock->Lblock->size);
    assert(ENV_R.front().size  == superblock->Rblock->size);
    assert(ENV2_L.back().size  == superblock->Lblock2->size);
    assert(ENV2_R.front().size == superblock->Rblock2->size);

    MPS_L.back()    = std::make_tuple(superblock->MPS->LB_left, superblock->MPS->GA);
    MPO_L.back()    = *superblock->HA;
    ENV_L.back()    = *superblock->Lblock;
    ENV2_L.back()   = *superblock->Lblock2;

    MPS_C           = superblock->MPS->LA;

    MPS_R.front()   = std::make_tuple(superblock->MPS->GB     , superblock->MPS->LB);
    MPO_R.front()   = *superblock->HB;
    ENV_R.front()   = *superblock->Rblock;
    ENV2_R.front()  = *superblock->Rblock2;
//    print_storage();
}

int class_finite_chain_sweeper::move(){
    //Take current MPS and generate an Lblock one larger and store it in list for later loading

//    std::cout << "Starting move in direction: " << direction << " MPS_L size: " << MPS_L.size() << " current e " << superblock->HA->get_site_energy()<< std::endl;
//    std::cout << "Starting move in direction: " << direction << " MPS_R size: " << MPS_R.size() << " current e " << superblock->HB->get_site_energy()<< std::endl;
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
        MPS_L.emplace_back(std::make_tuple(superblock->MPS->LB_left, superblock->MPS->GA));
        MPO_L.emplace_back (*superblock->HA);
        ENV_L.emplace_back (*superblock->Lblock);
        ENV2_L.emplace_back (*superblock->Lblock2);


        MPS_R.pop_front();
        MPO_R.pop_front();
        ENV_R.pop_front();
        ENV2_R.pop_front();

        superblock->MPS->LB_left    = superblock->MPS->LA;
        superblock->MPS->GA         = superblock->MPS->GB;
        superblock->MPS->LA         = superblock->MPS->LB;
        superblock->MPS->GB         = std::get<0>(MPS_R.front());
        superblock->MPS->LB         = std::get<1>(MPS_R.front());
        MPS_C                       = superblock->MPS->LA;
        *superblock->HA = *superblock->HB;
        *superblock->HB = MPO_R.front();

        *superblock->Rblock  = ENV_R.front();
        *superblock->Rblock2 = ENV2_R.front();

    }else{
        //Note that Rblock must just have grown!!
//        assert(MPS_L.size() > 1);
        assert(ENV_R.front().size  + 1  == superblock->Rblock->size);
        assert(ENV2_R.front().size + 1  == superblock->Rblock2->size);
        assert(ENV_L.back().size        == superblock->Lblock->size);
        assert(ENV2_L.back().size       == superblock->Lblock2->size);
        MPS_R.emplace_front (std::make_tuple(superblock->MPS->GB, superblock->MPS->LB));
        MPO_R.emplace_front (*superblock->HB);
        ENV_R.emplace_front (*superblock->Rblock);
        ENV2_R.emplace_front (*superblock->Rblock2);


        MPS_L.pop_back();
        MPO_L.pop_back();
        ENV_L.pop_back();
        ENV2_L.pop_back();

        superblock->MPS->LB         = superblock->MPS->LA;
        superblock->MPS->GB         = superblock->MPS->GA;
        superblock->MPS->LA         = superblock->MPS->LB_left;
        superblock->MPS->GA         = std::get<1>(MPS_L.back());
        superblock->MPS->LB_left    = std::get<0>(MPS_L.back());
        MPS_C                       = superblock->MPS->LA;
        *superblock->HB = *superblock->HA;
        *superblock->HA = MPO_L.back();

        *superblock->Lblock  = ENV_L.back();
        *superblock->Lblock2 = ENV2_L.back();
    }
    assert(superblock->MPS->GA.dimension(1) == superblock->MPS->LB_left.dimension(0));
    assert(superblock->MPS->GA.dimension(2) == superblock->MPS->LA.dimension(0));
    assert(superblock->MPS->GB.dimension(1) == superblock->MPS->LA.dimension(0));
    assert(superblock->MPS->GB.dimension(2) == superblock->MPS->LB.dimension(0));

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

    write_list_to_file<1>(get_MPS_L(),sim_name + "/chain/MPS/G", counter);
    write_list_to_file<0>(get_MPS_R(),sim_name + "/chain/MPS/G", counter);

    counter = 0;
    write_list_to_file<0>(get_MPS_L(),sim_name + "/chain/MPS/L", counter);
    hdf5->write_dataset(get_MPS_C(),  sim_name + "/chain/MPS/L_" + "_" + to_string(counter++));
    write_list_to_file<1>(get_MPS_R(),sim_name + "/chain/MPS/L", counter);

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
        std::cout << "L[" << setw(3) << i  <<  "]: " << std::get<0>(it).dimensions()  << " "
                  << "G[" << setw(3) << i  <<  "]: " << std::get<1>(it).dimensions()  << " ";
        if(&it == &MPS_L.back()){
            std::cout << " <--- Position A";
        }
        std::cout << std::endl;
        i++;
    }
    std::cout << "L[" << setw(3) << '*'  <<  "]: " << MPS_C.dimensions() << "                    <--- Center" << std::endl;
    for(auto &it : MPS_R){
        std::cout << "G[" << setw(3) << i  <<  "]: " << std::get<0>(it).dimensions()  << " "
                  << "L[" << setw(3) << i  <<  "]: " << std::get<1>(it).dimensions()  << " ";
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
    if(!MPS_L.empty()){std::cout << "MPS_L[" <<setw(3) << MPS_L.size()-1 << "]: " << std::get<1>(MPS_L.back()).dimensions() <<  "   <--- Also current iteration A" << std::endl;}
    std::cout << "L[" << setw(3) << '*'  <<  "]: " << MPS_C.dimensions() << "                    <--- Center" << std::endl;
    if(!MPS_R.empty()){std::cout << "MPS_R[" <<setw(3) << MPS_R.size()-1 << "]: " << std::get<0>(MPS_R.front()).dimensions() << "   <--- Also current iteration B" << std::endl;}
    if(!ENV_R.empty()){std::cout << "ENV_R[" <<setw(3) << ENV_R.size()-1 << "]: " << ENV_R.front().block.dimensions() << " Particles: " << ENV_R.front().size << " <--- Also current environment R"  << std::endl;}
 }

void class_finite_chain_sweeper::print_hamiltonians() {
    int i = 0;
    std::cout << setprecision(10);
    std::cout << setw(16) << left << "MPO"
              << setw(16) << left << "J"
              << setw(16) << left << "g"
              << setw(16) << left << "e"
              << setw(16) << left << "r"
              << std::endl;
    for(auto &it : MPO_L){
        std::cout << setw(16) << left << i
                  << setw(16) << left << it.get_site_coupling()
                  << setw(16) << left << it.get_site_field()
                  << setw(16) << left << it.get_site_energy()
                  << setw(16) << left << it.get_site_random_field();

//        std::cout << "MPO[" << setw(3) << i  << "]: e = " << setw(12) << it.get_site_energy();
        if(&it == &MPO_L.back()){
            std::cout << " <--- Also current HA";
        }
        std::cout << std::endl;
        i++;
    }
    for(auto &it : MPO_R){
        std::cout << setw(16) << left << i
                  << setw(16) << left << it.get_site_coupling()
                  << setw(16) << left << it.get_site_field()
                  << setw(16) << left << it.get_site_energy()
                  << setw(16) << left << it.get_site_random_field();
//        std::cout << "e[" << setw(3) << i  << "] = " << setw(12) << it.get_site_energy();
        if(&it == &MPO_R.front()){
            std::cout << " <--- Also current HB";
        }
        std::cout << std::endl;
        i++;
    }
    std::cout << std::endl;
}


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

