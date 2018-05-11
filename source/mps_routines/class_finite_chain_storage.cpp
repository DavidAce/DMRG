//
// Created by david on 2017-11-13.
//

#include "class_finite_chain_storage.h"
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_mps.h>
#include <mps_routines/class_mpo.h>
#include <IO/class_hdf5_file.h>

using namespace std;
using namespace Textra;
using Scalar = class_finite_chain_storage::Scalar;

class_finite_chain_storage::class_finite_chain_storage(
        int max_length_,
        std::shared_ptr<class_superblock> superblock_,
        std::shared_ptr<class_hdf5_file>  hdf5_)
{
    set_max_length(max_length_);
    set_superblock(std::move(superblock_));
    set_hdf5_file(std::move(hdf5_));
};




void class_finite_chain_storage::set_max_length(int max_length_) {
    max_length = max_length_;
    max_length_is_set = true;
}
void class_finite_chain_storage::set_superblock(std::shared_ptr<class_superblock> superblock_) {
    superblock = superblock_;
    superblock_is_set = true;
}

void class_finite_chain_storage::set_hdf5_file(std::shared_ptr<class_hdf5_file> hdf5_){
    hdf5 = hdf5_;
    hdf5_file_is_set = true;
}

void class_finite_chain_storage::update_current_length() {
    current_length = ENV_L.back().size + ENV_R.front().size;
    assert(current_length == superblock->chain_length);
}

int class_finite_chain_storage::insert_edges(){
    assert(ENV_L.empty() and ENV_R.empty() and ENV2_L.empty() and ENV2_R.empty());
    ENV_L.emplace_back(*superblock->Lblock);
    ENV_R.emplace_front(*superblock->Rblock);
    ENV2_L.emplace_back(*superblock->Lblock2);
    ENV2_R.emplace_front(*superblock->Rblock2);
//    std::cout << "Inserted edge: " << std::endl;
//    print_storage();
//    print_hamiltonian_energies();
    return 0;
};




int class_finite_chain_storage::insert(){

    if(!max_length_is_set){print_error_and_exit(1);}
    if(!superblock_is_set){print_error_and_exit(1);}
    if(!hdf5_file_is_set){print_error_and_exit(1);}

    assert(ENV_L.back().size   + 1  == superblock->Lblock->size);
    assert(ENV_R.front().size  + 1  == superblock->Rblock->size);
    assert(ENV2_L.back().size  + 1  == superblock->Lblock2->size);
    assert(ENV2_R.front().size + 1  == superblock->Rblock2->size);


    MPS_L.emplace_back (std::make_tuple(superblock->MPS->LB_left, superblock->MPS->GA));
    MPS_R.emplace_front(std::make_tuple(superblock->MPS->GB     , superblock->MPS->LB));
    LA = superblock->MPS->LA;

    ENV_L.emplace_back(*superblock->Lblock);
    ENV_R.emplace_front(*superblock->Rblock);
    ENV2_L.emplace_back(*superblock->Lblock2);
    ENV2_R.emplace_front(*superblock->Rblock2);
    MPO_L.emplace_back      (*superblock->HA);
    MPO_R.emplace_front     (*superblock->HB);
//    std::cout << "Inserted -- New state reflects current superblock: " << std::endl;
//    print_storage();
//    print_hamiltonian_energies();
    assert(ENV_L.back().size + ENV_R.front().size == superblock->chain_length);
    assert(ENV_L.back().size   == superblock->Lblock->size);
    assert(ENV_R.front().size  == superblock->Rblock->size);
    assert(ENV2_L.back().size  == superblock->Lblock2->size);
    assert(ENV2_R.front().size == superblock->Rblock2->size);
    return (int)MPS_L.size();
}



int class_finite_chain_storage::load(){
    return (int)MPS_L.size();
}

void class_finite_chain_storage::overwrite_MPS() {
//    std::cout << "Overwriting positions: "  << MPS_L.size() << " and " <<  MPS_R.size()  << std::endl;
    assert(!MPS_L.empty() and !MPS_R.empty());
    assert(ENV_L.back().size   == superblock->Lblock->size);
    assert(ENV_R.front().size  == superblock->Rblock->size);
    assert(ENV2_L.back().size  == superblock->Lblock2->size);
    assert(ENV2_R.front().size == superblock->Rblock2->size);

    MPS_L.back()    = std::make_tuple(superblock->MPS->LB_left, superblock->MPS->GA);
    MPO_L.back()    = *superblock->HA;
    ENV_L.back()    = *superblock->Lblock;
    ENV2_L.back()   = *superblock->Lblock2;

    LA              = superblock->MPS->LA;

    MPS_R.front()   = std::make_tuple(superblock->MPS->GB     , superblock->MPS->LB);
    MPO_R.front()   = *superblock->HB;
    ENV_R.front()   = *superblock->Rblock;
    ENV2_R.front()  = *superblock->Rblock2;
//    print_storage();
}

int class_finite_chain_storage::move(int &direction, int &sweep){
    //Take current MPS and generate an Lblock one larger and store it in list for later loading

//    std::cout << "Starting move in direction: " << direction << " MPS_L size: " << MPS_L.size() << " current e " << superblock->HA->get_site_energy()<< std::endl;
//    std::cout << "Starting move in direction: " << direction << " MPS_R size: " << MPS_R.size() << " current e " << superblock->HB->get_site_energy()<< std::endl;
    assert(!MPS_L.empty() and !MPS_R.empty());
    assert(MPS_L.size() + MPS_R.size() == max_length);
    assert(ENV_L.size() + ENV_R.size() == max_length + 2);
    assert(ENV_L.back().size + ENV_R.front().size == superblock->chain_length);

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

        LA                          = superblock->MPS->LA;
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

        LA                          = superblock->MPS->LA;
        *superblock->HB = *superblock->HA;
        *superblock->HA = MPO_L.back();

        *superblock->Lblock  = ENV_L.back();
        *superblock->Lblock2 = ENV2_L.back();
    }
    assert(superblock->MPS->GA.dimension(1) == superblock->MPS->LB_left.dimension(0));
    assert(superblock->MPS->GA.dimension(2) == superblock->MPS->LA.dimension(0));
    assert(superblock->MPS->GB.dimension(1) == superblock->MPS->LA.dimension(0));
    assert(superblock->MPS->GB.dimension(2) == superblock->MPS->LB.dimension(0));

//    std::cout << "Finished move in direction: " << direction << " MPS_L size: " << MPS_L.size() << " empty: "<< MPS_L.empty() << " current e " << superblock->HA->get_site_energy()<< std::endl;
//    std::cout << "Finished move in direction: " << direction << " MPS_R size: " << MPS_R.size() << " empty: "<< MPS_R.empty() << " current e " << superblock->HB->get_site_energy()<< std::endl;
//    std::cout << "Current state: " << std::endl;
//    print_storage();
//    print_hamiltonian_energies();
//    std::cout << std::endl;

//    Check edge
    if (MPS_L.size() == 1 or MPS_R.size() == 1) {
        direction *= -1;
    }

    //Check if the middle is passed
    if(direction == 1 && (int)MPS_L.size() == max_length/2  -1){
        sweep++;
    }
    return (int)MPS_L.size();

}


void class_finite_chain_storage::print_storage(){
//    auto it1 = MPS_L.begin();
//    auto it2 = MPS_
    int i = 0;
    std::cout << setprecision(10);
    std::cout << "Current chain length: " << superblock->chain_length
              << " | Storage length: " << MPS_L.size() + MPS_R.size()
              << " | Particles in environment: " << ENV_L.back().size + ENV_R.front().size
              << std::endl;
    if(!ENV_L.empty()){std::cout << "ENV_L[" <<setw(3) << ENV_L.size()-1 << "]: " << ENV_L.back().block.dimensions() << " Particles: " << ENV_L.back().size << "  <--- Also current environment L" << std::endl;}

    for(auto &it : MPS_L){
        std::cout << "L[" << setw(3) << i  <<  "]: " << std::get<0>(it).dimensions()  << " "
                  << "G[" << setw(3) << i  <<  "]: " << std::get<1>(it).dimensions()  << " ";
        if(&it == &MPS_L.back()){
            std::cout << " <--- Position A";
//            std::cout << " -- ENV_L[" << i << "] size: " << ENV_L.back().size ;
        }
        std::cout << std::endl;
        i++;
    }
    std::cout << "L[" << setw(3) << '*'  <<  "]: " << LA.dimensions() << "                    <--- Center" << std::endl;
    for(auto &it : MPS_R){
        std::cout << "G[" << setw(3) << i  <<  "]: " << std::get<0>(it).dimensions()  << " "
                  << "L[" << setw(3) << i  <<  "]: " << std::get<1>(it).dimensions()  << " ";
        if(&it == &MPS_R.front()){
            std::cout << " <--- Position B" ;
//            std::cout << " -- ENV_R[" << i << "] size: " << ENV_R.front().size ;
        }
        std::cout << std::endl;
        i++;
    }
    if(!ENV_R.empty()){std::cout << "ENV_R[" <<setw(3) << ENV_R.size()-1 << "]: " << ENV_R.front().block.dimensions() << " Particles: " << ENV_R.front().size << " <--- Also current environment R"  << std::endl;}

}


void class_finite_chain_storage::print_storage_compact(){
//    auto it1 = MPS_L.begin();
//    auto it2 = MPS_
//    int i = 0;
    std::cout << setprecision(10);
    std::cout << "Current chain length: " << superblock->chain_length
              << " | Storage length: " << MPS_L.size() + MPS_R.size()
              << " | Particles in environment: " << ENV_L.back().size + ENV_R.front().size
              << std::endl;
    if(!ENV_L.empty()){std::cout << "ENV_L[" <<setw(3) << ENV_L.size()-1 << "]: " << ENV_L.back().block.dimensions() << " Particles: " << ENV_L.back().size << "  <--- Also current environment L" << std::endl;}
    if(!MPS_L.empty()){std::cout << "MPS_L[" <<setw(3) << MPS_L.size()-1 << "]: " << std::get<1>(MPS_L.back()).dimensions() <<  "   <--- Also current position A" << std::endl;}
    std::cout << "L[" << setw(3) << '*'  <<  "]: " << LA.dimensions() << "                    <--- Center" << std::endl;
    if(!MPS_R.empty()){std::cout << "MPS_R[" <<setw(3) << MPS_R.size()-1 << "]: " << std::get<0>(MPS_R.front()).dimensions() << "   <--- Also current position B" << std::endl;}
    if(!ENV_R.empty()){std::cout << "ENV_R[" <<setw(3) << ENV_R.size()-1 << "]: " << ENV_R.front().block.dimensions() << " Particles: " << ENV_R.front().size << " <--- Also current environment R"  << std::endl;}
 }

void class_finite_chain_storage::print_hamiltonian_energies() {
    int i = 0;
    std::cout << setprecision(10);
    for(auto &it : MPO_L){
        std::cout << "e[" << setw(3) << i  << "] = " << setw(12) << it.get_site_energy();
        if(&it == &MPO_L.back()){
            std::cout << " <--- Also current HA";
        }
        std::cout << std::endl;
        i++;
    }


//    std::cout << "e[" << setw(3) << i++  << "] = " << setw(12) << superblock->HA->get_site_energy() << std::endl;
//    std::cout << "e[" << setw(3) << i++  << "] = " << setw(12) << superblock->HB->get_site_energy() << std::endl;

    for(auto &it : MPO_R){
        std::cout << "e[" << setw(3) << i  << "] = " << setw(12) << it.get_site_energy();
        if(&it == &MPO_R.front()){
            std::cout << " <--- Also current HB";
        }
        std::cout << std::endl;
        i++;
    }
    std::cout << std::endl;
}


void class_finite_chain_storage::print_error_and_exit(int error_type){
    cerr << "Pointers for sweep env_storage has not been set!" << '\n';
    cerr << "try calling class_fDMRG_storage::set_XXXXX(__) before using it." << '\n';
    exit(error_type);
}

