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
    set_length(max_length_);
    set_superblock(std::move(superblock_));
    set_hdf5_file(std::move(hdf5_));
};




void class_finite_chain_storage::set_length(int max_length_) {
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


int class_finite_chain_storage::insert(){
    if(!max_length_is_set){print_error_and_exit(1);}
    if(!superblock_is_set){print_error_and_exit(1);}
    if(!hdf5_file_is_set){print_error_and_exit(1);}

    MPS_L.emplace_back (std::make_tuple(superblock->MPS->LB_left, superblock->MPS->GA));
    MPS_R.emplace_front(std::make_tuple(superblock->MPS->GB     , superblock->MPS->LB));
    LA = superblock->MPS->LA;
    ENV_L.emplace_back(*superblock->Lblock);
    ENV_R.emplace_front(*superblock->Rblock);
    ENV2_L.emplace_back(*superblock->Lblock2);
    ENV2_R.emplace_front(*superblock->Rblock2);
    MPO_L.emplace_back      (*superblock->HA);
    MPO_R.emplace_front     (*superblock->HB);
    print_hamiltonian_energies();
    print_storage();
    return (int)MPS_L.size();
}



int class_finite_chain_storage::load(){
    return (int)MPS_L.size();
}

void class_finite_chain_storage::overwrite_MPS() {
    std::cout << "Overwriting positions: " << MPS_L.size() << " and " << MPS_R.size()  << std::endl;

    MPS_L.back()    = std::make_tuple(superblock->MPS->LB_left, superblock->MPS->GA);
    MPS_R.front()   = std::make_tuple(superblock->MPS->GB     , superblock->MPS->LB);
    LA = superblock->MPS->LA;
    ENV_L.back()   = *superblock->Lblock;
    ENV_R.front()  = *superblock->Rblock;
    ENV2_L.back()  = *superblock->Lblock2;
    ENV2_R.front() = *superblock->Rblock2;
    MPO_L.back()    = *superblock->HA;
    MPO_R.front()   = *superblock->HB;
}

int class_finite_chain_storage::move(int &direction, int &sweep){
    //Take current MPS and generate an Lblock one larger and store it in list for later loading

    std::cout << "Starting move in direction: " << direction << " MPS_L size: " << MPS_L.size() << " empty: "<< MPS_L.empty() << " current e " << superblock->HA->get_site_energy()<< std::endl;
    std::cout << "Starting move in direction: " << direction << " MPS_R size: " << MPS_R.size() << " empty: "<< MPS_R.empty() << " current e " << superblock->HB->get_site_energy()<< std::endl;

//    if (MPS_L.size() <= 1 or MPS_R.size() <= 1) {
////        direction *= -1;
//        MPO_L.back() = *superblock->HA;
//        MPO_R.front() = *superblock->HB;
//    }
    if (direction == 1){
        //Note that Lblock must just have grown!!

        MPS_L.emplace_back(std::make_tuple(superblock->MPS->LB_left, superblock->MPS->GA));
        MPO_L.emplace_back (*superblock->HA);
        ENV_L.emplace_back (*superblock->Lblock);
        ENV2_L.emplace_back (*superblock->Lblock2);


        MPS_R.pop_front();
        MPO_R.pop_front();
        ENV_R.pop_front();
        ENV2_R.pop_front();

        std::tie(superblock->MPS->LB_left,superblock->MPS->GA)  = std::make_tuple(superblock->MPS->LA, superblock->MPS->GB);
        std::tie(superblock->MPS->GB     , superblock->MPS->LB) = MPS_R.front();

        LA  = superblock->MPS->LB;
        superblock->MPS->LA = LA;
        *superblock->HA = *superblock->HB;
        *superblock->HB = MPO_R.front();

//        MPO_L.back() = *superblock->HA;

        *superblock->Rblock  = ENV_R.front();
        *superblock->Rblock2 = ENV2_R.front();



    }else{
        //Note that Rblock must just have grown!!
        MPS_R.emplace_front (std::make_tuple(superblock->MPS->GB, superblock->MPS->LB));
        MPO_R.emplace_front (*superblock->HB);
        ENV_R.emplace_front (*superblock->Rblock);
        ENV2_R.emplace_front (*superblock->Rblock2);


        MPS_L.pop_back();
        MPO_L.pop_back();
        ENV_L.pop_back();
        ENV2_L.pop_back();

        std::tie(superblock->MPS->GB     , superblock->MPS->LB) = std::make_tuple(superblock->MPS->GB, superblock->MPS->LB);
        std::tie(superblock->MPS->LB_left,superblock->MPS->GA)  = MPS_L.back();


        LA  = superblock->MPS->LB_left;
        superblock->MPS->LA = LA;
        *superblock->HB = *superblock->HA;
        *superblock->HA = MPO_L.back();
//        MPO_L.back() = *superblock->HB;

        *superblock->Lblock  = ENV_L.back();
        *superblock->Lblock2 = ENV2_L.back();


    }

    std::cout << "Finished move in direction: " << direction << " MPS_L size: " << MPS_L.size() << " empty: "<< MPS_L.empty() << " current e " << superblock->HA->get_site_energy()<< std::endl;
    std::cout << "Finished move in direction: " << direction << " MPS_R size: " << MPS_R.size() << " empty: "<< MPS_R.empty() << " current e " << superblock->HB->get_site_energy()<< std::endl;
    print_hamiltonian_energies();
    print_storage();
//    Check edge
    if (MPS_L.size() <= 1 or MPS_R.size() <= 1) {
        direction *= -1;
//        MPO_L.back() = *superblock->HA;
//        MPO_R.front() = *superblock->HB;
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
    for(auto &it : MPS_L){
        std::cout << "L[" << setw(3) << i  <<  "]: " << std::get<0>(it).dimensions()  << " "
                  << "G[" << setw(3) << i  <<  "]: " << std::get<1>(it).dimensions()  << " ";
        if(&it == &MPS_L.back()){
            std::cout << " <--- Position A";
            std::cout << " -- Left block size: " << ENV_L.size() ;
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
            std::cout << " -- Right block size: " << ENV_R.size();
        }
        std::cout << std::endl;
        i++;
    }

}

void class_finite_chain_storage::print_hamiltonian_energies() {
    int i = 0;
    std::cout << setprecision(10);
    for(auto &it : MPO_L){
        std::cout << "e[" << setw(3) << i  << "] = " << setw(12) << it.get_site_energy();
        if(&it == &MPO_L.back()){
            std::cout << " <--- Position A  -- Actual e = " << superblock->HA->get_site_energy();
        }
        std::cout << std::endl;
        i++;
    }
    for(auto &it : MPO_R){
        std::cout << "e[" << setw(3) << i  << "] = " << setw(12) << it.get_site_energy();
        if(&it == &MPO_R.front()){
            std::cout << " <--- Position B  -- Actual e = " << superblock->HB->get_site_energy(); ;
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

