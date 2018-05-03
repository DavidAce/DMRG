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

    MPS_L.push_back (std::make_tuple(superblock->MPS->GA, superblock->MPS->LA));
    MPS_R.push_front(std::make_tuple(superblock->MPS->GB, superblock->MPS->LB));
//    Lblock.push_back (superblock->Lblock->block);
//    Rblock.push_front(superblock->Rblock->block);
    Lblock_list.push_back(*superblock->Lblock);
    Rblock_list.push_front(*superblock->Rblock);
    Lblock2_list.push_back(*superblock->Lblock2);
    Rblock2_list.push_front(*superblock->Rblock2);
    MPO_L.push_back      (*superblock->HA);
    MPO_R.push_front     (*superblock->HB);
    return (int)MPS_L.size();
}



int class_finite_chain_storage::load(){
    return (int)MPS_L.size();
}

void class_finite_chain_storage::overwrite_MPS() {
    MPS_L.back()    = std::make_tuple(superblock->MPS->GA, superblock->MPS->LA);
    MPS_R.front()   = std::make_tuple(superblock->MPS->GB, superblock->MPS->LB);
    Lblock_list.back()   = *superblock->Lblock;
    Rblock_list.front()  = *superblock->Rblock;
    Lblock2_list.back()  = *superblock->Lblock2;
    Rblock2_list.front() = *superblock->Rblock2;

}

int class_finite_chain_storage::move(int &direction, int &sweep){
    //Take current MPS and generate an Lblock one larger and store it in list for later loading
    if (direction == 1){
        //Note that Lblock must just have grown!!
        MPS_L.push_back(std::make_tuple(superblock->MPS->GB, superblock->MPS->LB));
        MPS_R.pop_front();
        Lblock_list.push_back (*superblock->Lblock);
        Rblock_list.pop_front();
        Lblock2_list.push_back (*superblock->Lblock2);
        Rblock2_list.pop_front();
        MPO_L.push_back (*superblock->HA);
        MPO_R.pop_front();


    }else{
        //Note that Rblock must just have grown!!
        MPS_R.push_front (std::make_tuple(superblock->MPS->GA, superblock->MPS->LA));
        MPS_L.pop_back();
        Rblock_list.push_front (*superblock->Rblock);
        Lblock_list.pop_back();
        Rblock2_list.push_front (*superblock->Rblock2);
        Lblock2_list.pop_back();
        MPO_R.push_front (*superblock->HB);
        MPO_L.pop_back();

    }

    std::tie(superblock->MPS->GA, superblock->MPS->LA) = MPS_L.back();
    std::tie(superblock->MPS->GB, superblock->MPS->LB) = MPS_R.front();
    *superblock->Lblock  = Lblock_list.back();
    *superblock->Rblock  = Rblock_list.front();
    *superblock->Lblock2 = Lblock2_list.back();
    *superblock->Rblock2 = Rblock2_list.front();
    *superblock->HA = MPO_L.back();
    *superblock->HB = MPO_R.front();

    if ( MPS_L.size() > 1 ){
        auto next_to_last = std::prev(MPS_L.end(),2);
        superblock->MPS->LB_left = std::get<1>(*next_to_last);
    }else{
        superblock->MPS->LB_left.resize(array1{1});
        superblock->MPS->LB_left.setConstant(1.0);
    }

//    Check edge
    if (MPS_L.size()  <= 1 || MPS_R.size() <= 1) {
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
    for(auto &it : MPS_L){
        std::cout << "G[" << i  << "]: " << std::get<0>(it).dimensions()  << " "
                  << "L[" << i  << "]: " << std::get<1>(it).dimensions()  << " ";
        if(&it == &MPS_L.back()){
            std::cout << " <--- Position A";
        }
        std::cout << std::endl;
        i++;
    }
    for(auto &it : MPS_R){
        std::cout << "G[" << i  << "]: " << std::get<0>(it).dimensions()  << " "
                  << "L[" << i  << "]: " << std::get<1>(it).dimensions()  << " ";
        if(&it == &MPS_R.front()){
            std::cout << " <--- Position B" ;
        }
        std::cout << std::endl;
        i++;
    }
}


void class_finite_chain_storage::print_error_and_exit(int error_type){
    cerr << "Pointers for sweep env_storage has not been set!" << '\n';
    cerr << "try calling class_fDMRG_storage::set_XXXXX(__) before using it." << '\n';
    exit(error_type);
}

