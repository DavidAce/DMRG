//
// Created by david on 2017-11-13.
//

#include "class_environment_storage.h"
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_mps.h>
#include <mps_routines/class_mpo.h>
#include <IO/class_hdf5_file.h>

using namespace std;
using namespace Textra;
using Scalar = class_environment_storage::Scalar;

class_environment_storage::class_environment_storage(
        int max_length_,
        std::shared_ptr<class_superblock> superblock_,
        std::shared_ptr<class_hdf5_file>  hdf5_)
{
    set_length(max_length_);
    set_superblock(std::move(superblock_));
    set_hdf5_file(std::move(hdf5_));
};




void class_environment_storage::set_length(int max_length_) {
    max_length = max_length_;
    max_length_is_set = true;
}
void class_environment_storage::set_superblock(std::shared_ptr<class_superblock> superblock_) {
    superblock = superblock_;
    superblock_is_set = true;
}

void class_environment_storage::set_hdf5_file(std::shared_ptr<class_hdf5_file> hdf5_){
    hdf5 = hdf5_;
    hdf5_file_is_set = true;
}


int class_environment_storage::insert(){
    if(!max_length_is_set){print_error_and_exit(1);}
    if(!superblock_is_set){print_error_and_exit(1);}
    if(!hdf5_file_is_set){print_error_and_exit(1);}

    MPS_L.push_back (std::make_tuple(superblock->MPS->GA, superblock->MPS->LA, superblock->Lblock->block));
    MPS_R.push_front(std::make_tuple(superblock->MPS->GB, superblock->MPS->LB, superblock->Rblock->block));

    MPO_L.push_back      (superblock->H->M);
    MPO_R.push_front     (superblock->H->M);
    block_Sq_L.push_back (superblock->Lblock2->block);
    block_Sq_R.push_front(superblock->Rblock2->block);
    return (int)MPS_L.size();
}



int class_environment_storage::load(){
    return (int)MPS_L.size();
}

void class_environment_storage::overwrite_MPS() {
    MPS_L.back()   = std::make_tuple(superblock->MPS->GA, superblock->MPS->LA, superblock->Lblock->block);
    MPS_R.front()  = std::make_tuple(superblock->MPS->GB, superblock->MPS->LB, superblock->Rblock->block);

    block_Sq_L.back()  = superblock->Lblock2->block;
    block_Sq_R.front() = superblock->Rblock2->block;

}

int class_environment_storage::move(int &direction, int &sweep){
    //Take current MPS and generate an Lblock one larger and store it in list for later loading
    if (direction == 1){
        //Note that Lblock must just have grown!!
        MPS_L.push_back (std::make_tuple(superblock->MPS->GB, superblock->MPS->LB, superblock->Lblock->block));
        MPS_R.pop_front();

        MPO_L.push_back (superblock->H->M);
        MPO_R.pop_front();
        block_Sq_L.push_back (superblock->Lblock2->block);
        block_Sq_R.pop_front();

    }else{
        //Note that Rblock must just have grown!!
        MPS_R.push_front (std::make_tuple(superblock->MPS->GA, superblock->MPS->LA, superblock->Rblock->block));
        MPS_L.pop_back();

        MPO_R.push_front (superblock->H->M);
        MPO_L.pop_back();
        block_Sq_R.push_front (superblock->Rblock2->block);
        block_Sq_L.pop_back();

    }
    std::tie(superblock->MPS->GA, superblock->MPS->LA, superblock->Lblock->block) = MPS_L.back();
    std::tie(superblock->MPS->GB, superblock->MPS->LB, superblock->Rblock->block) = MPS_R.front();

    superblock->H->M = MPO_L.back();
    superblock->H->M = MPO_R.front();

    superblock->Lblock2->block = block_Sq_L.back();
    superblock->Rblock2->block = block_Sq_R.front();

    if ( MPS_L.size() > 1 ){
        auto next_to_last = std::prev(MPS_L.end(),2);
        superblock->MPS->L_tail = std::get<1>(*next_to_last);
    }else{
        superblock->MPS->L_tail.resize(array1{1});
        superblock->MPS->L_tail.setConstant(1.0);
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


void class_environment_storage::print_storage(){

//    for (auto &G : G_list){
//        cout << G.first  << ": " << G.second.dimensions() << '\n';
//    }
}


void class_environment_storage::print_error_and_exit(int error_type){
    cerr << "Pointers for sweep env_storage has not been set!" << '\n';
    cerr << "try calling class_fDMRG_storage::set_XXXXX(__) before using it." << '\n';
    exit(error_type);
}

