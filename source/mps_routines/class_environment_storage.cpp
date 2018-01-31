//
// Created by david on 2017-11-13.
//

#include "class_environment_storage.h"
#include <mps_routines/class_superblock.h>
#include <IO/class_hdf5_file.h>

using namespace std;
using namespace Textra;

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


void class_environment_storage::insert(){
    if(!max_length_is_set){print_error_and_exit(1);}
    if(!superblock_is_set){print_error_and_exit(1);}
    if(!hdf5_file_is_set){print_error_and_exit(1);}

    position_L = superblock->Lblock.size;
    position_R = max_length - superblock->Rblock.size - 1;

//    cout << "Storing position L: [" << position_L << " R: " << position_R << "] \n";
    G_list.insert(std::make_pair(position_L,superblock->MPS.GA));
    G_list.insert(std::make_pair(position_R,superblock->MPS.GB));
    L_list.insert(std::make_pair(position_L,superblock->MPS.LA));
    L_list.insert(std::make_pair(position_R,superblock->MPS.LB));

    Lblock_list[position_L]  = superblock->Lblock;
    Rblock_list[position_R]  = superblock->Rblock;
}


void class_environment_storage::load(){

//    cout << "Loading position [" << position_L << " " << position_R << "]\n";
    superblock->MPS.GA = G_list.at(position_L);
    superblock->MPS.GB = G_list.at(position_R);

    superblock->MPS.LA = L_list.at(position_L);
    superblock->MPS.LB = L_list.at(position_R);

    superblock->Lblock = Lblock_list.at(position_L);
    superblock->Rblock = Rblock_list.at(position_R);

    superblock->MPS.L_tail =  L_list.at(position_L-1);

}

void class_environment_storage::overwrite_MPS() {

    G_list[position_L] = superblock->MPS.GA;
    G_list[position_R] = superblock->MPS.GB;

    L_list[position_L] = superblock->MPS.LA;
    L_list[position_R] = superblock->MPS.LB;

}

void class_environment_storage::move(int &direction, int &sweep){
    //Take current MPS and generate an Lblock one larger and store it in disk for later loading

    position_L += direction;
    position_R += direction;
    if (direction == 1){
        Lblock_list[position_L]  = superblock->Lblock;
    }else{
        Rblock_list[position_R]  = superblock->Rblock;
    }

    //Check edge
    if (position_L <= 1 || position_R >= max_length - 1) {
        direction *= -1;
    }

    //Check if the middle is passed
    if(direction == 1 && position_L == max_length/2 -1 && position_R == max_length/2){
        sweep++;
    }
}


void class_environment_storage::print_storage(){

    for (auto &G : G_list){
        cout << G.first  << ": " << G.second.dimensions() << '\n';
    }
}


void class_environment_storage::print_error_and_exit(int error_type){
    cerr << "Pointers for sweep env_storage has not been set!" << '\n';
    cerr << "try calling class_fDMRG_storage::set_XXXXX(__) before using it." << '\n';
    exit(error_type);
}

