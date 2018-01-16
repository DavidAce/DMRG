//
// Created by david on 2017-11-13.
//

#include "class_sweep_storage.h"
#include "class_superblock.h"

using namespace std;
using namespace Textra;

void class_sweep_storage::set_length(int max_length_) {
    max_length = max_length_;
    length_is_set = true;
}

void class_sweep_storage::insert(const class_superblock &superblock){
    if(!length_is_set){print_error_and_exit(1);}

    position_L = superblock.Lblock.size;
    position_R = max_length - superblock.Rblock.size - 1;

//    cout << "Storing position L: [" << position_L << " R: " << position_R << "] \n";
    G_list.insert(std::make_pair(position_L,superblock.MPS.GA));
    G_list.insert(std::make_pair(position_R,superblock.MPS.GB));
    L_list.insert(std::make_pair(position_L,superblock.MPS.LA));
    L_list.insert(std::make_pair(position_R,superblock.MPS.LB));

    Lblock_list[position_L]  = superblock.Lblock;
    Rblock_list[position_R]  = superblock.Rblock;
}


void class_sweep_storage::load(class_superblock &superblock){

//    cout << "Loading position [" << position_L << " " << position_R << "]\n";
    superblock.MPS.GA = G_list.at(position_L);
    superblock.MPS.GB = G_list.at(position_R);

    superblock.MPS.LA = L_list.at(position_L);
    superblock.MPS.LB = L_list.at(position_R);

    superblock.Lblock = Lblock_list.at(position_L);
    superblock.Rblock = Rblock_list.at(position_R);

    superblock.MPS.L_tail =  L_list.at(position_L-1);

}

void class_sweep_storage::overwrite_MPS(const class_superblock &superblock) {

    G_list[position_L] = superblock.MPS.GA;
    G_list[position_R] = superblock.MPS.GB;

    L_list[position_L] = superblock.MPS.LA;
    L_list[position_R] = superblock.MPS.LB;

}

void class_sweep_storage::move(class_superblock &superblock, int &direction, int &sweep){
    //Take current MPS and generate an Lblock one larger and store it in disk for later loading

    position_L += direction;
    position_R += direction;
    if (direction == 1){
        Lblock_list[position_L]  = superblock.Lblock;
    }else{
        Rblock_list[position_R]  = superblock.Rblock;
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


void class_sweep_storage::print_storage(){

    for (auto &G : G_list){
        cout << G.first  << ": " << G.second.dimensions() << '\n';
    }
}


void class_sweep_storage::print_error_and_exit(int error_type){
    cout << "Maximum chain chain_length has not been set!" << '\n';
    cout << "try calling class_fDMRG_storage::set_length(int) before using it." << '\n';
    exit(error_type);
}

