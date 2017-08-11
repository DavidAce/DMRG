//
// Created by david on 7/18/17.
//

#include <class_storage.h>
#include <class_superblock.h>


void class_storage::set_length(int max_length_) {
    max_length = max_length_;
    length_is_set = true;
}


void class_storage::store_insert(const class_superblock &superblock){
    if(!length_is_set){print_error_and_exit(1);}

    position_L = superblock.Lblock.size;
    position_R = max_length - superblock.Rblock.size - 1;

//    cout << "Storing position L: [" << position_L << " R: " << position_R << "] \n";
//    cout << "Lbl: " << superblock.Lblock.block.dimensions() << '\n';
//    cout << "G.A: " << superblock.MPS.GA.dimensions() << '\n';
//    cout << "L.A: " << superblock.MPS.LA.dimensions() << '\n';
//    cout << "G.B: " << superblock.MPS.GB.dimensions() << '\n';
//    cout << "L.B: " << superblock.MPS.LB.dimensions() << '\n';
//    cout << "Rbl: " << superblock.Rblock.block.dimensions() << '\n';

    G_list.insert(std::make_pair(position_L,superblock.MPS.GA));
    G_list.insert(std::make_pair(position_R,superblock.MPS.GB));
    L_list.insert(std::make_pair(position_L,superblock.MPS.LA));
    L_list.insert(std::make_pair(position_R,superblock.MPS.LB));

    Lblock_list[position_L]  = superblock.Lblock;
    Rblock_list[position_R]  = superblock.Rblock;

//    cout << "Block sizes: [" << superblock.Lblock.size << " " << superblock.Rblock.size << "]\n";

}

void class_storage::load(class_superblock &superblock){

//    cout << "Loading position [" << position_L << " " << position_R << "]\n";
    superblock.MPS.GA = G_list.at(position_L);
    superblock.MPS.GB = G_list.at(position_R);

    superblock.MPS.LA = L_list.at(position_L);
    superblock.MPS.LB = L_list.at(position_R);

    superblock.Lblock = Lblock_list.at(position_L);
    superblock.Rblock = Rblock_list.at(position_R);

    superblock.MPS.L_tail =  L_list.at(position_L-1);

//    cout << "Loaded sizes: [" << superblock.Lblock.size << " " << superblock.Rblock.size << "])\n";
//
//    cout << "G.A: " << superblock.MPS.GA.dimensions() << '\n';
//    cout << "G.B: " << superblock.MPS.GB.dimensions() << '\n';
//    cout << "L.A: " << superblock.MPS.LA.dimensions() << '\n';
//    cout << "L.B: " << superblock.MPS.LB.dimensions() << '\n';
//
}

void class_storage::overwrite_MPS(const class_superblock &superblock) {
//    cout << "Overwriting position [" << position_L << " " << position_R << "]\n";
//    cout << "G.A: " << superblock.MPS.GA.dimensions() << '\n';
//    cout << "G.B: " << superblock.MPS.GB.dimensions() << '\n';
//    cout << "L.A: " << superblock.MPS.LA.dimensions() << '\n';
//    cout << "L.B: " << superblock.MPS.LB.dimensions() << '\n';

    G_list[position_L] = superblock.MPS.GA;
    G_list[position_R] = superblock.MPS.GB;

    L_list[position_L] = superblock.MPS.LA;
    L_list[position_R] = superblock.MPS.LB;

}


void class_storage::move(class_superblock &superblock, int &direction, int &sweep){
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


void class_storage::print_storage(){

    for (auto &G : G_list){
        cout << G.first  << ": " << G.second.dimensions() << '\n';
    }
}


void  class_storage::print_error_and_exit(int error_type){
    cout << "Maximum chain length has not been set!" << '\n';
    cout << "try calling class_storage::set_length(int) before using it." << '\n';
    exit(error_type);
}
