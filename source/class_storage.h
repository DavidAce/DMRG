//
// Created by david on 7/18/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_STORAGE_H
#define FINITE_DMRG_EIGEN_CLASS_STORAGE_H

#include <list>
#include <vector>
#include <n_tensor_extra.h>
#include <class_superblock.h>
#include <map>

using namespace Textra;
using namespace std;



class class_storage {
private:
    template <typename list_type, typename inType>
    void insert_middle(list_type &target_list, const inType &elem){
        typename list_type::iterator it = target_list.begin();
        std::advance(it, std::distance(target_list.begin(), target_list.end())/2);
        target_list.insert(it, elem);
    }
    template <typename list_type, typename inType>
    void replace(list_type &target_list, const inType &elem, const long at){
        typename list_type::iterator it = target_list.begin();
        std::advance(it, std::distance(target_list.begin(), target_list.begin()+at));
        *it = elem;
    }
public:

    std::map<size_t, Tensor3> G_list;
    std::map<size_t, Tensor1> L_list;

//    std::vector<Tensor3> Block_list_old;
//    std::vector<std::tuple<Tensor3,long ,std::string>> Block_list;
    std::map<size_t,class_environment_L> Lblock_list;
    std::map<size_t,class_environment_R> Rblock_list;

    const size_t max_length;
    size_t position_L = 0;
    size_t position_R = max_length - 1;
    class_storage(size_t L):max_length(L){
    };



    void store_insert(const class_superblock &superblock){


        position_L = superblock.Lblock.size;
        position_R = max_length - superblock.Lblock.size - 1;

        cout << "Storing position L: " << position_L << " R: " << position_R << "]" << endl;
        cout << "Lbl: " << superblock.Lblock.block.dimensions() << endl;
        cout << "G.A: " << superblock.MPS.G.A.dimensions() << endl;
        cout << "L.A: " << superblock.MPS.L.A.dimensions() << endl;
        cout << "G.B: " << superblock.MPS.G.B.dimensions() << endl;
        cout << "L.B: " << superblock.MPS.L.B.dimensions() << endl;
        cout << "Rbl: " << superblock.Rblock.block.dimensions() << endl;

        G_list.insert(std::make_pair(position_L,superblock.MPS.G.A));
        G_list.insert(std::make_pair(position_R,superblock.MPS.G.B));
        L_list.insert(std::make_pair(position_L,superblock.MPS.L.A));
        L_list.insert(std::make_pair(position_R,superblock.MPS.L.B));

        Lblock_list[position_L]  = superblock.Lblock;
        Rblock_list[position_R]  = superblock.Rblock;

        cout << "Block sizes: [" << superblock.Lblock.size << " " << superblock.Rblock.size << "]" << endl;

    }

    void load(class_superblock &superblock){

        cout << "Loading position [" << position_L << " " << position_R << "]" << endl;
        superblock.MPS.G.A = G_list.at(position_L);
        superblock.MPS.G.B = G_list.at(position_R);

        superblock.MPS.L.A = L_list.at(position_L);
        superblock.MPS.L.B = L_list.at(position_R);

        superblock.Lblock = Lblock_list.at(position_L);
        superblock.Rblock = Rblock_list.at(position_R);

        superblock.MPS.L_tail =  L_list.at(position_L-1);

        cout << "Loaded sizes: [" << superblock.Lblock.size << " " << superblock.Rblock.size << "])" << endl;

        cout << "G.A: " << superblock.MPS.G.A.dimensions() << endl;
        cout << "G.B: " << superblock.MPS.G.B.dimensions() << endl;
        cout << "L.A: " << superblock.MPS.L.A.dimensions() << endl;
        cout << "L.B: " << superblock.MPS.L.B.dimensions() << endl;


    }

    void overwrite(const class_superblock &superblock){
        cout << "Overwriting position [" << position_L << " " << position_R << "]" << endl;
        cout << "G.A: " << superblock.MPS.G.A.dimensions() << endl;
        cout << "G.B: " << superblock.MPS.G.B.dimensions() << endl;
        cout << "L.A: " << superblock.MPS.L.A.dimensions() << endl;
        cout << "L.B: " << superblock.MPS.L.B.dimensions() << endl;

        G_list[position_L] = superblock.MPS.G.A;
        G_list[position_R] = superblock.MPS.G.B;

        L_list[position_L] = superblock.MPS.L.A;
        L_list[position_R] = superblock.MPS.L.B;

        cout << "Written block sizes: [" << superblock.Lblock.size << " " << superblock.Rblock.size << "]" << endl;
    }


    void move(class_superblock &superblock, const long direction){
        //Take current MPS and generate an Lblock one larger and store it in disk for later loading

        superblock.enlarge_environment(direction);
        position_L += direction;
        position_R += direction;
        if (direction == 1){
            Lblock_list[position_L]  = superblock.Lblock;
        }else{
            Rblock_list[position_R]  = superblock.Rblock;
        }





    }


    void print_storage(){

        for (auto &G : G_list){
            cout << G.first  << ": " << G.second.dimensions() << endl;
        }
    }
};


#endif //FINITE_DMRG_EIGEN_CLASS_STORAGE_H
