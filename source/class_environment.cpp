//
// Created by david on 7/21/17.
//

#include "class_environment.h"

class_environment_L::class_environment_L(){
    size = 0;
    block.resize(1, 1, 3);
    block.setZero();
    block(0, 0, 0) = 1;
}


class_environment_R::class_environment_R(){
    size = 0;
    block.resize(1, 1, 3);
    block.setZero();
    block(0, 0, 2) = 1;
}


void class_environment_L::enlarge(const class_MPS &MPS, const Tensor4 &W) {

    picture.append(single_picture);
    size++;
    Tensor3 block_enlarged = block.contract(asDiagonal(MPS.L_tail), idxlist1{idx2(0,0)})
            .contract(MPS.G.A, idxlist1{idx2(2,1)})
            .contract(W, idxlist2{idx2(1,0), idx2(2,2)})
            .contract(asDiagonal(MPS.L_tail), idxlist1{idx2(0,0)})
            .contract(MPS.G.A.conjugate(), idxlist2{idx2(3,1),idx2(2,0)})
            .shuffle(array3{0,2,1});
    block = block_enlarged;

}

void class_environment_R::enlarge(const class_MPS &MPS, const Tensor4 &W){
    picture.append(single_picture);
    size++;
    Tensor3 block_enlarged = block.contract(asDiagonal(MPS.L.B), idxlist1{idx2(0,0)})
            .contract(MPS.G.B, idxlist1{idx2(2,2)})
            .contract(W, idxlist2{idx2(1,1), idx2(2,3)})
            .contract(asDiagonal(MPS.L.B), idxlist1{idx2(0,1)})
            .contract(MPS.G.B.conjugate(), idxlist2{idx2(3,2),idx2(2,0)})
            .shuffle(array3{0,2,1});
    block = block_enlarged;
}











