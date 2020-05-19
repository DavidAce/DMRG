//
// Created by david on 2020-05-12.
//

#include "class_env_var.h"
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>

class_env_var::class_env_var(std::string side_, const class_mps_site & MPS, const class_mpo_site &MPO):
    class_env_base(side_,MPS,MPO){set_edge_dims(MPS,MPO);}

void class_env_var::clear() {
    block = Eigen::Tensor<Scalar,4>();
}

void class_env_var::assert_validity() const {
    if(Textra::hasNaN(block, "env2" + side)){
        throw std::runtime_error("enrironment2 " + side +  " at position " + std::to_string(get_position()) + " has NAN's");
    }
}

bool class_env_var::is_real() const {
    return Textra::isReal(block, "env2" + side);
}

bool class_env_var::has_nan() const {
    return Textra::hasNaN(block, "env2" + side);
}


class_env_var class_env_var::enlarge(const class_mps_site & MPS, const class_mpo_site &MPO){
    if(MPS.get_position() != MPO.get_position()) throw std::logic_error("MPS and MPO not at the same position!");
    class_env_var env = *this;
    if(env.sites == 0 and not env.edge_has_been_set){
        env.set_edge_dims(MPS,MPO);
        env.position = MPS.get_position();
        return env;
    }
    env.enlarge(MPS.get_M_bare(),MPO.MPO());
    if (env.side == "L"){
        env.position = MPS.get_position() + 1;
    }else if (env.side == "R"){
        env.position = MPS.get_position() - 1;
    }else{
        throw std::logic_error("Expected environment side L or R, got: " + side);
    }

    return env;
}


void class_env_var::enlarge(const Eigen::Tensor<Scalar,3>  &MPS, const Eigen::Tensor<Scalar,4> &MPO){
    /*!< Contracts a site into the block. */
//    if(sites == 0 and not edge_has_been_set){set_edge_dims(MPS,MPO);}
    Eigen::Tensor<Scalar,4> block_enlarged;
    if (side == "L"){

        /*! # Left environment contraction
         * [      ]--0 0--[LB_diag]--1 1--[  GA   ]--2
         * [      ]                      |
         * [      ]                      0
         * [      ]
         * [      ]                      2
         * [      ]                      |
         * [      ]--2            0--[   M    ]--1
         * [      ]                      |
         * [      ]                      3
         * [ left ]
         * [      ]                      2
         * [      ]                      |
         * [      ]--3            0--[   M    ]--1
         * [      ]                      |
         * [      ]                      3
         * [      ]
         * [      ]                      0
         * [      ]                      |
         * [      ]--1 0--[LB_diag]--1  1--[GA conj ]--2
         */
        if(MPS.dimension(0) != MPO.dimension  (2)) throw std::runtime_error(fmt::format("ENV2 L pos {} dimension mismatch: MPS dim[{}]:{} != MPO   dim[{})]:{}",position.value(),0,MPS.dimension(0) ,2,MPO.dimension  (2)));
        if(MPS.dimension(1) != block.dimension(0)) throw std::runtime_error(fmt::format("ENV2 L pos {} dimension mismatch: MPS dim[{}]:{} != block dim[{})]:{}",position.value(),1,MPS.dimension(1) ,0,block.dimension(0)));
        if(MPO.dimension(0) != block.dimension(2)) throw std::runtime_error(fmt::format("ENV2 L pos {} dimension mismatch: MPO dim[{}]:{} != block dim[{})]:{}",position.value(),0,MPO.dimension(0) ,2,block.dimension(2)));

        sites++;
        block_enlarged =
            block.contract(MPS,                Textra::idx({0},{1}))
                .contract(MPO,                 Textra::idx({1,3},{0,2}))
                .contract(MPO,                 Textra::idx({1,4},{0,2}))
                .contract(MPS.conjugate(),     Textra::idx({0,4},{1,0}))
                .shuffle(Textra::array4{0,3,1,2});
        block = block_enlarged;
    }
    if (side == "R"){
        /*! # Right environment contraction
         *  1--[   GB    ]--2 0--[LB_diag]--1  0--[      ]
         *          |                        [      ]
         *          0                        [      ]
         *                                   [      ]
         *          2                        [      ]
         *          |                        [      ]
         *   0--[   M    ]--1             2--[      ]
         *          |                        [      ]
         *          3                        [      ]
         *                                   [ right]
         *          2                        [      ]
         *          |                        [      ]
         *   0--[   M    ]--1             2--[      ]
         *          |                        [      ]
         *          3                        [      ]
         *                                   [      ]
         *          0                        [      ]
         *          |                        [      ]
         *  1--[ GB conj ]--2 0--[LB_diag]--1  1--[      ]
        */
        if(MPS.dimension(0) != MPO.dimension  (2)) throw std::runtime_error(fmt::format("ENV R pos {} dimension mismatch: MPS dim[{}]:{} != MPO   dim[{})]:{}",position.value(),0,MPS.dimension(0) ,2,MPO.dimension  (2)));
        if(MPS.dimension(2) != block.dimension(0)) throw std::runtime_error(fmt::format("ENV R pos {} dimension mismatch: MPS dim[{}]:{} != block dim[{})]:{}",position.value(),2,MPS.dimension(2) ,0,block.dimension(0)));
        if(MPO.dimension(1) != block.dimension(2)) throw std::runtime_error(fmt::format("ENV R pos {} dimension mismatch: MPO dim[{}]:{} != block dim[{})]:{}",position.value(),1,MPO.dimension(1) ,2,block.dimension(2)));

        sites++;
        block_enlarged =
            block.contract(MPS,            Textra::idx({0},{2}))
                .contract(MPO,             Textra::idx({1,3},{1,2}))
                .contract(MPO,             Textra::idx({1,4},{1,2}))
                .contract(MPS.conjugate(), Textra::idx({0,4},{2,0}))
                .shuffle(Textra::array4{0, 3, 1, 2});
        block = block_enlarged;
    }
}

void class_env_var::set_edge_dims(const class_mps_site & MPS, const class_mpo_site &MPO) {
    if(edge_has_been_set) return;
    if (side == "L") {
        long mpsDim = MPS.get_chiL();
        long mpoDim = MPO.MPO().dimension(0);
        block.resize(Textra::array4{mpsDim,mpsDim, mpoDim, mpoDim});
        block.setZero();
        Eigen::Tensor<Scalar,2> double_edge = MPO.get_MPO_edge_left().contract(MPO.get_MPO_edge_left(),Textra::idx());
        for (long i = 0; i < mpsDim; i++){
            Eigen::array<long, 2> extent2 = {mpoDim,mpoDim};
            Eigen::array<long, 4> offset4 = {i,i,0,0};
            Eigen::array<long, 4> extent4 = {1,1,mpoDim,mpoDim};
            block.slice(offset4, extent4).reshape(extent2) =  double_edge;
        }
    }
    if(side == "R"){
        long mpsDim = MPS.get_chiR();
        long mpoDim = MPO.MPO().dimension(1);
        block.resize(Textra::array4{mpsDim,mpsDim, mpoDim, mpoDim});
        block.setZero();
        Eigen::Tensor<Scalar,2> double_edge = MPO.get_MPO_edge_right().contract(MPO.get_MPO_edge_right(),Textra::idx());
        for (long i = 0; i < mpsDim; i++){
            Eigen::array<long, 2> extent2 = {mpoDim,mpoDim};
            Eigen::array<long, 4> offset4 = {i,i,0,0};
            Eigen::array<long, 4> extent4 = {1,1,mpoDim,mpoDim};
            block.slice(offset4, extent4).reshape(extent2) =  double_edge;
        }
    }
    sites = 0;
    edge_has_been_set = true;
}


