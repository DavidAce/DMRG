//
// Created by david on 2020-05-12.
//

#include "class_env_ene.h"
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>

#include <utility>

class_env_ene::class_env_ene(std::string side_, const class_mps_site &MPS, const class_mpo_site &MPO) : class_env_base(std::move(side_), MPS, MPO) {
    set_edge_dims(MPS, MPO);
}

void class_env_ene::clear() {
    // Do not clear the empty edge
    if(sites > 0) block.resize(0, 0, 0); // = Eigen::Tensor<Scalar,3>();
}

bool class_env_ene::has_block() const { return block.size() != 0; }

void class_env_ene::assert_validity() const {
    if(Textra::hasNaN(block, "env " + side)) {
        throw std::runtime_error(fmt::format("Environment {} at position {} has NAN's",side,get_position()));
    }
}

bool class_env_ene::is_real() const { return Textra::isReal(block, "env " + side); }

bool class_env_ene::has_nan() const { return Textra::hasNaN(block, "env " + side); }

class_env_ene class_env_ene::enlarge(const class_mps_site &MPS, const class_mpo_site &MPO) {
    if(MPS.get_position() != MPO.get_position())
        throw std::logic_error(fmt::format("MPS and MPO have different positions: {} != {}", MPS.get_position(), MPO.get_position()));

    if(not edge_has_been_set) throw std::logic_error("Have to set edge dimensions first!");

    class_env_ene env = *this;

    env.enlarge(MPS.get_M_bare(), MPO.MPO());
    // Update positions assuming this is a finite chain.
    // This needs to be corrected (on the right side) on infinite chains
    if(env.side == "L") {
        env.position = MPS.get_position() + 1;
    } else if(env.side == "R") {
        env.position = MPS.get_position() - 1;
    } else {
        throw std::logic_error("Expected environment side L or R, got: " + side);
    }

    return env;
}

void class_env_ene::enlarge(const Eigen::Tensor<Scalar, 3> &MPS, const Eigen::Tensor<Scalar, 4> &MPO) {
    /*!< Contracts a site into the block. */
    //    if(sites == 0 and not edge_has_been_set){set_edge_dims(MPS,MPO);}
    if(not has_block()) throw std::runtime_error(fmt::format("env_ene pos {} has no block to enlarge", get_position()));
    if(side == "L") {
        /*! # Left environment contraction
         * [      ]--0 0--[LB]--1 1--[  GA    ]--2
         * [      ]                      |
         * [      ]                      0
         * [      ]
         * [      ]                      2
         * [      ]                      |
         * [ left ]--2            0--[   M    ]--1
         * [      ]                      |
         * [      ]                      3
         * [      ]
         * [      ]                      0
         * [      ]                      |
         * [      ]--1 0--[LB]--1  1--[GA conj ]--2
         */

        if(MPS.dimension(0) != MPO.dimension(2))
            throw std::runtime_error(
                fmt::format("ENV L pos {} dimension mismatch: MPS dim[{}]:{} != MPO   dim[{}]:{}", position.value(), 0, MPS.dimension(0), 2, MPO.dimension(2)));
        if(MPS.dimension(1) != block.dimension(0))
            throw std::runtime_error(fmt::format("ENV L pos {} dimension mismatch: MPS dim[{}]:{} != block dim[{}]:{}", position.value(), 1, MPS.dimension(1),
                                                 0, block.dimension(0)));
        if(MPO.dimension(0) != block.dimension(2))
            throw std::runtime_error(fmt::format("ENV L pos {} dimension mismatch: MPO dim[{}]:{} != block dim[{}]:{}", position.value(), 0, MPO.dimension(0),
                                                 2, block.dimension(2)));

        sites++;
        Eigen::Tensor<Scalar, 3> block_enlarged = block.contract(MPS, Textra::idx({0}, {1}))
                                                      .contract(MPO, Textra::idx({1, 2}, {0, 2}))
                                                      .contract(MPS.conjugate(), Textra::idx({0, 3}, {1, 0}))
                                                      .shuffle(Textra::array3{0, 2, 1});
        block = block_enlarged;
    } else if(side == "R") {
        /*! # Right environment contraction
         *  1--[ GB conj ]--2 0--[LB]--1  0--[      ]
         *          |                        [      ]
         *          0                        [      ]
         *                                   [      ]
         *          2                        [      ]
         *          |                        [      ]
         *   0--[   M    ]--1             2--[ right]
         *          |                        [      ]
         *          3                        [      ]
         *                                   [      ]
         *          0                        [      ]
         *          |                        [      ]
         *    1--[  GB   ]--2 0--[LB]--1  1--[      ]
         */

        if(MPS.dimension(0) != MPO.dimension(2))
            throw std::runtime_error(
                fmt::format("ENV R pos {} dimension mismatch: MPS dim[{}]:{} != MPO   dim[{}]:{}", position.value(), 0, MPS.dimension(0), 2, MPO.dimension(2)));
        if(MPS.dimension(2) != block.dimension(0))
            throw std::runtime_error(fmt::format("ENV R pos {} dimension mismatch: MPS dim[{}]:{} != block dim[{}]:{}", position.value(), 2, MPS.dimension(2),
                                                 0, block.dimension(0)));
        if(MPO.dimension(1) != block.dimension(2))
            throw std::runtime_error(fmt::format("ENV R pos {} dimension mismatch: MPO dim[{}]:{} != block dim[{}]:{}", position.value(), 1, MPO.dimension(1),
                                                 2, block.dimension(2)));
        sites++;
        Eigen::Tensor<Scalar, 3> block_enlarged = block.contract(MPS, Textra::idx({0}, {2}))
                                                      .contract(MPO, Textra::idx({1, 2}, {1, 2}))
                                                      .contract(MPS.conjugate(), Textra::idx({0, 3}, {2, 0}))
                                                      .shuffle(Textra::array3{0, 2, 1});
        block = block_enlarged;
    }
}

void class_env_ene::set_edge_dims(const class_mps_site &MPS, const class_mpo_site &MPO) {
    if(edge_has_been_set) return;
    if(side == "L") {
        long mpsDim = MPS.get_chiL();
        long mpoDim = MPO.MPO().dimension(0);
        block.resize(Textra::array3{mpsDim, mpsDim, mpoDim});
        block.setZero();
        for(long i = 0; i < mpsDim; i++) {
            Eigen::array<long, 1> extent1                  = {mpoDim};
            Eigen::array<long, 3> offset3                  = {i, i, 0};
            Eigen::array<long, 3> extent3                  = {1, 1, mpoDim};
            block.slice(offset3, extent3).reshape(extent1) = MPO.get_MPO_edge_left();
        }
    }
    if(side == "R") {
        long mpsDim = MPS.get_chiR();
        long mpoDim = MPO.MPO().dimension(1);
        block.resize(Textra::array3{mpsDim, mpsDim, mpoDim});
        block.setZero();
        for(long i = 0; i < mpsDim; i++) {
            Eigen::array<long, 1> extent1                  = {mpoDim};
            Eigen::array<long, 3> offset3                  = {i, i, 0};
            Eigen::array<long, 3> extent3                  = {1, 1, mpoDim};
            block.slice(offset3, extent3).reshape(extent1) = MPO.get_MPO_edge_right();
        }
    }
    sites             = 0;
    edge_has_been_set = true;
}
