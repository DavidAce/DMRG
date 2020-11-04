//
// Created by david on 2020-05-12.
//

#include "class_env_base.h"
#include <general/nmspc_tensor_extra.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>
#include <tools/common/fmt.h>
#include <utility>

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_env_base::class_env_base() : block(std::make_unique<Eigen::Tensor<Scalar, 3>>()) { assert_block(); } // default ctor
class_env_base::~class_env_base()                               = default;                                 // default dtor
class_env_base::class_env_base(class_env_base &&other) noexcept = default;                                 // default move ctor
class_env_base &class_env_base::operator=(class_env_base &&other) noexcept = default;                      // default move assign

class_env_base::class_env_base(const class_env_base &other)
    : edge_has_been_set(other.edge_has_been_set), block(std::make_unique<Eigen::Tensor<Scalar, 3>>(*other.block)), sites(other.sites), position(other.position),
      side(other.side), tag(other.tag) {
    assert_block();
}

class_env_base::class_env_base(std::string side_, size_t position_)
    : block(std::make_unique<Eigen::Tensor<Scalar, 3>>()), position(position_), side(std::move(side_)) {
    assert_block();
}
class_env_base::class_env_base(std::string side_, const class_mps_site &MPS, const class_mpo_site &MPO)
    : block(std::make_unique<Eigen::Tensor<Scalar, 3>>()), side(std::move(side_)) {
    if(MPS.get_position() != MPO.get_position())
        throw std::logic_error(fmt::format("MPS and MPO have different positions: {} != {}", MPS.get_position(), MPO.get_position()));
    position = MPS.get_position();
    assert_block();
}

class_env_base &class_env_base::operator=(const class_env_base &other) {
    if(this != &other) {
        edge_has_been_set = other.edge_has_been_set;
        block             = std::make_unique<Eigen::Tensor<Scalar, 3>>(*other.block);
        sites             = other.sites;
        position          = other.position;
        side              = other.side;
        tag               = other.tag;
    }
    assert_block();
    return *this;
}

using Scalar = class_env_base::Scalar;

void class_env_base::assert_block() const {
    if(block == nullptr) throw std::runtime_error(fmt::format("env {} {} at pos {} is null", tag, side, get_position()));
}

void class_env_base::enlarge(const Eigen::Tensor<Scalar, 3> &MPS, const Eigen::Tensor<Scalar, 4> &MPO) {
    /*!< Contracts a site into the block-> */
    //    if(sites == 0 and not edge_has_been_set){set_edge_dims(MPS,MPO);}
    assert_block();
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
            throw std::runtime_error(fmt::format("env {} {} pos {} dimension mismatch: MPS dim[{}]:{} != MPO   dim[{}]:{}", tag, side, position.value(), 0,
                                                 MPS.dimension(0), 2, MPO.dimension(2)));
        if(MPS.dimension(1) != block->dimension(0))
            throw std::runtime_error(fmt::format("env {} {} pos {} dimension mismatch: MPS dim[{}]:{} != block dim[{}]:{}", tag, side, position.value(), 1,
                                                 MPS.dimension(1), 0, block->dimension(0)));
        if(MPO.dimension(0) != block->dimension(2))
            throw std::runtime_error(fmt::format("env {} {} pos {} dimension mismatch: MPO dim[{}]:{} != block dim[{}]:{}", tag, side, position.value(), 0,
                                                 MPO.dimension(0), 2, block->dimension(2)));

        sites++;
        Eigen::Tensor<Scalar, 3> block_enlarged = block->contract(MPS, Textra::idx({0}, {1}))
                                                      .contract(MPO, Textra::idx({1, 2}, {0, 2}))
                                                      .contract(MPS.conjugate(), Textra::idx({0, 3}, {1, 0}))
                                                      .shuffle(Textra::array3{0, 2, 1});
        *block = block_enlarged;
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
            throw std::runtime_error(fmt::format("env {} {} pos {} dimension mismatch: MPS dim[{}]:{} != MPO   dim[{}]:{}", tag, side, position.value(), 0,
                                                 MPS.dimension(0), 2, MPO.dimension(2)));
        if(MPS.dimension(2) != block->dimension(0))
            throw std::runtime_error(fmt::format("env {} {} pos {} dimension mismatch: MPS dim[{}]:{} != block dim[{}]:{}", tag, side, position.value(), 2,
                                                 MPS.dimension(2), 0, block->dimension(0)));
        if(MPO.dimension(1) != block->dimension(2))
            throw std::runtime_error(fmt::format("env {} {} pos {} dimension mismatch: MPO dim[{}]:{} != block dim[{}]:{}", tag, side, position.value(), 1,
                                                 MPO.dimension(1), 2, block->dimension(2)));
        sites++;
        Eigen::Tensor<Scalar, 3> block_enlarged = block->contract(MPS, Textra::idx({0}, {2}))
                                                      .contract(MPO, Textra::idx({1, 2}, {1, 2}))
                                                      .contract(MPS.conjugate(), Textra::idx({0, 3}, {2, 0}))
                                                      .shuffle(Textra::array3{0, 2, 1});
        *block = block_enlarged;
    }
}

void class_env_base::clear() {
    // Do not clear the empty edge
    assert_block();
    if(sites > 0) block->resize(0, 0, 0); // = Eigen::Tensor<Scalar,3>();
}

const Eigen::Tensor<Scalar, 3> &class_env_base::get_block() const {
    assert_block();
    return *block;
}
Eigen::Tensor<Scalar, 3> &class_env_base::get_block() {
    assert_block();
    return *block;
}

bool class_env_base::has_block() const { return block != nullptr and block->size() != 0; }

void class_env_base::assert_validity() const {
    assert_block();
    if(Textra::hasNaN(*block, fmt::format("env {} {}", tag, side))) {
        throw std::runtime_error(fmt::format("Environment {} side {} at position {} has NAN's", tag, side, get_position()));
    }
}

bool class_env_base::is_real() const {
    assert_block();
    return Textra::isReal(*block, fmt::format("env {} {}", tag, side));
}

bool class_env_base::has_nan() const {
    assert_block();
    return Textra::hasNaN(*block, fmt::format("env {} {}", tag, side));
}

size_t class_env_base::get_position() const {
    if(position) return position.value();
    else
        throw std::runtime_error(fmt::format("Position hasn't been set on env side {}", side));
}

size_t class_env_base::get_sites() const { return sites; }

void class_env_base::set_edge_dims(const Eigen::Tensor<Scalar, 3> &MPS, const Eigen::Tensor<Scalar, 4> &MPO, const Eigen::Tensor<Scalar, 1> &edge) {
    if(edge_has_been_set) return;
    assert_block();
    if(side == "L") {
        long mpsDim = MPS.dimension(1);
        long mpoDim = MPO.dimension(0);
        block->resize(Textra::array3{mpsDim, mpsDim, mpoDim});
        block->setZero();
        for(long i = 0; i < mpsDim; i++) {
            Eigen::array<long, 1> extent1                   = {mpoDim};
            Eigen::array<long, 3> offset3                   = {i, i, 0};
            Eigen::array<long, 3> extent3                   = {1, 1, mpoDim};
            block->slice(offset3, extent3).reshape(extent1) = edge;
        }
    }
    if(side == "R") {
        long mpsDim = MPS.dimension(2);
        long mpoDim = MPO.dimension(1);
        block->resize(Textra::array3{mpsDim, mpsDim, mpoDim});
        block->setZero();
        for(long i = 0; i < mpsDim; i++) {
            Eigen::array<long, 1> extent1                   = {mpoDim};
            Eigen::array<long, 3> offset3                   = {i, i, 0};
            Eigen::array<long, 3> extent3                   = {1, 1, mpoDim};
            block->slice(offset3, extent3).reshape(extent1) = edge;
        }
    }
    sites             = 0;
    edge_has_been_set = true;
}
