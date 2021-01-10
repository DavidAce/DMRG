//
// Created by david on 2020-05-12.
//

#include "class_env_base.h"
#include <config/debug.h>
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>
#include <math/hash.h>
#include <math/num.h>
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
        unique_id         = other.unique_id;
        unique_id_env     = other.unique_id_env;
        unique_id_mps     = other.unique_id_mps;
        unique_id_mpo     = other.unique_id_mpo;
    }
    assert_block();
    return *this;
}

using Scalar = class_env_base::Scalar;

void class_env_base::assert_block() const {
    if(not block) throw std::runtime_error(fmt::format("env {} {} at pos {} is null", tag, side, get_position()));
}

void class_env_base::build_block(Eigen::Tensor<Scalar, 3> &otherblock, const Eigen::Tensor<Scalar, 3> &mps, const Eigen::Tensor<Scalar, 4> &mpo) {
    /*!< Contracts a site into the block-> */
    // Note that otherblock, mps and mpo should correspond to the same site! I.e. their "get_position()" are all equal.
    // This can't be checked here though, so do that before calling this function.
    unique_id     = std::nullopt;
    unique_id_env = std::nullopt;
    unique_id_mps = std::nullopt;
    unique_id_mpo = std::nullopt;

    if(not block) block = std::make_unique<Eigen::Tensor<Scalar, 3>>();
    if(side == "L") {
        /*! # Left environment block contraction
         *   [       ]--0         [       ]--0 0--[LB]--1 1--[  GA    ]--2
         *   [       ]            [       ]                      |
         *   [       ]            [       ]                      0
         *   [       ]            [       ]
         *   [       ]            [       ]                      2
         *   [       ]            [       ]                      |
         *   [ block ]--2     =   [lblock ]--2            0--[   M    ]--1
         *   [       ]            [       ]                      |
         *   [       ]            [       ]                      3
         *   [       ]            [       ]
         *   [       ]            [       ]                      0
         *   [       ]            [       ]                      |
         *   [       ]--1         [       ]--1 0--[LB]--1  1--[GA conj ]--2
         */

        if(mps.dimension(0) != mpo.dimension(2))
            throw std::runtime_error(fmt::format("env{} {} pos {} dimension mismatch: mps dim[{}]:{} != mpo   dim[{}]:{}", side, tag, position.value(), 0,
                                                 mps.dimension(0), 2, mpo.dimension(2)));
        if(mps.dimension(1) != otherblock.dimension(0))
            throw std::runtime_error(fmt::format("env{} {} pos {} dimension mismatch: mps dim[{}]:{} != left-block dim[{}]:{}", side, tag, position.value(), 1,
                                                 mps.dimension(1), 0, otherblock.dimension(0)));
        if(mpo.dimension(0) != otherblock.dimension(2))
            throw std::runtime_error(fmt::format("env{} {} pos {} dimension mismatch: mpo dim[{}]:{} != left-block dim[{}]:{}", side, tag, position.value(), 0,
                                                 mpo.dimension(0), 2, otherblock.dimension(2)));

        block->resize(mps.dimension(2), mps.dimension(2), mpo.dimension(1));
        block->device(Textra::omp::getDevice()) = otherblock.contract(mps, Textra::idx({0}, {1}))
                                                      .contract(mpo, Textra::idx({1, 2}, {0, 2}))
                                                      .contract(mps.conjugate(), Textra::idx({0, 3}, {1, 0}))
                                                      .shuffle(Textra::array3{0, 2, 1});
    } else if(side == "R") {
        /*! # Right environment contraction
         *   0--[       ]          1--[   GB   ]--2  0--[LB]--1  0--[       ]
         *      [       ]                  |                        [       ]
         *      [       ]                  0                        [       ]
         *      [       ]                                           [       ]
         *      [       ]                  2                        [       ]
         *      [       ]                  |                        [       ]
         *   2--[ block ]     =     0--[   M   ]--1              2--[ rblock]
         *      [       ]                  |                        [       ]
         *      [       ]                  3                        [       ]
         *      [       ]                                           [       ]
         *      [       ]                  0                        [       ]
         *      [       ]                  |                        [       ]
         *   1--[       ]           1--[GB conj]--2  0--[LB]--1  1--[       ]
         */

        if(mps.dimension(0) != mpo.dimension(2))
            throw std::runtime_error(fmt::format("env{} {} pos {} dimension mismatch: mps dim[{}]:{} != mpo dim[{}]:{}", side, tag, position.value(), 0,
                                                 mps.dimension(0), 2, mpo.dimension(2)));
        if(mps.dimension(2) != otherblock.dimension(0))
            throw std::runtime_error(fmt::format("env{} {} pos {} dimension mismatch: mps dim[{}]:{} != right-block dim[{}]:{}", side, tag, position.value(), 2,
                                                 mps.dimension(2), 0, otherblock.dimension(0)));
        if(mpo.dimension(1) != otherblock.dimension(2))
            throw std::runtime_error(fmt::format("env{} {} pos {} dimension mismatch: mpo dim[{}]:{} != right-block dim[{}]:{}", side, tag, position.value(), 1,
                                                 mpo.dimension(1), 2, otherblock.dimension(2)));
        block->resize(mps.dimension(1), mps.dimension(1), mpo.dimension(0));
        block->device(Textra::omp::getDevice()) = otherblock.contract(mps, Textra::idx({0}, {2}))
                                                      .contract(mpo, Textra::idx({1, 2}, {1, 2}))
                                                      .contract(mps.conjugate(), Textra::idx({0, 3}, {2, 0}))
                                                      .shuffle(Textra::array3{0, 2, 1});
    }
}

void class_env_base::enlarge(const Eigen::Tensor<Scalar, 3> &mps, const Eigen::Tensor<Scalar, 4> &mpo) {
    /*!< Contracts a site into the current block-> */

    // NOTE:
    // If side == "L", the mps and mpo should correspond to this env position+1
    // If side == "R", the mps and mpo should correspond to this env position-1
    // There is no way to check the positions here, so checks should be done before calling this function

    assert_block();
    Eigen::Tensor<Scalar, 3> thisblock = *block;
    build_block(thisblock, mps, mpo);
    sites++;
    if(position) {
        if(side == "L")
            position = position.value() + 1;
        else if(position.value() > 0)
            position = position.value() - 1;
    }
}

void class_env_base::clear() {
    // Do not clear the empty edge
    assert_block();
    if(sites > 0) {
        block->resize(0, 0, 0); // = Eigen::Tensor<Scalar,3>();
        tools::log->trace("Ejected env{} pos {}", side, get_position());
        unique_id = std::nullopt;
    }
    unique_id_env = std::nullopt;
    unique_id_mps = std::nullopt;
    unique_id_mpo = std::nullopt;
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

void class_env_base::assert_unique_id(const class_env_base &env, const class_mps_site &mps, const class_mpo_site &mpo) const {
    std::string msg;
    if(env.get_unique_id() != unique_id_env) {
        msg.append(fmt::format("| env({}) {} !=", env.get_position(), env.get_unique_id()));
        if(unique_id_env) msg.append(fmt::format(" {} ", unique_id_env.value()));
    }
    if(mps.get_unique_id() != unique_id_mps) {
        msg.append(fmt::format("| mps({}) {} !=", mps.get_position(), mps.get_unique_id()));
        if(unique_id_mps) msg.append(fmt::format(" {} ", unique_id_mps.value()));
    }
    auto mpo_unique_id = tag == "ene" ? mpo.get_unique_id() : mpo.get_unique_id_sq();
    if(mpo_unique_id != unique_id_mpo) {
        msg.append(fmt::format("| mpo({}) {} !=", mpo.get_position(), mpo_unique_id));
        if(unique_id_mpo) msg.append(fmt::format(" {} ", unique_id_mpo.value()));
    }
    if(not msg.empty()) throw std::runtime_error(fmt::format("Environment {} side {}: unique id mismatch: {}", tag, side, msg));
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
    if(position)
        return position.value();
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
    unique_id         = std::nullopt;
    unique_id_env     = std::nullopt;
    unique_id_mps     = std::nullopt;
    unique_id_mpo     = std::nullopt;
}

std::size_t class_env_base::get_unique_id() const {
    if(unique_id) return unique_id.value();
    unique_id = hash::hash_buffer(get_block().data(), static_cast<size_t>(get_block().size()));
    return unique_id.value();
}

std::optional<std::size_t> class_env_base::get_unique_id_env() const { return unique_id_env; }
std::optional<std::size_t> class_env_base::get_unique_id_mps() const { return unique_id_mps; }
std::optional<std::size_t> class_env_base::get_unique_id_mpo() const { return unique_id_mpo; }
